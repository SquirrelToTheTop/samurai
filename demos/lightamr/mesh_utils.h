#pragma once

#include <array>
#include <cstdint>
#include <cassert>

#define LVLMAX 32

struct Logical_Pos_t {
    uint32_t i, j, k;
};

struct Pos_t {
    double x, y , z;
};

namespace MeshUtils {

    struct Cell{
        Logical_Pos_t coord;
        uint32_t dataIndex;
        uint8_t level;
        bool isCoarse; // True -> coarse

        // used for unit test !
        inline
        friend bool operator==(const Cell &c1, const Cell &c2){
        return c1.dataIndex == c2.dataIndex && c1.level == c2.level && c1.isCoarse == c2.isCoarse
                && c1.coord.i == c2.coord.i && c1.coord.j == c2.coord.j && c1.coord.k == c2.coord.k; 
        }
    };

    // ordres des cellules feuilles pour le calcul des positions logiques pour un octree de raffinement 2
    const uint32_t _IJK_xorder[8] = {0, 1, 0, 1, 0, 1, 0, 1};
    const uint32_t _IJK_yorder[8] = {0, 0, 1, 1, 0, 0, 1, 1};
    const uint32_t _IJK_zorder[8] = {0, 0, 0, 0, 1, 1, 1, 1};

    // ordre des cellules feuilles d'une cellule coarse pour le calcul des centres 3D cartesiens
    const int32_t _XYZ_xorder[8] = {-1, 1,-1, 1,-1, 1,-1, 1};
    const int32_t _XYZ_yorder[8] = {-1,-1, 1, 1,-1,-1, 1, 1};
    const int32_t _XYZ_zorder[8] = {-1,-1,-1,-1, 1, 1, 1, 1};

    /**
     * Precompute at compile time cell sizes for box of lenght [0,1]
     */
    constexpr auto _csize{[]() constexpr{
        std::array<double, LVLMAX> result{};
        result[0] = 1.0;
        for (int i = 1; i < LVLMAX; ++i){
            result[ i ] = result[i - 1] * 0.5;
        }
        return result;
    }()};

    /**
     * return the length of a cell (x axis)
     */
    inline
    constexpr double getCellSize( uint8_t level ) {
        assert( level < LVLMAX );
        return _csize[ level ];
    }

    /**
     * Return cell volume for a cubic box [0,1]x[0,1]x[0,1]
     * 
     */
    inline
    double getCellVolume( uint8_t level ) {
        assert( level < LVLMAX );
        return _csize[level] * _csize[level] * _csize[level]; 
    }

    /* ---------------------------------------------------------------------------------------- */
    /* Function _r2 are hard-coded for a refinement factor of 2 ------------------------------- */
    /* ---------------------------------------------------------------------------------------- */

    /**
     * Get cell center in 3D cartesian coordinates,
     * 
     * Valid for a refinement factor of 2.
     * 
     * Thread-safe
     */
    inline
    Pos_t getCellCenter_r2( const Logical_Pos_t & cijk, uint8_t level, const int dim ) {
        assert( dim == 2 || dim == 3 );
        assert( level < LVLMAX );

        // set to 0 last position if 2D !
        return Pos_t { _csize[ level ] * ( static_cast<double>( cijk.i ) + 0.5 ), 
                       _csize[ level ] * ( static_cast<double>( cijk.j ) + 0.5 ),
                       _csize[ level ] * ( static_cast<double>( cijk.k ) + 0.5 ) * static_cast<double>( dim -2 ) }; 
    }

    inline
    Pos_t getCellCenter_r2( const Cell & c, const int dim ) {
        assert( dim == 2 || dim == 3 );
        assert( c.level < LVLMAX );

        // set to 0 last position if 2D !
        return Pos_t { _csize[ c.level ] * ( static_cast<double>( c.coord.i ) + 0.5 ), 
                       _csize[ c.level ] * ( static_cast<double>( c.coord.j ) + 0.5 ),
                       _csize[ c.level ] * ( static_cast<double>( c.coord.k ) + 0.5 ) * static_cast<double>( dim - 2 ) }; 
    }

    /**
     * Get leaf center for coarse 3d position. No possible check, so be sure to use correct 3D MeshUtils::Cell 
     * center. This function allow to increase speed by avoiding computing coarse 3D MeshUtils::Cell center  
     * 
     * Valid for a refinement factor of 2.
     * 
     * Thread-safe
     * 
     * @param cc        (in): coarse cell center 3D
     * @param ileaf     (in): leaf index [0,3] or [0,7] for refinement factor of 2
     * @param leafLevel (in): level of the leaf
     * @param dim       (in): dimension
     * 
     * @return 3D coordinates of center of the ileaf cell
     */
    inline
    Pos_t getLeafCenter_r2( const Pos_t & cc, int ileaf, uint8_t leafLevel, const int dim ) {
        assert( ileaf >= 0 && ileaf < 8 );
        assert( leafLevel < LVLMAX );
        assert( cc.x >= 0.0 && cc.y >= 0.0 && cc.z >= 0.0 );

        // set to 0 last position if 2D !
        return Pos_t { cc.x + _XYZ_xorder[ileaf] * _csize[ leafLevel ] * 0.5, 
                       cc.y + _XYZ_yorder[ileaf] * _csize[ leafLevel ] * 0.5,
                       cc.z + _XYZ_zorder[ileaf] * _csize[ leafLevel ] * 0.5 * (dim-2) };
    }

    /**
     * Get leaf center for coarse MeshUtils::Cell
     * 
     * Valid for a refinement factor of 2.
     * 
     * Thread-safe
     */
    inline
    Pos_t getLeafCenter_r2( const MeshUtils::Cell & coarse, int ileaf, const int dim ) {
        assert( coarse.level + 1 < LVLMAX );
        assert( ileaf >= 0 && ileaf < 8 );
        assert( coarse.isCoarse );

        Pos_t cc = getCellCenter_r2( coarse, dim );

        // set to 0 last position if 2D !
        return Pos_t { cc.x + _XYZ_xorder[ileaf] * _csize[ coarse.level + 1 ] * 0.5, 
                       cc.y + _XYZ_yorder[ileaf] * _csize[ coarse.level + 1 ] * 0.5,
                       cc.z + _XYZ_zorder[ileaf] * _csize[ coarse.level + 1 ] * 0.5 * (dim-2) }; 
    }

    /**
     * Compute the logical position of the @param ileaf from coarse logical position (at level + 1)
     * 
     * Valid for a refinement factor of 2.
     */
    inline
    Logical_Pos_t getLeafPosition_r2( const Logical_Pos_t & coarse, int ileaf, const int dim ) {
        assert( ileaf >= 0 && ileaf < 8 );

        return Logical_Pos_t { ( coarse.i << 1 ) + _IJK_xorder[ileaf],
                               ( coarse.j << 1 ) + _IJK_yorder[ileaf],
                               ( ( coarse.k << 1 ) + _IJK_zorder[ileaf]) * (dim-2) };
    }

    inline void
    getLeafPosition_r2( const Logical_Pos_t & coarse, int ileaf, const int dim, Logical_Pos_t &ijk ) {
        assert( ileaf >= 0 && ileaf < 8 );
        ijk.i = ( coarse.i << 1 ) + _IJK_xorder[ileaf];
        ijk.j = ( coarse.j << 1 ) + _IJK_yorder[ileaf];
        ijk.k = ( coarse.k << 1 ) + _IJK_zorder[ileaf] * (dim-2);

    }

    /**
     * Compute coarse logical coordinate from leaf logical position for f_d = 2 (refinement factor)
     * for 2D expect leaf to be 2D ... 
    */
    inline Logical_Pos_t getCoarsePosition_r2( const Logical_Pos_t & leaf ){
        Logical_Pos_t cijk;

        cijk.i = static_cast<uint32_t>( std::floor( leaf.i / 2 ) );
        cijk.j = static_cast<uint32_t>( std::floor( leaf.j / 2 ) );
        cijk.k = static_cast<uint32_t>( std::floor( leaf.k / 2 ) );

        return cijk;
    }

    /*
     * Compute in which child the point @param pt belongs. This suppose the point lies
     * inside the cell ! This must validated when calling this function
     *
     * @param cc : coarse cell center
     * @param pt : point
     * @param dim : dimension
     *
     */
    inline
    uint8_t getChildIndex( const Pos_t &cc, const Pos_t &pt, const int dim ){
        
        uint8_t ichild = 0;
        
        if( pt.z >= cc.z && dim == 3 ){ ichild = 4; }

        if( pt.x <= cc.x ){
            if( pt.y >= cc.y ){
                ichild += 2;
            } // else { ichild += 0; }
        }else{
            if( pt.y >= cc.y ){
                ichild += 3;
            }else{
                ichild += 1;
            }
        }

        return ichild;
    }

}
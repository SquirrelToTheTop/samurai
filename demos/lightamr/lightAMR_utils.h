/**
 * Some usefull function designed for lightAMR format
 * 
 * 
 */

#pragma once

#include <vector>
#include <cstdint>
#include <algorithm>

#include "lightAMR.h"
#include "mesh_utils.h"

namespace LightAMR_Utils{

    /**
     * Precompte the global offset per level where children starts.
     * 
     * ie: sum_ncpl[i] -> offset where to find the very first child of the 
     *                    very first coarse cell at level i
     *
     * @param ncpl (in) : lightAMR number of cells per level
     */
    inline std::vector<uint64_t>
    precomputeSumIndex(const std::vector<uint64_t> &ncpl){
        std::vector<uint64_t> sum_ncpl(ncpl.size());
        sum_ncpl[0] = ncpl[0];
        for( size_t ilvl=1; ilvl < ncpl.size(); ++ilvl ){
            sum_ncpl[ilvl] = ncpl[ilvl] + sum_ncpl[ilvl-1];
        }
        return sum_ncpl;
    } 

    /**
     * Compute level of cell based on its index in the lightAMR desc
     * 
     * @param ncplSumId : vector of summed cell index
     * @param index     : cell index in the tree/mask array
     */
    inline uint8_t
    findLevelFromIndex(const std::vector<uint64_t> &ncplSumId, uint32_t index){
        for(size_t ilvl=0; ilvl<ncplSumId.size(); ++ilvl){
            if( ncplSumId[ilvl] > index ) 
                return static_cast<uint8_t>( ilvl );
        }
        return static_cast<uint8_t>( ncplSumId.size() );
    }

    /*
    * Find parent cell index in "local", ie: level - 1 [0-n]. Retourne la position de la cellule coarse
    * vis-a-vis des autres cellules coarses, ie: la i-th coarse du niveau - 1
    */
    inline uint32_t
    findParentCellLocIndex( const std::vector<uint64_t> &ncplSumId, const std::vector<uint64_t> &nCCpl, 
                            uint32_t index, int nchildren){
        uint8_t clvl = findLevelFromIndex(ncplSumId, index);

        if( clvl <= 1 )
            return 0;

        uint32_t loc_index = index - ncplSumId[clvl-1]; // child local pos
        uint32_t loc_pindex = (loc_index / nchildren) + 1; // division entiere

        return loc_pindex;
    }

    /*
     * For some optimization,
     * do not use this function directly unless you know very exactly what you want
     */
    inline uint32_t
    _findParentCellLocIndex_nocheck( const std::vector<uint64_t> &ncplSumId, const std::vector<uint64_t> &nCCpl, 
                                    uint32_t index, uint8_t clvl, int nchildren){
        
        if( clvl <= 1 ) return 0;

        uint32_t loc_index = index - ncplSumId[clvl-1]; // child local pos
        uint32_t loc_pindex = ( loc_index / nchildren ) + 1; // division entiere

        return loc_pindex;
    }

    /**
     * Compute the number of coarse cell per level
     * 
     * @param tree (in): lightAMR grid refinement array
     * @param ncpl (in): lightAMR number of cells per level
     */
    inline std::vector<uint64_t>
    precomputeCoarsePerLevel(const std::vector<uint8_t> &tree, const std::vector<uint64_t> &ncpl){

        std::vector<uint64_t> nccpl(ncpl.size(), 0);

        uint32_t id = 0;
        for(size_t ilvl=0; ilvl<ncpl.size(); ++ilvl){
            for(uint32_t icell=0; icell<ncpl[ilvl]; ++icell){
                nccpl[ilvl] += tree[id];
                id ++;
            }
        }

        return nccpl;
    }

    /**
     * Compute the max number of coarse in order to reserve buffer
     */
    inline size_t
    getMaxCoarseInLevel(const std::vector<uint8_t> &tree, const std::vector<uint64_t> &ncpl){

        size_t maxc = 0;
        uint32_t id = 0;
        for(size_t ilvl=0; ilvl<ncpl.size(); ++ilvl){

            size_t maxc_loc=0;
            for(uint32_t icell=0; icell<ncpl[ilvl]; ++icell){
                 maxc_loc += tree[id];
                id ++;
            }
            maxc = std::max( maxc, maxc_loc );
        }

        return maxc;
    }

    /**
     * Compute the indexes of coarse cell per level
     */
    inline std::vector<std::vector<uint64_t>>
    precomputeCoarseIndexPerLevel(const std::vector<uint8_t> &tree, const std::vector<uint64_t> &ncpl){

        std::vector<std::vector<uint64_t>> nccplidx(ncpl.size());

        uint32_t id = 0;
        for(size_t ilvl=0; ilvl<ncpl.size(); ++ilvl){
            for(uint32_t icell=0; icell<ncpl[ilvl]; ++icell){
                if( tree[id] == 1 ) nccplidx[ilvl].emplace_back( id );
                id ++;
            }
        }

        return nccplidx;
    }

    /*
     * Compute the logical coordinate (i,j,k) based on a lightAMR cell global index.
     *
     * Params:
     *         nCCpl     (in): number of coarse cell per level
     *         nCCplIdx  (in): coarse cell global index per level
     *         ncplsum   (in): precomputed global index where level starts in lightAMR
     *         clvl      (in): level of global index ( if unknown use `findLevelFromIndex(...)` )
     *         index     (in): global index
     *         nchildren (in): number of children per coarse
     *          
     */
    // inline Logical_Pos_t 
    // getLogicalIJKFromLightAMRIndex(const std::vector<uint64_t> &nCCpl, const std::vector<std::vector<uint64_t>> & nCCplIdx, 
    //                                const std::vector<uint64_t> &ncplsum, int clvl, uint64_t index, int nchildren) {

    //     assert( nchildren == 4 || nchildren == 8 ); // if not power of 2 then modulo cannot be replaced
        
    //     Logical_Pos_t ijk = {0, 0, 0}; // ijk coordinate of current cell

    //     uint64_t plocidxpl[clvl+1]; // -> given in [1:n] based
    //     uint64_t pglobidxpl[clvl+1];

    //     if( clvl >= 1 ){

    //         plocidxpl[0] = 1;
    //         plocidxpl[clvl] = 0;

    //         pglobidxpl[0] = 0;
    //         pglobidxpl[clvl] = index;
    //         for( uint8_t lvl=clvl-1; lvl>0; --lvl){
    //             // plocidxpl[lvl] = findParentCellLocIndex(ncplsum, nCCpl, pglobidxpl[lvl+1], nchildren);
    //             plocidxpl[lvl] = _findParentCellLocIndex_nocheck(ncplsum, nCCpl, pglobidxpl[lvl+1], lvl+1, nchildren);               
    //             pglobidxpl[lvl] = nCCplIdx[lvl][plocidxpl[lvl]-1]; 

    //             // size_t flocidx = (pglobidxpl[lvl] - ncplsum[lvl-1]) % nchildren;
    //         }
            
    //         size_t flocidx = 0;
    //         for(uint8_t lvl=1; lvl<=clvl; ++lvl){
    //             flocidx = (pglobidxpl[lvl] - ncplsum[lvl]) % nchildren;
    //             ijk.i = ( ijk.i << 1 ) + IJKLOOKUP_ZCURVE[flocidx].i;
    //             ijk.j = ( ijk.j << 1 ) + IJKLOOKUP_ZCURVE[flocidx].j;
    //             ijk.k = ( ijk.k << 1 ) + IJKLOOKUP_ZCURVE[flocidx].k;
    //         }

    //     }

    //     return ijk;
    // }

    /*
     * Compute first child index for coarse cells
     */
    static std::vector<size_t>
    computeFirstChildIndex( const std::vector<uint8_t> &tree, const std::vector<uint64_t> &ncpl, int nchildren ) {
        
        auto ncplsum = LightAMR_Utils::precomputeSumIndex( ncpl );

        std::vector<size_t> firstChildIndex( tree.size(), 0 );

        size_t offset = 0;
        for(size_t ilvl=0; ilvl<ncpl.size(); ++ilvl){

            size_t bufferIndirectArrayFirstChild = 0;
            for(size_t icell=0; icell<ncpl[ ilvl ]; ++icell) {
                
                if( tree[ offset ] == 1 ){
                    uint32_t idFirstChild = bufferIndirectArrayFirstChild * nchildren;
                    idFirstChild += ncplsum[ ilvl ];
                    bufferIndirectArrayFirstChild ++;

                    firstChildIndex[ offset ] = idFirstChild;
                }

                offset ++;
            }

        }

        return firstChildIndex;

    }

    /**
     * Find the index in the lightAMR to which the point belongs to, entry point for point sampler. The fact that
     * the point must belong to the inital grid (level 0) must be checked before. Actually if the point does not
     * lies inside the the return index would be undefined behavior (ie: left most grid down to pt_lvl)
     * 
     * @param pt     (in) : 3d cartesian coordinate of point
     * @param pt_lvl (in) : level max for the point 
     * @param lamr   (in) : lightAMR data structure
     * @param fci    (in) : vector of first child index
    */
    static 
    size_t getLightAMRGlobalIndex( const Pos_t &pt, int pt_lvl, const LightAMR &lamr, const std::vector<size_t> &fci ){

        size_t the_index;
         
        int ilvl, ichild;
        int order[ 3 ] = {0, 0, 0};

        double c_cs;
        Pos_t c_cc, tmp;

        // std::cout <<  "    > [findPoint] # ", ilvl, " cell {", pt(1), ", ", pt(2), ", ", pt(3)

        // level 0 cell in lightAMR
        const auto &tree = lamr.getConstReferenceAMRtree();
        const auto &mask = lamr.getConstReferenceMask();
        const auto &ncpl = lamr.getConstReferenceNCellsPerLevel();

        assert( pt_lvl < ncpl.size() );

        // at least for now 
        assert( ncpl[0] == 1 );
        
        c_cs = 0.5;
        c_cc = {0.5, 0.5, 0.5};

        if( lamr.getDim() == 3 )
            c_cc.z = 0.0;

        the_index = 0;
        for( int ilvl=0; ilvl < pt_lvl; ++ilvl ){
            // std::cout << "    > Level # ", ilvl, " cell {", ijk(1), ", ", ijk(2), ", ", ijk(3), "}, csize : ", c_cs

            // find child index
            pt.x > c_cc.x ? order[0] = 1 : order[0] = 0; 
            pt.y > c_cc.y ? order[1] = 1 : order[1] = 0;

            if( lamr.getDim() == 3 ){
                pt.z > c_cc.z ? order[2] = 1 : order[2] = 0;
            }

            // update child info
            c_cs = c_cs * 0.5;

            c_cc.x = c_cc.x + ( 2 * order[0] - 1 ) * c_cs;
            c_cc.y = c_cc.y + ( 2 * order[1] - 1 ) * c_cs;

            if( lamr.getDim() == 3 ){
                c_cc.z = c_cc.z + ( 2 * order[2] - 1 ) * c_cs;
            }

            // child index
            ichild = order[0] + order[ 1 ] * 2;

            if( lamr.getDim() == 3 ) ichild += order[ 2 ] * 4;

            // "        > Child : ", ichild, " ijk (", ijk(1), ", ", ijk(2), ", ", ijk(3), ")"
            // "        > Child size : ", c_cs, " center (", c_cc(1), ", ", c_cc(2), ", ", c_cc(3), ")"

            assert( the_index < fci.size() );
            the_index = fci[ the_index ] + ichild;
            
            //  leaf cell
            assert( the_index < tree.size() );
            if( tree[ the_index ] == 0 ){
                return the_index;
            }

        }

        assert( the_index < mask.size() );
        if( mask[ the_index ] == 1 ){
            std::cerr << "[findPoint] Ownership mask panic !" << std::endl;
        }

        return the_index;
    }

}

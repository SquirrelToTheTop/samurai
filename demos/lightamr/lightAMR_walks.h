/**
 * Some functions that walks through a lightAMR data structure using different
 * kind of algorithms. Some of them are GPU compatible.
 * 
 * WA_uc_* : use the uncompressed lightAMR stream
 * WA_c_*  : use the compressed lightAMR stream (CPS52) directly * 
 * 
 * WA_*_ijk_based : logical_pos_t descent (i,j,k) according to grid refinement
 * WA_*_ra_based  : random access "descent" (there is no real descend in the tree, just loops over cells)
 * 
 * WA_*_ijk_* : are NOT GPU compatible
 * WA_uc_ra_based : is GPU compatible (do more computation than ijk_based)
 *  
 * Written by L. Strafella 
 */

#pragma once

#include <vector>
#include <cstdint>
#include <algorithm>

#include "lightAMR.h"
#include "mesh_utils.h"
#include "lightAMR_utils.h"

namespace LightAMR_Utils{

    /* ---------------------------------------------------------------------------------------------------------- */
    /* IJK based lightAMR walks let's call it "monkey walks" ---------------------------------------------------- */
    /* Use decompressed lightAMR stream, move throught AMR by computing (i,j,k) coordinates and checking tree --- */
    /* ---------------------------------------------------------------------------------------------------------- */
    static std::vector<MeshUtils::Cell>
    WA_uc_ijk_based( const LightAMR &lamr, uint8_t maxLevel ) {

        std::vector<MeshUtils::Cell> vmeshDom;
        vmeshDom.reserve( lamr.getNumberOfCells() * 0.9 );

        int nchildren = 1 << lamr.getDim();

        assert( lamr.getDim() == 3 || lamr.getDim() == 2 );
        assert( lamr.getRefinementFactor() == 2 );
        assert( nchildren == 8  || nchildren == 4 );

        std::vector<Logical_Pos_t> coarse, next_coarse;

        const auto &tree = lamr.getConstReferenceAMRtree();
        const auto &mask = lamr.getConstReferenceMask();
        const auto &cijk = lamr.getConstReferenceCoarseIJK();
        const auto &ncpl = lamr.getConstReferenceNCellsPerLevel();

        uint8_t loc_maxLevel = std::min( maxLevel, static_cast<uint8_t>( ncpl.size() ) );

        std::vector<uint32_t> bufferIndirectArrayFirstChild( ncpl.size(), 0);
        std::vector<uint32_t> sumCellperLevel( ncpl.size(), 0 );

        for(size_t ilvl=0; ilvl<ncpl.size(); ++ilvl){
            for(size_t jlvl=0; jlvl<=ilvl; ++jlvl)
                sumCellperLevel[ ilvl ] += ncpl[ jlvl ];
        }

        assert( !tree.empty() );
        assert( !mask.empty() );
        assert( !cijk.empty() );
        assert( !ncpl.empty() );

        assert( loc_maxLevel > 0 && loc_maxLevel <= ncpl.size() );

        // build coarse array @ level 0
        for( size_t ic=0; ic<cijk.size(); ic+=lamr.getDim() ){
            Logical_Pos_t cij;
            cij.i = cijk[ic];
            cij.j = cijk[ic+1];
            lamr.getDim() == 3 ? cij.k = cijk[ic+2] : cij.k = 0;

            coarse.push_back( cij );
        }

        assert( coarse.size() == cijk.size()/lamr.getDim() );
        assert( coarse.size() == ncpl[0] );

        for(uint8_t ilvl=0; ilvl<loc_maxLevel; ++ilvl){

            for(size_t icoarse=0; icoarse<coarse.size(); ++icoarse) {

                uint32_t idFirstChild = bufferIndirectArrayFirstChild[ilvl] * nchildren;
                idFirstChild += sumCellperLevel[ilvl];
                bufferIndirectArrayFirstChild[ilvl]++;
                    
                for( uint8_t ileaf=0; ileaf < nchildren; ++ileaf ){

                    Logical_Pos_t lij = MeshUtils::getLeafPosition_r2( coarse[icoarse], ileaf, lamr.getDim() );

                    // a leaf that belong to the current domain
                    if( tree[idFirstChild+ileaf] == 0 && mask[idFirstChild+ileaf] == 0 ){

                        // set use coarse IJK as key and position in refAMRTree and first child as value
                        MeshUtils::Cell a = { lij, idFirstChild+ileaf, static_cast<uint8_t>(ilvl+1), false };
                        vmeshDom.push_back( a );

                    }else{

                        // a coarse node ( which does not necessarily belong to the current domain)
                        // even if it does not, we need to go inside, child node might belong to current domain
                        if( tree[ idFirstChild+ileaf ] == 1 ) {

                            if (ilvl == (loc_maxLevel - 1) && mask[ idFirstChild+ileaf ] == 0 ) {
                                MeshUtils::Cell a = { lij, idFirstChild + ileaf, static_cast<uint8_t>(ilvl+1), false };
                                vmeshDom.push_back( a );

                            }else {
                                next_coarse.push_back( lij );
                            }
                        }

                    }

                }

            }

            coarse.clear();
            coarse = next_coarse;
            next_coarse.clear();

        }
        
        vmeshDom.shrink_to_fit();

        // expect RVO
        return vmeshDom;
    }

    // structure for reducing the number of args for recursive functions
    struct _walk_recursive_struct {
        std::vector<uint64_t> _indirectArray;
        std::vector<uint64_t> _ncplsum;
        const uint8_t *tree, *mask;
        int _dim;
        int _nchildren;
        uint8_t _maxLevel;
        bool _found;
    };

    /* ---------------------------------------------------------------------------------------------------------- */
    /* recursive IJK based lightAMR walks let's call it "black squirrel walks" ---------------------------------- */
    /* Use of recursivity with uncompressed tree & mask stream                                                    */
    /* ---------------------------------------------------------------------------------------------------------- */
    static void
    recursive_advanced_uncomp( std::vector<MeshUtils::Cell> &vdom, Logical_Pos_t &ijk, uint32_t id_amr, uint8_t lvl,
                               _walk_recursive_struct & wrs ){

        // reduce tree to max level asked by user
        if( lvl >= wrs._maxLevel ){
            return;
        }

        if( wrs.tree[ id_amr] == 1 ){

            assert( lvl < wrs._indirectArray.size() );
            uint32_t indexFirstChild = wrs._indirectArray[ lvl ] * wrs._nchildren + wrs._ncplsum[ lvl ];
            wrs._indirectArray[ lvl ] ++;

            // uint64_t indexFirstChild = wrs._indirectArray[ id_amr ];

            // check if it a T-node
            int8_t sumRafAndMaskLeaves = 0;
            // int8_t sumMaskLeaves = 0;
            for(int ichild=0; ichild<wrs._nchildren; ++ichild ){
                sumRafAndMaskLeaves += wrs.tree[ indexFirstChild+ichild ] + wrs.mask[ indexFirstChild+ichild ];
                // sumMaskLeaves += mask[ indexFirstChild+ichild ];             
            }

            // if ( sumRafLeaves == 0  && sumMaskLeaves == 0 ){
            if( sumRafAndMaskLeaves == 0 ){
                    // set use coarse IJK as key and position in refAMRTree of first child as value
                    MeshUtils::Cell a = { ijk, indexFirstChild, lvl, true };
                    vdom.push_back( std::move( a ) );
            }else{

                Logical_Pos_t l_ijk[wrs._nchildren];
                for( int ichild=0; ichild<wrs._nchildren; ++ichild ){
                    l_ijk[ ichild ] = MeshUtils::getLeafPosition_r2( ijk, ichild, wrs._dim );
                }

                for( int ichild=0; ichild<wrs._nchildren; ++ichild ){
                    recursive_advanced_uncomp(vdom, l_ijk[ichild], indexFirstChild+ichild, lvl+1, wrs );
                }

            }

        }else{

            if( wrs.mask[ id_amr] == 0 ){
                // add it because it is a leaf that belongs to current domain
                MeshUtils::Cell a = { ijk, id_amr, lvl, false };
                vdom.push_back( std::move(a) );
            }

        }

    }

    static std::vector<MeshUtils::Cell>
    WA_uc_recursive( const LightAMR &lamr, uint8_t maxLevel ){

        const auto &ncpl = lamr.getConstReferenceNCellsPerLevel();
        assert( ! ncpl.empty() );

        std::vector<uint64_t> ncplsum = LightAMR_Utils::precomputeSumIndex(ncpl);
        
        struct _walk_recursive_struct _wrs;
        _wrs._indirectArray.resize( ncpl.size(), 0 );
        _wrs._ncplsum = LightAMR_Utils::precomputeSumIndex( ncpl );
        _wrs._maxLevel = maxLevel;
        _wrs._dim = lamr.getDim();
        _wrs._nchildren = 1 << lamr.getDim();
        _wrs.tree = lamr.getConstReferenceAMRtree().data();
        _wrs.mask = lamr.getConstReferenceMask().data();

        // _wrs._indirectArray = LightAMR_Utils::computeFirstChildIndex( lamr.getConstReferenceAMRtree(), ncpl, _wrs._nchildren );
        
        Logical_Pos_t ijk = { 0, 0, 0 };

        std::vector<MeshUtils::Cell> vmeshDom;
        recursive_advanced_uncomp( vmeshDom, ijk, 0, 0, _wrs );

        // expect RVO
        return vmeshDom;

    }

}


    // > [Master] Timers for Batch #0
    //                        Name            Total Elapsed (s)               Fraction (%) 
    //               WA_c_ijk_based                     39.8146                      83.145
    //               WA_c_recursive                      1.5599                      3.2576
    //              WA_uc_ijk_based                     0.43622                     0.91096
    //   --->       WA_uc_ra_based                      5.4326                      11.345
    //              WA_uc_recursive                     0.38328                     0.80041
    //                     read_amr                     0.21241                     0.44358
    //                   read_hydro                   6.753e-05                  0.00014102
    //            uncompress_amr_th                    0.046852                    0.097841
    //     ------------------------    ------------------------    ------------------------
    //                        Total                      47.886                         100

    // Adding interface of findParentCellLocIndex(...) taking level

    // > [Master] Timers for Batch #0
	//                        Name            Total Elapsed (s)               Fraction (%) 
	//               WA_c_ijk_based                       40.15                      84.448
	//               WA_c_recursive                      1.4901                      3.1342
	//              WA_uc_ijk_based                     0.43284                     0.91039
	//   --->        WA_uc_ra_based                      4.8427                      10.186
	//              WA_uc_recursive                     0.37879                     0.79671
	//                     read_amr                     0.20626                     0.43384
	//                   read_hydro                   6.389e-05                  0.00013438
	//            uncompress_amr_th                    0.043283                    0.091039
	//     ------------------------    ------------------------    ------------------------
	//                        Total                      47.544                         100

    // replacing modulo operation in getLogicalIJKFromLightAMRIndex by & if power of two

    //  > [Master] Timers for Batch #0
    //                            Name            Total Elapsed (s)               Fraction (%) 
    //                   WA_c_ijk_based                     40.0208                       89.34
    //                   WA_c_recursive                      1.5107                      3.3725
    //                  WA_uc_ijk_based                     0.43061                     0.96126
    //     --->          WA_uc_ra_based                      2.2001                      4.9114
    //                  WA_uc_recursive                     0.38011                     0.84854
    //                         read_amr                     0.20452                     0.45655
    //                       read_hydro                   6.535e-05                  0.00014588
    //                uncompress_amr_th                    0.048928                     0.10922
    //         ------------------------    ------------------------    ------------------------
    //                            Total                      44.796                         100

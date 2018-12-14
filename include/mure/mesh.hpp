#pragma once

#include <iostream>
#include <algorithm>

#include <xtensor/xfixed.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xio.hpp>

#include "box.hpp"
#include "cell_list.hpp"
#include "cell_array.hpp"
#include "static_algorithm.hpp"

namespace mure
{
    template<class MRConfig>
    class Mesh
    {
    public:

        static constexpr auto dim = MRConfig::dim;
        static constexpr auto max_refinement_level = MRConfig::max_refinement_level;
        static constexpr auto ghost_width = std::max(std::max(2*(int)MRConfig::graduation_width - 1,
                                                     (int)MRConfig::max_stencil_width),
                                                     (int)MRConfig::default_s_for_prediction);
        using index_t = typename MRConfig::index_t;
        using coord_index_t = typename MRConfig::coord_index_t;
        using point_t = typename Box<int, dim>::point_t;
        using interval_t = typename MRConfig::interval_t;

        Mesh(Box<double, dim> b, size_t init_level)
        {
            point_t start = b.min_corner()*std::pow(2, init_level);
            point_t end = b.max_corner()*std::pow(2, init_level);

            CellList<MRConfig> dcl;
            dcl[init_level].extend(xt::view(start, xt::drop(0)), xt::view(end, xt::drop(0)));
            dcl[init_level].fill({start[0], end[0]});

            m_cells = {dcl};
            update_ghost_nodes();
        }
    
        void to_stream(std::ostream &os) const
        {
            os << "Cells\n";
            m_cells.to_stream(os);
            os << "\nCells and ghosts\n";
            m_cells_and_ghosts.to_stream(os);
        }
    private:
        void update_ghost_nodes();
        void add_ng_ghosts_and_get_nb_leaves(CellList<MRConfig>& cell_list);
        void add_ghosts_for_level_m1(CellList<MRConfig>& cell_list);

        CellArray<MRConfig> m_cells;
        CellArray<MRConfig> m_cells_and_ghosts;
        index_t m_nb_local_cells;
        index_t m_nb_local_cells_and_ghosts;
    };


    template<class MRConfig>
    void Mesh<MRConfig>::update_ghost_nodes() {
        // FIXME: m_nb_local_cells is set in add_ng_ghosts_and_get_nb_leaves
        //        Then, do we need this line??
        // m_nb_local_cells = _leaves.nb_cells();

        // +/- w ghosts in level + 0 and 1, computation of _nb_local_leaf_cells
        CellList<MRConfig> cell_list;
        add_ng_ghosts_and_get_nb_leaves(cell_list);

        // optionnal +/- w nodes in level - 1
        if (MRConfig::need_pred_from_proj)
            add_ghosts_for_level_m1(cell_list);

        // // get external (remote or replicated) cells needed to compute the ghosts.
        // // Updates also _leaves_to_recv and _leaves_to_send (with leaf indices at this stage)
        // std::vector<std::vector<Cell<Carac>>> needed_nodes( mpi->size() );
        // _add_needed_or_replicated_cells( dcl, needed_nodes );

        // compaction
        m_cells_and_ghosts = {cell_list};

        // update of x0_indices, _leaf_to_ghost_indices, _nb_local_ghost_and_leaf_cells
        // update_x0_and_nb_ghosts(m_nb_local_cells_and_ghoss, m_cells_and_ghosts, m_cells);

        // // update of _leaves_to_recv and _leaves_to_send (using needed_nodes)
        // _update_leaves_to_send_and_recv( needed_nodes );

        // // update of _extended_leavesv (union of leaves + remote leaves + symmetry leaves, with correct ghost indices)
        // _update_of_extended_leaves( needed_nodes );
    }

    template<class MRConfig>
    void Mesh<MRConfig>::add_ng_ghosts_and_get_nb_leaves(CellList<MRConfig>& cell_list)
    {
        // +/- w nodes in the current level
        m_nb_local_cells = 0;
        for(int level = 0; level <= max_refinement_level; ++level)
        {
            const LevelCellArray<MRConfig> &level_cell_array = m_cells[level];
            if (level_cell_array.empty() == false)
            {
                LevelCellList<MRConfig> &level_cell_list = level_cell_array[level];
                level_cell_list.extend(level_cell_array.min_corner_yz() - ghost_width,
                                       level_cell_array.max_corner_yz() + ghost_width);
                level_cell_array.for_each_interval_in_x([&](xt::xtensor_fixed<coord_index_t, xt::xshape<dim-1>> const& index_yz,
                                                            interval_t const& interval)
                {
                    static_nested_loop<dim-1, -ghost_width, ghost_width+1>([&](auto stencil)
                    {
                        level_cell_list[index_yz + stencil].add_interval({interval.start - ghost_width,
                                                            interval.end + ghost_width}
                                                     );
                    });
                    m_nb_local_cells += interval.size();
                } );
            }
        }
    }

    template<class MRConfig>
    void Mesh<MRConfig>::add_ghosts_for_level_m1(CellList<MRConfig>& cell_list)
    {
        for(int level = 1; level <= max_refinement_level; ++level)
        {
            const LevelCellArray<MRConfig> &level_cell_array = m_cells[level];
            if (level_cell_array.empty() == false)
            {
                LevelCellList<MRConfig> &level_cell_list = level_cell_array[level - 1];
                constexpr coord_index_t s = MRConfig::default_s_for_prediction;
                
                level_cell_list.extend((level_cell_array.min_corner_yz() >> 1) - s,
                                       (level_cell_array.max_corner_yz() >> 1) + s);
                
                level_cell_array.for_each_interval_in_x([&](xt::xtensor_fixed<coord_index_t, xt::xshape<dim-1>> const& index_yz,
                                                            interval_t const& interval)
                {
                    static_nested_loop<dim-1, -s, s+1>([&](auto stencil)
                    {
                        level_cell_list[(index_yz >> 1) + stencil].add_interval({(interval.start >> 1) - s,
                                                                                 (interval.end >> 1) + s});
                    });
                });
            }
        }
    }

    // void Mesh<MRConfig>::update_x0_and_nb_ghosts(index_t& target_nb_local_cells_and_ghosts,
    //                                              CellList<MRConfig>& target_cells_and_ghosts,
    //                                              CellList<MRConfig>& target_cells)
    // {
    //     // get x0_index for each ghost_and_leaf node
    //     index_t ghost_index = 0;
    //     for(int level = 0; level <= (int)Carac::max_refinement_level; ++level )
    //     {
    //         target_cells_and_ghosts[ level ].for_each_interval_in_x([&](const interval_t& interval)
    //         {
    //             // const_cast<interval_t&>(interval).index = ghost_index - interval.start;
    //             interval.index = ghost_index - interval.start;
    //             ghost_index += interval.size();
    //         } );
    //     }
    //     target_nb_local_cells_and_ghosts = ghost_index;

    //     // get x0_index for each leaf node (ranges are not the same than in target)
    //     index_t cpt_leaf = 0;
    //     m_cell_to_ghost_indices.resize(nb_local_cells());
    //     for( int level = 0; level <= (int)Carac::max_refinement_level; ++level )
    //     {
    //         const LevelData<Carac> &ld_ghost = target_ghost_and_leaves[ level ];
    //         const LevelData<Carac> &ld_plain = target_leaves[ level ];
    //         for_each_x_range( ld_ghost && ld_plain, [&]( const Array<TC,dim-1> &yz, TC beg_x, TC end_x, Array<const Range<TC,Index> *,2> ranges ) {
    //             const_cast<Range<TC,Index> *>( ranges[ 1 ] )->x0_index = ranges[ 0 ]->x0_index;
    //             for( int pos_x = beg_x; pos_x < end_x; ++pos_x )
    //                 _leaf_to_ghost_indices[ cpt_leaf++ ] = ranges[ 0 ]->x0_index + pos_x;
    //         } );
    //     }
    // }
    template<class MRConfig>
    std::ostream& operator<<(std::ostream& out, const Mesh<MRConfig>& mesh)
    {
        mesh.to_stream(out);
        return out;
    }
}
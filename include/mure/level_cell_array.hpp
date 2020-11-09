#pragma once

#include <array>
#include <iterator>
#include <vector>

#include "mure/algorithm.hpp"
#include "mure/box.hpp"
#include "mure/interval.hpp"
#include "mure/level_cell_list.hpp"

namespace mure
{

    template<class LCA, bool is_const>
    class LevelCellArray_iterator;

    template<class iterator>
    class LevelCellArray_reverse_iterator: public std::reverse_iterator<iterator>
    {
    public:
        using base_type = std::reverse_iterator<iterator>;
        using coord_type = typename iterator::coord_type;

        LevelCellArray_reverse_iterator(iterator&& it)
        : base_type(std::move(it))
        {}

        const coord_type& index() const
        {
            iterator it = this->base();
            return (--it).index();
        }
    };

    ///////////////////////////////
    // LevelCellArray definition //
    ///////////////////////////////
    template<std::size_t Dim, class TInterval = Interval<int>>
    class LevelCellArray
    {
    public:
        constexpr static auto dim = Dim;
        using interval_t = TInterval;
        using index_t = typename interval_t::index_t;
        using coord_index_t = typename interval_t::value_t;

        using iterator = LevelCellArray_iterator<LevelCellArray<Dim, TInterval>, false>;
        using reverse_iterator = LevelCellArray_reverse_iterator<iterator>;
        using const_iterator = LevelCellArray_iterator<const LevelCellArray<Dim, TInterval>, true>;
        using const_reverse_iterator = LevelCellArray_reverse_iterator<const_iterator>;

        using index_type = std::array<coord_index_t, dim>;

        LevelCellArray() = default;

        LevelCellArray(const LevelCellArray&) = default;
        LevelCellArray& operator=(const LevelCellArray&) = default;

        LevelCellArray(LevelCellArray&&) = default;
        LevelCellArray& operator=(LevelCellArray&&) = default;

        LevelCellArray(const LevelCellList<Dim, TInterval> &lcl);
        LevelCellArray(std::size_t level, Box<coord_index_t, dim> box);
        LevelCellArray(std::size_t level);

        iterator begin();
        iterator end();

        const_iterator cbegin() const;
        const_iterator cend() const;

        const_reverse_iterator rcend() const;
        const_reverse_iterator rcbegin() const;

        /// Display to the given stream
        void to_stream(std::ostream &out) const;

        //// checks whether the container is empty
        bool empty() const;

        //// Gives the number of intervals in each dimension
        auto shape() const;

        //// Gives the number of cells
        std::size_t nb_cells() const;

        const std::vector<interval_t>& operator[](std::size_t d) const;
        std::vector<interval_t>& operator[](std::size_t d);

        const std::vector<std::size_t>& offsets(std::size_t d) const;
        std::vector<std::size_t>& offsets(std::size_t d);

        std::size_t level() const;

    private:
        /// Recursive construction from a level cell list along dimension > 0
        template<typename TGrid, std::size_t N>
        void initFromLevelCellList(TGrid const &grid,
                                   std::array<coord_index_t, dim - 1> index,
                                   std::integral_constant<std::size_t, N>);

        /// Recursive construction from a level cell list for the dimension 0
        template<typename TIntervalList>
        void initFromLevelCellList(TIntervalList const &interval_list,
                                   const std::array<coord_index_t, dim - 1>& index,
                                   std::integral_constant<std::size_t, 0>);

    private:
        std::array<std::vector<interval_t>, dim> m_cells; ///< All intervals in every direction
        std::array<std::vector<std::size_t>, dim - 1> m_offsets; ///< Offsets in interval list for each dim > 1
        std::size_t m_level;
    };

    ////////////////////////////////////////
    // LevelCellArray_iterator definition //
    ////////////////////////////////////////
    namespace detail
    {
        template <class LCA, bool is_const>
        struct LevelCellArray_iterator_types;

        template <class LCA>
        struct LevelCellArray_iterator_types<LCA, true>
        {
            using value_type = typename LCA::interval_t;
            using index_type = std::vector<value_type>;
            using index_type_iterator = typename std::vector<value_type>::const_iterator;
            using const_index_type_iterator = typename std::vector<value_type>::const_iterator;
            using reference = const value_type&;
            using pointer = const value_type*;
            using difference_type = std::ptrdiff_t;
        };

        template <class LCA>
        struct LevelCellArray_iterator_types<LCA, false>
        {
            using value_type = typename LCA::interval_t;
            using index_type = std::vector<value_type>;
            using index_type_iterator = typename std::vector<value_type>::iterator;
            using const_index_type_iterator = typename std::vector<value_type>::const_iterator;
            using reference = value_type&;
            using pointer = value_type*;
            using difference_type = std::ptrdiff_t;
        };
    }

    template <class LCA, bool is_const>
    class LevelCellArray_iterator: public xtl::xrandom_access_iterator_base3<LevelCellArray_iterator<LCA, is_const>,
                                                                             detail::LevelCellArray_iterator_types<LCA, is_const>>
    {
    public:

        static constexpr std::size_t dim = LCA::dim;
        using self_type = LevelCellArray_iterator<LCA, is_const>;
        using iterator_types = detail::LevelCellArray_iterator_types<LCA, is_const>;
        using value_type = typename iterator_types::value_type;
        using reference = typename iterator_types::reference;
        using pointer = typename iterator_types::pointer;
        using coord_index_t = typename value_type::coord_index_t;

        using index_type = typename iterator_types::index_type;
        using index_type_iterator = std::array<typename iterator_types::index_type_iterator, dim>;
        using const_index_type_iterator = typename iterator_types::const_index_type_iterator;

        using offset_type = std::vector<std::size_t>;
        using offset_type_iterator = std::array<typename offset_type::const_iterator, dim - 1>;

        using coord_type = xt::xtensor_fixed<coord_index_t, xt::xshape<dim-1>>;

        using difference_type = typename iterator_types::difference_type;
        using iterator_category = std::random_access_iterator_tag;

        LevelCellArray_iterator(LCA& lca,
                                offset_type_iterator&& offset_index,
                                index_type_iterator&& coord_index,
                                coord_type&& index);

        self_type& operator++();
        self_type& operator--();

        self_type& operator+=(difference_type n);
        self_type& operator-=(difference_type n);

        difference_type operator-(const self_type& rhs) const;

        reference operator*() const;
        pointer operator->() const;
        const coord_type& index() const;

        bool equal(const self_type& rhs) const;
        bool less_than(const self_type& rhs) const;

    private:
        offset_type_iterator m_offset_index;
        index_type_iterator m_current_index;
        mutable coord_type m_index;
        LCA* p_lca;
    };

    ///////////////////////////////////
    // LevelCellArray implementation //
    ///////////////////////////////////
    template<std::size_t Dim, class TInterval>
    inline LevelCellArray<Dim, TInterval>::LevelCellArray(const LevelCellList<Dim, TInterval>& lcl)
    {
        /* Estimating reservation size
         *
         * NOTE: the estimation takes time, more than the time needed for
         * reallocating the vectors... Maybe 2 other solutions:
         * - (highly) overestimating the needed size since the memory will be
         * actually allocated only when touched (at least under Linux)
         * - cnt_x and cnt_yz updated in LevelCellList during the filling
         * process
         *
         * NOTE2: in fact, hard setting the optimal values for cnt_x and cnt_yz
         * doesn't speedup things, strang...
         */
        m_level = lcl.level();
        if (!lcl.empty())
        {
            // Filling cells and offsets from the level cell list
            initFromLevelCellList(lcl.grid_yz(), {}, std::integral_constant<std::size_t, dim - 1>{});
            // Additionnal offset so that [m_offset[i], m_offset[i+1][ is always valid.
            for (std::size_t d = 0; d < dim - 1; ++d)
            {
                m_offsets[d].emplace_back(m_cells[d].size());
            }
        }
    }

    template<std::size_t Dim, class TInterval>
    inline LevelCellArray<Dim, TInterval>::LevelCellArray(std::size_t level, Box<coord_index_t, dim> box)
    : m_level{level}
    {
        auto dimensions = box.length();
        auto start = box.min_corner();
        auto end = box.max_corner();

        std::size_t size = 1;
        for (std::size_t d = dim - 1; d > 0; --d)
        {
            m_offsets[d - 1].resize((dimensions[d] * size) + 1);
            for (std::size_t i = 0; i < (dimensions[d] * size) + 1; ++i)
            {
                m_offsets[d - 1][i] = i;
            }
            m_cells[d].resize(size);
            for (std::size_t i = 0; i < size; ++i)
            {
                m_cells[d][i] = { start[d], end[d],
                                  static_cast<index_t>(m_offsets[d - 1][i * dimensions[d]]) - start[d]
                                };
            }
            size *= dimensions[d];
        }

        m_cells[0].resize(size);
        for (std::size_t i = 0; i < size; ++i)
        {
            m_cells[0][i] = { start[0], end[0],
                              static_cast<index_t>(i * dimensions[0]) - start[0]
                            };
        }
    }

    template<std::size_t Dim, class TInterval>
    inline LevelCellArray<Dim, TInterval>::LevelCellArray(std::size_t level)
        : m_level{level}
    {}

    template<std::size_t Dim, class TInterval>
    inline auto LevelCellArray<Dim, TInterval>::begin() -> iterator
    {
        typename iterator::offset_type_iterator offset_index;
        typename iterator::index_type_iterator current_index;
        typename iterator::coord_type index;

        for(std::size_t d = 0; d < dim; ++d)
        {
            current_index[d] = m_cells[d].begin();
        }

        for(std::size_t d = 0; d < dim-1; ++d)
        {
            offset_index[d] = m_offsets[d].cbegin();
            index[d] = current_index[d + 1]->start;
        }
        return iterator(*this, std::move(offset_index), std::move(current_index), std::move(index));
    }

    template<std::size_t Dim, class TInterval>
    inline auto LevelCellArray<Dim, TInterval>::end() -> iterator
    {
        typename iterator::offset_type_iterator offset_index;
        typename iterator::index_type_iterator current_index;
        typename iterator::coord_type index;

        for(std::size_t d = 0; d < dim; ++d)
        {
            current_index[d] = m_cells[d].end() - 1;
        }
        ++current_index[0];

        for(std::size_t d = 0; d < dim-1; ++d)
        {
            offset_index[d] = m_offsets[d].cend() - 2;
            index[d] = current_index[d + 1]->end - 1;
        }

        return iterator(*this, std::move(offset_index), std::move(current_index), std::move(index));
    }

    template<std::size_t Dim, class TInterval>
    inline auto LevelCellArray<Dim, TInterval>::cbegin() const -> const_iterator
    {
        typename const_iterator::offset_type_iterator offset_index;
        typename const_iterator::index_type_iterator current_index;
        typename const_iterator::coord_type index;

        for(std::size_t d = 0; d < dim; ++d)
        {
            current_index[d] = m_cells[d].cbegin();
        }

        for(std::size_t d = 0; d < dim-1; ++d)
        {
            offset_index[d] = m_offsets[d].cbegin();
            index[d] = current_index[d + 1]->start;
        }
        return const_iterator(*this, std::move(offset_index), std::move(current_index), std::move(index));
    }

    template<std::size_t Dim, class TInterval>
    inline auto LevelCellArray<Dim, TInterval>::cend() const -> const_iterator
    {
        typename const_iterator::offset_type_iterator offset_index;
        typename const_iterator::index_type_iterator current_index;
        typename const_iterator::coord_type index;

        for(std::size_t d = 0; d < dim; ++d)
        {
            current_index[d] = m_cells[d].cend() - 1;
        }
        ++current_index[0];

        for(std::size_t d = 0; d < dim-1; ++d)
        {
            offset_index[d] = m_offsets[d].cend() - 2;
            index[d] = current_index[d + 1]->end - 1;
        }

        return const_iterator(*this, std::move(offset_index), std::move(current_index), std::move(index));
    }

    template<std::size_t Dim, class TInterval>
    inline auto LevelCellArray<Dim, TInterval>::rcend() const -> const_reverse_iterator
    {
        return const_reverse_iterator(cbegin());
    }

    template<std::size_t Dim, class TInterval>
    inline auto  LevelCellArray<Dim, TInterval>::rcbegin() const -> const_reverse_iterator
    {
        return const_reverse_iterator(cend());
    }

    template<std::size_t Dim, class TInterval>
    inline bool LevelCellArray<Dim, TInterval>::empty() const
    {
        return m_cells[0].empty();
    }

    template<std::size_t Dim, class TInterval>
    inline auto LevelCellArray<Dim, TInterval>::shape() const
    {
        std::array<coord_index_t, dim> output;
        for (std::size_t d = 0; d < dim; ++d)
        {
            output[d] = m_cells[d].size();
        }
        return output;
    }

    template<std::size_t Dim, class TInterval>
    inline std::size_t LevelCellArray<Dim, TInterval>::nb_cells() const
    {
        auto op = [](auto &&init, auto const &interval)
        {
            return std::move(init) + interval.size();
        };

        return std::accumulate(m_cells[0].cbegin(), m_cells[0].cend(), std::size_t(0), op);
    }

    template<std::size_t Dim, class TInterval>
    inline std::size_t LevelCellArray<Dim, TInterval>::level() const
    {
        return m_level;
    }

    template<std::size_t Dim, class TInterval>
    inline auto LevelCellArray<Dim, TInterval>::operator[](std::size_t d) const -> const std::vector<interval_t>&
    {
        return m_cells[d];
    }

    template<std::size_t Dim, class TInterval>
    inline auto LevelCellArray<Dim, TInterval>::operator[](std::size_t d) -> std::vector<interval_t>&
    {
        return m_cells[d];
    }

    template<std::size_t Dim, class TInterval>
    inline const std::vector<std::size_t>& LevelCellArray<Dim, TInterval>::offsets(std::size_t d) const
    {
        assert(d > 0);
        return m_offsets[d - 1];
    }

    template<std::size_t Dim, class TInterval>
    inline std::vector<std::size_t>& LevelCellArray<Dim, TInterval>::offsets(std::size_t d)
    {
        assert(d > 0);
        return m_offsets[d - 1];
    }

    template<std::size_t Dim, class TInterval>
    template<typename TGrid, std::size_t N>
    inline void LevelCellArray<Dim, TInterval>::initFromLevelCellList(TGrid const &grid,
                                                                      std::array<coord_index_t, dim - 1> index,
                                                                      std::integral_constant<std::size_t, N>)
    {
        // Working interval
        interval_t curr_interval(0, 0, 0);

        // For each position along the Nth dimension
        for (const auto& point : grid)
        {
            // Coordinate along the Nth dimension
            const auto i = point.first;

            // Recursive call on the current position for the (N-1)th dimension
            index[N - 1] = i;
            const std::size_t previous_offset = m_cells[N - 1].size();
            initFromLevelCellList(point.second, index, std::integral_constant<std::size_t, N - 1>{});

            /* Since we move on a sparse storage, each coordinate have non-empty
             * co-dimensions So the question is, are we continuing an existing
             * interval or have we jump to another one.
             *
             * WARNING: we are supposing that the sparse array of dimension
             * dim-1 has no empty entry. Otherwise, we should check that the
             * recursive call has do something by comparing previous_offset
             * with the size of m_cells[N-1].
             */
            if (curr_interval.is_valid())
            {
                // If the coordinate has jump out of the current interval
                if (i > curr_interval.end)
                {
                    // Adding the previous interval...
                    m_cells[N].emplace_back(curr_interval);

                    // ... and creating a new one.
                    curr_interval = interval_t(i, i + 1, static_cast<index_t>(m_offsets[N - 1].size()) - i);
                }
                else
                {
                    // Otherwise, we are just continuing the current interval
                    ++curr_interval.end;
                }
            }
            else
            {
                // If there is no current interval (at the beginning of the
                // loop) we create a new one.
                curr_interval = interval_t(i, i + 1, static_cast<index_t>(m_offsets[N - 1].size()) - i);
            }

            // Updating m_offsets (at each iteration since we are always
            // updating an interval)
            m_offsets[N - 1].emplace_back(previous_offset);
        }

        // Adding the working interval if valid
        if (curr_interval.is_valid())
        {
            m_cells[N].emplace_back(curr_interval);
        }
    }

    template<std::size_t Dim, class TInterval>
    template<typename TIntervalList>
    inline void LevelCellArray<Dim, TInterval>::initFromLevelCellList(TIntervalList const &interval_list,
                                                                      const std::array<coord_index_t, dim - 1>& /* index */,
                                                                      std::integral_constant<std::size_t, 0>)
    {
        // Along the X axis, simply copy the intervals in cells[0]
        std::copy(interval_list.begin(), interval_list.end(), std::back_inserter(m_cells[0]));
    }

    template<std::size_t Dim, class TInterval>
    inline void LevelCellArray<Dim, TInterval>::to_stream(std::ostream& os) const
    {
        for (std::size_t d = 0; d < dim; ++d)
        {
            fmt::print(os, fmt::format(fmt::emphasis::bold, "{:>10}", fmt::format("dim {}\n", d)));

            fmt::print(os, fmt::format("{:>20}", "cells = "));
            fmt::print(os, fmt::format("{}\n\n", fmt::join(m_cells[d], " ")));

            if (d > 0)
            {
                fmt::print(os, fmt::format("{:>20}", "offsets = "));
                fmt::print(os, fmt::format("{}\n\n", fmt::join(m_offsets[d-1], " ")));
            }
        }
    }

    template<std::size_t Dim, class TInterval>
    inline bool operator==(const LevelCellArray<Dim, TInterval>& lca_1,
                           const LevelCellArray<Dim, TInterval>& lca_2)
    {
        if (lca_1.level() != lca_2.level())
        {
            return false;
        }

        if (lca_1.shape() != lca_2.shape())
        {
            return false;
        }

        for (std::size_t i = 0; i < Dim; ++i)
        {
            if (lca_1[i] != lca_2[i])
            {
                return false;
            }
        }

        for (std::size_t i = 1; i < Dim; ++i)
        {
            if (lca_1.offsets(i) != lca_2.offsets(i))
            {
                return false;
            }
        }
        return true;
    }

    template<std::size_t Dim, class TInterval>
    inline std::ostream& operator<<(std::ostream &out, LevelCellArray<Dim, TInterval> const &level_cell_array)
    {
        level_cell_array.to_stream(out);
        return out;
    }

    ////////////////////////////////////////////
    // LevelCellArray_iterator implementation //
    ////////////////////////////////////////////

    template <class LCA, bool is_const>
    inline LevelCellArray_iterator<LCA, is_const>::LevelCellArray_iterator(LCA& lca,
                                                                           offset_type_iterator&& offset_index,
                                                                           index_type_iterator&& current_index,
                                                                           coord_type&& index)
    : p_lca(&lca)
    , m_offset_index(std::move(offset_index))
    , m_current_index(std::move(current_index))
    , m_index(std::move(index))
    {}

    template <class LCA, bool is_const>
    inline auto LevelCellArray_iterator<LCA, is_const>::operator++() -> self_type&
    {
        if (m_current_index[0] == (*p_lca)[0].end())
        {
            return *this;
        }
        ++m_current_index[0];

        for (std::size_t d = 0; d < m_current_index.size() - 1; ++d)
        {
            auto dst = static_cast<std::size_t>(std::distance((*p_lca)[d].cbegin(), static_cast<const_index_type_iterator>(m_current_index[d])));
            if (dst == *(m_offset_index[d] + 1))
            {
                ++m_offset_index[d];
                ++m_index[d];
                if (m_index[d] == m_current_index[d+1]->end)
                {
                    ++m_current_index[d+1];
                    m_index[d] = m_current_index[d+1]->start;
                }
            }
            else
            {
                break;
            }
        }
        return *this;
    }

    template <class LCA, bool is_const>
    inline auto LevelCellArray_iterator<LCA, is_const>::operator--() -> self_type&
    {
        if (m_current_index[0] == (*p_lca)[0].begin())
        {
            --m_current_index[0];
            return *this;
        }
        --m_current_index[0];

        for (std::size_t d = 0; d < m_current_index.size() - 1; ++d)
        {
            auto dst = static_cast<std::size_t>(std::distance((*p_lca)[d].cbegin(), static_cast<const_index_type_iterator>(m_current_index[d])));
            if (dst == *m_offset_index[d] - 1)
            {
                --m_offset_index[d];
                if (m_index[d] == m_current_index[d+1]->start)
                {
                    --m_current_index[d+1];
                    m_index[d] = m_current_index[d+1]->end - 1;
                }
                else
                {
                    --m_index[d];
                }
            }
            else
            {
                break;
            }
        }
        return *this;
    }

    template <class LCA, bool is_const>
    inline auto LevelCellArray_iterator<LCA, is_const>::operator+=(difference_type n) -> self_type&
    {
        for (difference_type i = 0; i < n; ++i)
        {
            ++(*this);
        }
        return *this;
    }

    template <class LCA, bool is_const>
    inline auto LevelCellArray_iterator<LCA, is_const>::operator-=(difference_type n) -> self_type&
    {
        for (difference_type i = 0; i < n; ++i)
        {
            --(*this);
        }
        return *this;
    }

    template <class LCA, bool is_const>
    inline auto LevelCellArray_iterator<LCA, is_const>::operator-(const self_type& rhs) const -> difference_type
    {
        return m_current_index[0] - rhs.m_current_index[0];
    }

    template <class LCA, bool is_const>
    inline auto LevelCellArray_iterator<LCA, is_const>::operator*() const -> reference
    {
        return *(m_current_index[0]);
    }

    template <class LCA, bool is_const>
    inline auto LevelCellArray_iterator<LCA, is_const>::operator->() const -> pointer
    {
        return &(this->operator*());
    }

    template <class LCA, bool is_const>
    inline auto LevelCellArray_iterator<LCA, is_const>::index() const -> const coord_type&
    {
        return m_index;
    }

    template <class LCA, bool is_const>
    inline bool LevelCellArray_iterator<LCA, is_const>::equal(const self_type& rhs) const
    {
        return p_lca == rhs.p_lca && m_current_index[0] == rhs.m_current_index[0];
    }

    template <class LCA, bool is_const>
    inline bool LevelCellArray_iterator<LCA, is_const>::less_than(const self_type& rhs) const
    {
        return p_lca == rhs.p_lca && m_current_index[0] < rhs.m_current_index[0];
    }

    template <class LCA, bool is_const>
    inline bool operator==(const LevelCellArray_iterator<LCA, is_const>& it1,
                           const LevelCellArray_iterator<LCA, is_const>& it2)
    {
        return it1.equal(it2);
    }

    template <class LCA, bool is_const>
    inline bool operator<(const LevelCellArray_iterator<LCA, is_const>& it1,
                          const LevelCellArray_iterator<LCA, is_const>& it2)
    {
        return it1.less_than(it2);
    }

    template <class LCA, bool is_const>
    inline bool operator==(const std::reverse_iterator<LevelCellArray_iterator<LCA, is_const>>& it1,
                           const std::reverse_iterator<LevelCellArray_iterator<LCA, is_const>>& it2)
    {
        return it1.base().equal(it2.base());
    }
}
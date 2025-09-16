//
// Created by Willem Lambooy on 15/04/2025.
//

#ifndef SIDB_SKELETON_HPP
#define SIDB_SKELETON_HPP

#include <fiction/layouts/bounding_box.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include <fiction/technology/fcn_gate_library.hpp>
#include <fiction/traits.hpp>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fiction
{

template <typename Lyt, uint16_t GateSizeX, uint16_t GateSizeY>
class sidb_skeleton : public Lyt
{
  public:
    template <typename... Args>
    explicit sidb_skeleton(Args&&... args) noexcept(noexcept(Lyt(std::forward<Args>(args)...))) :
            Lyt(std::forward<Args>(args)...)
    {}

    void determine_bounding_box() noexcept
    {
        bbox.update_bounding_box();
    }

    void collect_range(const coordinate<Lyt>& min, const coordinate<Lyt>& max) noexcept
    {
        std::vector<cell<Lyt>> collected_cells{};

        this->foreach_cell(
            [&](const cell<Lyt>& c)
            {
                // std::cout << "x: " << c.x << ", y: " << c.y << std::endl;
                // std::cout << "max_x: " << max.x << ", max_y: " << max.y << std::endl;
                // std::cout << "min_x: " << min.x << ", min_y: " << min.y << std::endl;
                if (c.x >= min.x && c.x <= max.x && c.y >= min.y && c.y <= max.y)
                {
                    // std::cout <<"yes"<< std::endl;
                    collected_cells.emplace_back(c);
                }
            });

        std::sort(collected_cells.begin(), collected_cells.end(),
                  [](const cell<Lyt>& c1, const cell<Lyt>& c2) noexcept { return c1.y < c2.y; });

        // for (const auto& c : collected_cells)
        // {
        //     std::cout << "x: " << c.x << ", y: " << c.y << std::endl;
        // }
        // std::cout << std::endl;

        typename std::vector<wire>::iterator it = find_wire_in_input_wires(collected_cells);

        if (it == input_wires.end())
        {
            // std::cout << "new inp" << std::endl;
            input_wires.emplace_back(std::move(collected_cells));
        }
        else
        {
            // std::cout << "new out" << std::endl;
            output_wires.emplace_back(std::move(collected_cells));
            input_wires.erase(it);
        }
    }

    void sort_wires() noexcept
    {
        assert(!input_wires.empty() && "no input wires");
        assert(!output_wires.empty() && "no output wires");

        std::sort(input_wires.begin(), input_wires.end(),
                  [](const std::vector<cell<Lyt>>& w1, const std::vector<cell<Lyt>>& w2) noexcept
                  { return w1.front().x < w2.front().x; });
        std::sort(output_wires.begin(), output_wires.end(),
                  [](const std::vector<cell<Lyt>>& w1, const std::vector<cell<Lyt>>& w2) noexcept
                  { return w1.back().x < w2.back().x; });

        if (input_wires.size() == 1)
        {
            if (input_wires.front().front().x < bbox.get_min().x + bbox.get_x_size() / 2)
            {
                port_direction_indices.erase(port_direction{port_direction::cardinal::NORTH_EAST});
            }
            else
            {
                port_direction_indices.erase(port_direction{port_direction::cardinal::NORTH_WEST});
                port_direction_indices[port_direction{port_direction::cardinal::NORTH_EAST}].second = 0;
            }
        }

        if (output_wires.size() == 1)
        {
            if (output_wires.front().front().x < bbox.get_min().x + bbox.get_x_size() / 2)
            {
                port_direction_indices.erase(port_direction{port_direction::cardinal::SOUTH_EAST});
            }
            else
            {
                port_direction_indices.erase(port_direction{port_direction::cardinal::SOUTH_WEST});
                port_direction_indices[port_direction{port_direction::cardinal::SOUTH_EAST}].second = 0;
            }
        }
    }

    typename fcn_gate_library<sidb_technology, GateSizeX, GateSizeY>::fcn_gate
    make_skeleton(const port_list<port_direction>& port_routing) const noexcept
    {
        assert(bbox.get_x_size() <= GateSizeX && "GateSizeX too small");
        assert(bbox.get_y_size() <= GateSizeY && "GateSizeY too small");

        assert(output_wires.front().size() >= 3 && "not enough cells in output wire");

        typename fcn_gate_library<sidb_technology, GateSizeX, GateSizeY>::fcn_gate g =
            fcn_gate_library<sidb_technology, GateSizeX, GateSizeY>::EMPTY_GATE;

        static_assert(std::is_same_v<coordinate<Lyt>, offset::ucoord_t>, "Needs to have fiction coordinates");

        decltype(coordinate<Lyt>{}.x) offset_x = (GateSizeX - bbox.get_x_size()) / 2;

        for (const auto& [pd, t] : port_direction_indices)
        {
            const auto& [w, ix] = t;

            if (w == false && port_routing.inp.count(pd) != 0)
            {
                assert(input_wires.at(ix).size() >= 2 && "not enough cells in input wire");

                for (const cell<Lyt>& c : input_wires.at(ix))
                {
                    const sidb_technology::cell_type cell_type =
                        c == input_wires.at(ix).front() || c == *std::next(input_wires.at(ix).cbegin(), 1) ?
                            sidb_technology::cell_type::INPUT :
                            sidb_technology::cell_type::NORMAL;

                    g[c.y - bbox.get_min().y][c.x - bbox.get_min().x + offset_x] = cell_type;
                }
            }

            if (w == true && port_routing.out.count(pd) != 0)
            {
                assert(output_wires.at(ix).size() >= 3 && "not enough cells in output wire");

                for (const cell<Lyt>& c : output_wires.at(ix))
                {
                    const sidb_technology::cell_type cell_type =
                        c == output_wires.at(ix).back() ? sidb_technology::cell_type::OUTPUT_PERTURBER :
                        c == *std::prev(output_wires.at(ix).cend(), 2) ||
                                c == *std::prev(output_wires.at(ix).cend(), 3) ?
                                                          sidb_technology::cell_type::OUTPUT :
                                                          sidb_technology::cell_type::NORMAL;

                    g[c.y - bbox.get_min().y][c.x - bbox.get_min().x + offset_x] = cell_type;
                }
            }
        }

        return g;
    }

    std::vector<std::vector<cell<Lyt>>> get_input_wires() const noexcept
    {
        return input_wires;
    }

    std::vector<std::vector<cell<Lyt>>> get_output_wires() const noexcept
    {
        return output_wires;
    }

  private:
    using wire = std::vector<cell<Lyt>>;

    bounding_box_2d<Lyt> bbox{*this};

    std::vector<wire> input_wires{};
    std::vector<wire> output_wires{};

    std::unordered_map<port_direction, std::pair<bool, uint64_t>> port_direction_indices = {
        {port_direction{port_direction::cardinal::NORTH_WEST}, {false, 0}},
        {port_direction{port_direction::cardinal::NORTH_EAST}, {false, 1}},
        {port_direction{port_direction::cardinal::SOUTH_WEST}, {true, 0}},
        {port_direction{port_direction::cardinal::SOUTH_EAST}, {true, 1}}};

    typename std::vector<wire>::iterator find_wire_in_input_wires(const wire& w) noexcept
    {
        typename std::vector<wire>::iterator it = input_wires.begin();

        for (; it != input_wires.end(); ++it)
        {
            if (w.size() != it->size())
            {
                continue;
            }

            bool failed = false;

            for (uint64_t i = 0; !failed && i < w.size(); i++)
            {
                if (w[i] != (*it)[i])
                {
                    failed = true;
                }
            }

            if (!failed)
            {
                return it;
            }
        }

        return input_wires.end();
    }
};

}  // namespace fiction

#endif  // SIDB_SKELETON_HPP

//
// Created by Jan Drewniok on 26.09.23.
//

#ifndef FICTION_SIDB_SKELETON_GATE_LIBRARY_HPP
#define FICTION_SIDB_SKELETON_GATE_LIBRARY_HPP

#include "fiction/technology/cell_ports.hpp"
#include "fiction/technology/cell_technologies.hpp"
#include "fiction/technology/fcn_gate_library.hpp"
#include "fiction/technology/sidb_skeletons/sidb_skeleton.hpp"
#include "fiction/traits.hpp"

namespace fiction
{

/**
 * This library contains SiDB I/O wires designed for both 1- and 2-input functions.
 * Each wire comprises 2 BDL pairs. The library contains all mirrored versions, a double wire and a crossing.
 */
template <uint16_t GateSizeX, uint16_t GateSizeY>
class sidb_skeleton_gate_library
        : public fcn_gate_library<sidb_technology, GateSizeX, GateSizeY>  // width and height of a hexagon
{
  public:
    using sidb_skeleton_t =
        sidb_skeleton<sidb_100_cell_clk_lyt, fcn_gate_library<sidb_technology, GateSizeX, GateSizeY>::gate_x_size(),
                      fcn_gate_library<sidb_technology, GateSizeX, GateSizeY>::gate_y_size()>;

    sidb_skeleton_gate_library(sidb_skeleton_t&& skeleton_prod) : skeleton_producer{std::move(skeleton_prod)} {}

    /**
     * Overrides the corresponding function in fcn_gate_library. Given a tile `t`, this function takes all necessary
     * information from the stored grid into account to choose the correct fcn_gate representation for that tile. May it
     * be a gate or wires. Rotation and special marks like input and output, const cells etc. are computed additionally.
     *
     * @tparam GateLyt Pointy-top hexagonal gate-level layout type.
     * @param lyt Layout that hosts tile `t`.
     * @param t Tile to be realized as a Bestagon skeleton gate.
     * @return Bestagon skeleton gate representation of `t` including mirroring.
     */
    template <typename GateLyt>
    [[nodiscard]] typename fcn_gate_library<sidb_technology, GateSizeX, GateSizeY>::fcn_gate
    set_up_gate(const GateLyt& lyt, const tile<GateLyt>& t) const
    {
        static_assert(is_gate_level_layout_v<GateLyt>, "GateLyt must be a gate-level layout");
        static_assert(is_hexagonal_layout_v<GateLyt>, "GateLyt must be a hexagonal layout");
        static_assert(has_pointy_top_hex_orientation_v<GateLyt>, "GateLyt must be a pointy-top hexagonal layout");

        const auto n = lyt.get_node(t);
        const auto p = determine_port_routing(lyt, t);

        try
        {
            if constexpr (fiction::has_is_buf_v<GateLyt>)
            {
                if (lyt.is_buf(n))
                {
                    if (lyt.is_ground_layer(t))
                    {
                        // crossing case
                        if (const auto at = lyt.above(t); (t != at) && lyt.is_wire_tile(at))
                        {
                            // two possible options: actual crossover and (parallel) hourglass wire
                            const auto pa = determine_port_routing(lyt, at);

                            port_list<port_direction> combo_p{};

                            for (const auto& pd : p.inp)
                            {
                                combo_p.inp.emplace(pd);
                            }

                            for (const auto& pd : pa.inp)
                            {
                                combo_p.inp.emplace(pd);
                            }

                            for (const auto& pd : p.out)
                            {
                                combo_p.out.emplace(pd);
                            }

                            for (const auto& pd : pa.out)
                            {
                                combo_p.out.emplace(pd);
                            }

                            return skeleton_producer.make_skeleton(combo_p);
                        }
                        // regular wire: look-up in the wire_map

                        return skeleton_producer.make_skeleton(p);
                    }

                    return fcn_gate_library<sidb_technology, GateSizeX, GateSizeY>::EMPTY_GATE;
                }
            }

            return skeleton_producer.make_skeleton(p);
        }
        catch (const std::out_of_range&)
        {
            throw unsupported_gate_orientation_exception(t, p);
        }

        throw unsupported_gate_type_exception(t);
    }

    template <typename Lyt>
    [[nodiscard]] static port_list<port_direction> determine_port_routing(const Lyt& lyt, const tile<Lyt>& t) noexcept
    {
        port_list<port_direction> p{};

        // determine incoming connector ports
        if (lyt.has_north_eastern_incoming_signal(t))
        {
            p.inp.emplace(port_direction::cardinal::NORTH_EAST);
        }
        if (lyt.has_north_western_incoming_signal(t))
        {
            p.inp.emplace(port_direction::cardinal::NORTH_WEST);
        }

        // determine outgoing connector ports
        if (lyt.has_south_eastern_outgoing_signal(t))
        {
            p.out.emplace(port_direction::cardinal::SOUTH_EAST);
        }
        if (lyt.has_south_western_outgoing_signal(t))
        {
            p.out.emplace(port_direction::cardinal::SOUTH_WEST);
        }

        // gates without connector ports

        // 1-input functions
        if (const auto n = lyt.get_node(t); lyt.is_pi(n) || lyt.is_po(n) || lyt.is_buf(n) || lyt.is_inv(n))
        {
            if (lyt.has_no_incoming_signal(t))
            {
                p.inp.emplace(port_direction::cardinal::NORTH_WEST);
            }
            if (lyt.has_no_outgoing_signal(t))
            {
                p.out.emplace(port_direction::cardinal::SOUTH_EAST);
            }
        }
        else  // 2-input functions
        {
            if (lyt.has_no_incoming_signal(t))
            {
                p.inp.emplace(port_direction::cardinal::NORTH_WEST);
                p.inp.emplace(port_direction::cardinal::NORTH_EAST);
            }
            if (lyt.has_no_outgoing_signal(t))
            {
                p.out.emplace(port_direction::cardinal::SOUTH_EAST);
            }
        }

        return p;
    }

  private:
    sidb_skeleton_t skeleton_producer;
};

}  // namespace fiction

#endif  // FICTION_SIDB_SKELETON_GATE_LIBRARY_HPP

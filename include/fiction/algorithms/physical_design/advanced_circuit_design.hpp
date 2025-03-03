//
// Created by Willem Lambooy on 27.01.2025.
//

#ifndef ADVANCED_CIRCUIT_DESIGN_HPP
#define ADVANCED_CIRCUIT_DESIGN_HPP

#include "fiction/algorithms/physical_design/exact.hpp"
#include "fiction/technology/sidb_defect_surface.hpp"
#include "fiction/technology/sidb_on_the_fly_gate_library.hpp"
#include "fiction/technology/sidb_skeleton_bestagon_library.hpp"
#include "fiction/technology/sidb_surface_analysis.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"
#include "fiction/utils/gate_design_utils.hpp"

#include <mockturtle/utils/stopwatch.hpp>

#include <cstdint>
#include <optional>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

#include <amxintrin.h>
#include <eigen3/Eigen/src/Core/util/IntegralConstant.h>

namespace fiction
{

/**
 * This struct stores the parameters to design an SiDB circuit on a defective surface.
 *
 * @tparam CellLyt SiDB cell-level layout type.
 */
template <typename CellLyt>
struct advanced_circuit_design_params
{
    /**
     * Parameters for the SiDB on-the-fly gate library.
     */
    sidb_on_the_fly_gate_library_params<CellLyt> sidb_on_the_fly_gate_library_parameters = {};
    /**
     * Parameters for the *exact* placement and routing algorithm.
     */
    exact_physical_design_params exact_design_parameters = {};
    uint64_t                     max_num_trials          = 500;
};

/**
 * Statistics for the on-the-fly defect-aware circuit design.
 */
template <typename GateLyt>
struct advanced_circuit_design_stats
{
    /**
     * The total runtime of the operational domain computation.
     */
    mockturtle::stopwatch<>::duration time_total{0};
    /**
     * The `stats` of the *exact* algorithm.
     */
    exact_physical_design_stats exact_stats{};
    /**
     * The gate-level layout after P&R.
     */
    std::optional<GateLyt> gate_layout{};
};

namespace detail
{

template <typename Ntk, typename CellLyt, typename GateLyt>
class advanced_circuit_design_impl
{
  public:
    advanced_circuit_design_impl(const Ntk& ntk, const advanced_circuit_design_params<CellLyt>& design_params,
                                 const GateLyt& tiling, advanced_circuit_design_stats<GateLyt>& st) :
            lattice_tiling{tiling},
            network{ntk},
            params{design_params},
            stats{st}
    {}

    [[nodiscard]] sidb_defect_surface<CellLyt> design_circuit_on_defective_surface()
    {
        const mockturtle::stopwatch stop{stats.time_total};

        std::optional<GateLyt> gate_lyt = std::nullopt;

        CellLyt lyt{};

        // generating the blacklist based on neutral defects. The long-range electrostatic influence of charged defects
        // is not considered as gates are designed on-the-fly.
        auto black_list = sidb_surface_analysis<sidb_skeleton_bestagon_library>(
            lattice_tiling, params.sidb_on_the_fly_gate_library_parameters.defect_surface, std::make_pair(0, 0));

        while (!gate_lyt.has_value())
        {
            // P&R with *exact* and the pre-determined blacklist
            gate_lyt =
                exact_with_blacklist<GateLyt>(network, black_list, params.exact_design_parameters, &stats.exact_stats);

            if (!gate_lyt.has_value())
            {
                // P&R was unsuccessful
                break;
            }

            std::unordered_map<mockturtle::node<Ntk>, std::vector<sidb_skeleton_bestagon_library::fcn_gate>>
                operational_gate_designs{};

            try
            {
                gate_lyt->foreach_node(
                    [&, this](const auto& n, [[maybe_unused]] auto i)
                    {
                        if (!gate_lyt->is_constant(n))
                        {
                            const auto t = gate_lyt->get_tile(n);

                            operational_gate_designs[n] = sidb_on_the_fly_gate_library::set_up_gates<
                                GateLyt, CellLyt, sidb_on_the_fly_gate_library_params<CellLyt>>(
                                gate_lyt, t, params.sidb_on_the_fly_gate_library_parameters);
                        }
                    });

                prune_gate_designs(*gate_lyt, operational_gate_designs);
            }

            // on-the-fly gate design was unsuccessful at a certain tile. Hence, this tile-gate pair is added to the
            // blacklist and the process is rerun.
            catch (const gate_design_exception<tt, GateLyt>& e)
            {
                black_list[e.which_tile()][e.which_truth_table()].push_back(e.which_port_list());
            }
        }

        stats.gate_layout = std::optional{gate_lyt};

        sidb_defect_surface<CellLyt> sidbs_and_defects{lyt};

        // add defects to the circuit.
        params.sidb_on_the_fly_gate_library_parameters.defect_surface.foreach_sidb_defect(
            [&sidbs_and_defects](const auto& defect)
            { sidbs_and_defects.assign_sidb_defect(defect.first, defect.second); });

        return sidbs_and_defects;
    }

  private:
    /**
     * Gate-level layout.
     */
    GateLyt lattice_tiling;
    /**
     * Network.
     */
    Ntk network;
    /**
     * Parameters for the on-the-fly circuit design.
     */
    const advanced_circuit_design_params<CellLyt> params{};
    /**
     * Statistics for the on-the-fly circuit design.
     */
    advanced_circuit_design_stats<GateLyt>& stats;
    /**
     *
     */
    void
    prune_gate_designs(const GateLyt& gate_lyt,
                       std::unordered_map<mockturtle::node<Ntk>, std::vector<sidb_skeleton_bestagon_library::fcn_gate>>&
                           operational_gate_designs) const noexcept
    {

        std::unordered_map<mockturtle::node<Ntk>, cell<CellLyt>> top_left_corner_abs{};

        gate_lyt.foreach_node(
            [&](const auto& n)
            {
                top_left_corner_abs.insert(
                    n,
                    relative_to_absolute_cell_position<sidb_on_the_fly_gate_library::gate_x_size(),
                                                       sidb_on_the_fly_gate_library::gate_y_size(), GateLyt, CellLyt>(
                        gate_lyt, gate_lyt.get_tile(n), cell<CellLyt>{0, 0}));
            });

        bool big_fixpoint = false;

        while (!big_fixpoint)
        {
            big_fixpoint = true;

            gate_lyt.foreach_node(
                [&](const auto& n)
                {
                    if (gate_lyt->is_constant(n))
                    {
                        return;
                    }

                    const auto t = gate_lyt->get_tile(n);

                    std::vector<mockturtle::node<Ntk>> tiles_selected_for_operational_assessment{};

                    gate_lyt.foreach_fanin(n,
                                           [&](const auto& before_n)
                                           {
                                               if (!gate_lyt->is_constant(before_n))  // necessary?
                                               {
                                                   tiles_selected_for_operational_assessment.push_back(before_n);
                                               }
                                           });

                    if (tiles_selected_for_operational_assessment.size() == 1)
                    {
                        return;
                    }

                    bool fixpoint = false;

                    while (!fixpoint)
                    {
                        fixpoint = true;

                        uint64_t current_iteration_gate = 0;

                        // std::vector<mockturtle::node<Ntk>> other_n_order{};
                        // other_n_order.reserve(tiles_selected_for_operational_assessment.size());
                        //
                        // std::vector<uint64_t> gate_design_indices{};
                        // gate_design_indices.reserve(tiles_selected_for_operational_assessment.size());
                        //
                        // for (const tile<GateLyt>& other_n : tiles_selected_for_operational_assessment)
                        // {
                        //     gate_design_indices.push_back(0);
                        //     other_n_order.push_back(other_n);
                        // }

                        const uint64_t num_trials = params.max_num_trials;

                        CellLyt cell_lyt{};

                        assign_gate<CellLyt, sidb_on_the_fly_gate_library, GateLyt>(
                            cell_lyt, top_left_corner_abs.at(n), operational_gate_designs.at(t).front(), gate_lyt, n);

                        uint64_t current_trial = 0;

                        while (current_trial < num_trials)
                        {
                            current_trial++;

                            CellLyt cell_lyt_clone = cell_lyt.clone();

                            for (const auto& other_n : tiles_selected_for_operational_assessment)
                            {
                                const tile<GateLyt>& other_t = gate_lyt.get_tile(other_n);

                                assign_gate<CellLyt, sidb_on_the_fly_gate_library, GateLyt>(
                                    cell_lyt, top_left_corner_abs.at(other_n),
                                    operational_gate_designs.at(other_t).at(std::uniform_int_distribution<>{
                                        0, operational_gate_designs.at(other_n).size() -
                                               1}(std::mt19937{std::random_device{}()})),
                                    gate_lyt, n);
                            }

                            if (!is_operational(cell_lyt /*function*/))
                            {
                                operational_gate_designs[n].erase(0);
                                ;
                            }
                        }

                        // while (fixpoint && current_iteration_gate < other_n_order.size())
                        // {
                        //     if (gate_design_indices.at(current_iteration_gate) ==
                        //     operational_gate_designs.at(other_n_order.at(current_iteration_gate)))
                        //     {
                        //
                        //     }
                        //
                        //     CellLyt cell_lyt{};
                        //
                        //     assign_gate<CellLyt, sidb_on_the_fly_gate_library, GateLyt>(
                        //         cell_lyt, c, operational_gate_designs.at(t).at(counter++), gate_lyt, n);
                        //
                        //     // if (gate_design_indices.at(current_iteration_gate) ==
                        //     other_n_order[current_iteration_gate]())
                        //     // if ()
                        //     // {
                        //     //     fixpoint = false;
                        //     //     continue;
                        //     // }
                        //
                        //
                        // }
                    }

                    // assign_gate<CellLyt, sidb_on_the_fly_gate_library, GateLyt>(cell_lyt, c,
                    // operational_gate_designs.at(t).at(counter++), gate_lyt, n);
                });
        }
    }
};

}  // namespace detail

/**
 *
 *
 * @tparam Ntk The type of the input network.
 * @tparam CellLyt SiDB cell-level layout type.
 * @tparam GateLyt Gate-level layout type.
 * @param ntk The input network to be mapped onto the defective surface.
 * @param lattice_tiling The lattice tiling used for the circuit design.
 * @param params The parameters used for designing the circuit, encapsulated in an
 * `advanced_circuit_design_params` object.
 * @param stats Pointer to a structure for collecting statistics. If nullptr, statistics are not collected.
 * @return A `sidb_defect_surface<CellLyt>` representing the designed circuit on the defective surface.
 */
template <typename Ntk, typename CellLyt, typename GateLyt>
[[nodiscard]] sidb_defect_surface<CellLyt>
advanced_circuit_design(const Ntk& ntk, const GateLyt& lattice_tiling,
                        const advanced_circuit_design_params<CellLyt>& params = {},
                        advanced_circuit_design_stats<GateLyt>*        stats  = nullptr)
{
    static_assert(is_gate_level_layout_v<GateLyt>, "GateLyt is not a gate-level layout");
    static_assert(is_hexagonal_layout_v<GateLyt>, "GateLyt is not a hexagonal");
    static_assert(is_cell_level_layout_v<CellLyt>, "CellLyt is not a cell-level layout");
    static_assert(has_sidb_technology_v<CellLyt>, "CellLyt is not an SiDB layout");
    static_assert(mockturtle::is_network_type_v<Ntk>, "Ntk is not a network type");

    advanced_circuit_design_stats<GateLyt> st{};

    detail::advanced_circuit_design_impl<Ntk, CellLyt, GateLyt> p{ntk, params, lattice_tiling, st};

    const auto result = p.design_circuit_on_defective_surface();

    if (stats)
    {
        *stats = st;
    }

    return result;
}

}  // namespace fiction

#endif  // advanced_CIRCUIT_DESIGN_HPP

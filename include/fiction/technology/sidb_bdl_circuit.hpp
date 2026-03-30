//
// Created by Willem Lambooy on 06/06/2025.
//

#ifndef SIDB_BDL_CIRCUIT_HPP
#define SIDB_BDL_CIRCUIT_HPP

#include "fiction/algorithms/iter/bdl_input_iterator.hpp"
#include "fiction/algorithms/physical_design/apply_gate_library.hpp"
#include "fiction/algorithms/simulation/sidb/detect_bdl_wires.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/technology/sidb_charge_state.hpp"
#include "fiction/technology/sidb_cluster_hierarchy.hpp"
#include "fiction/technology/sidb_defects.hpp"
#include "fiction/traits.hpp"

#include <fiction/io/write_sqd_layout.hpp>

#include <kitty/print.hpp>
#include <phmap.h>

#include <cmath>
#include <limits>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include <oneapi/tbb/partitioner.h>
#include <signal.h>

namespace fiction
{

template <typename CellLyt, typename GateLyt>
class sidb_bdl_circuit
{
  public:  // todo fix ALL GATES
    explicit sidb_bdl_circuit(
        const GateLyt& gate_lyt, const sidb_simulation_parameters& simulation_parameters,
        const std::unordered_map<mockturtle::node<GateLyt>, std::pair<uint8_t, uint8_t>>& sidb_counts,
        const bool print_skeleton = true) noexcept :
            gate_layout{gate_lyt.clone()},
            sim_params{simulation_parameters},
            canvasses{make_absolute_canvasses(gate_layout, sidb_counts)},
            canvasses_lyt{make_canvas_layout(gate_layout, canvasses, print_skeleton)}
    {
        write_sqd_layout(canvasses_lyt, "test.sqd");
        set_initial_gate_design_influence_bounds();
    }
    /**
     * SiDB gate-level layout.
     */
    const GateLyt gate_layout;

    const sidb_simulation_parameters sim_params;

    struct canvas
    {
        std::vector<cell<CellLyt>> positions{};
        uint8_t                    min_count;
        uint8_t                    max_count;
    };

    template <typename T>
    using foreach_node = std::unordered_map<mockturtle::node<GateLyt>, T>;

    const foreach_node<canvas> canvasses;

    const CellLyt canvasses_lyt;

    /**
     * A canvas combination is a combination of canvas positions as a vector of canvas position indices.
     */
    using canvas_combination = std::vector<std::size_t>;

    // struct cell_design
    // {
    //     canvas_combination positions;
    //
    // };

    foreach_node<std::vector<canvas_combination>> gate_designs{};

    [[nodiscard]] static foreach_node<canvas>
    make_absolute_canvasses(const GateLyt&                                   gate_lyt,
                            const foreach_node<std::pair<uint8_t, uint8_t>>& sidb_counts) noexcept
    {
        foreach_node<canvas> canvas_per_node{};

        gate_lyt.foreach_node(
            [&](const auto& n)
            {
                if (gate_lyt.is_constant(n))
                {
                    return;
                }

                const tile<GateLyt>& t = gate_lyt.get_tile(n);

                const cell<CellLyt> nw_position = relative_to_absolute_cell_position<3, 4, GateLyt, CellLyt>(
                    gate_lyt, t, {gate_lyt.is_in_odd_row(t) ? 1 : 0, 0});
                const cell<CellLyt> se_position =
                    relative_to_absolute_cell_position<3, 4, GateLyt, CellLyt>(gate_lyt, t, {2, 2});

                canvas_per_node[n].positions = all_coordinates_in_spanned_area(nw_position, se_position);
                canvas_per_node[n].min_count = sidb_counts.at(n).first;
                canvas_per_node[n].max_count = sidb_counts.at(n).second;
            });

        return canvas_per_node;
    }

    static CellLyt make_canvas_layout(
        const GateLyt& gate_lyt, const std::unordered_map<mockturtle::node<GateLyt>, canvas>& abs_canvas_at_tile,
        const bool                                                   print          = true,
        const std::optional<std::vector<mockturtle::node<GateLyt>>>& node_whitelist = std::nullopt) noexcept
    {
        CellLyt canvas_lyt{};

        gate_lyt.foreach_node(
            [&](const auto& n)
            {
                if (gate_lyt.is_constant(n) ||
                    (node_whitelist.has_value() &&
                     std::find(node_whitelist->cbegin(), node_whitelist->cend(), n) == node_whitelist->cend()))
                {
                    return;
                }

                for (const cell<CellLyt>& c : abs_canvas_at_tile.at(n).positions)
                {
                    canvas_lyt.assign_cell_type(c, sidb_technology::cell_type::LOGIC);
                    canvas_lyt.assign_cell_tile(c, {gate_lyt.get_tile(n).x, gate_lyt.get_tile(n).y, 0});
                }
            });

        if (print)
        {
            std::cout << "Canvas layout:" << std::endl;
            print_layout(canvas_lyt);
            std::cout << std::endl;
            std::cout << "tiles: ";
            gate_lyt.foreach_node(
                [&](const auto& n)
                {
                    if (gate_lyt.is_constant(n) ||
                        (node_whitelist.has_value() &&
                         std::find(node_whitelist->cbegin(), node_whitelist->cend(), n) == node_whitelist->cend()))
                    {
                        return;
                    }
                    std::cout << gate_lyt.get_tile(n) << '\t';
                });
            std::cout << std::endl;
        }

        return canvas_lyt;
    }

    [[nodiscard]] std::reference_wrapper<const std::array<double, 2>>
    get_gate_design_influence_bounds(const uint64_t input_index, const mockturtle::node<GateLyt>& n,
                                     const cell<CellLyt>& c, const mockturtle::node<GateLyt>& other_n) const noexcept
    {
        return std::cref(gate_design_influence_bounds_per_input.at(input_index).at(n).at(c).at(other_n));
    }

    void tighten_gate_influence_bounds_until_fixpoint(
        const uint64_t available_threads = 1)  // std::thread::hardware_concurrency()) noexcept
    {
        foreach_node<uint64_t> gate_design_counts{};
        gate_design_counts.reserve(gate_designs.size());

        for (const auto& [n, designs] : gate_designs)
        {
            gate_design_counts[n] = designs.size();
        }

        std::vector<uint64_t> input_indices{};
        input_indices.reserve(1 << gate_layout.num_pis());
        for (uint64_t i = 0; i < 1 << gate_layout.num_pis(); input_indices.push_back(i++))
        {}

        // tighten bounds until fixed point

        bool big_fixpoint = false;

        while (!big_fixpoint)
        {
            std::cout << "\nSTARTING FIXED POINT ITERATION " << std::endl;

            big_fixpoint = true;

            // std::shuffle(input_indices.begin(), input_indices.end(), std::mt19937(std::random_device()()));

            for (const uint64_t input_index : input_indices)
            {
                std::cout << "NARROWING INFLUENCE BOUNDS FOR INPUT INDEX " << input_index << std::endl;

                bool fixpoint = false;

                while (!fixpoint)
                {
                    fixpoint = true;

                    gate_layout.foreach_node(
                        [&](const auto& n)
                        {
                            if (gate_layout.is_constant(n))
                            {
                                return;
                            }

                            // std::cout << "\n\nnow pruning for tile " << gate_layout.get_tile(n) << std::endl;

                            foreach_node<std::set<uint64_t>> gate_nums_to_prune{};
                            std::mutex                       mutex_to_protect_gate_nums_to_prune{};

                            const uint64_t num_cells = canvasses.at(n).positions.size();

                            const uint64_t num_threads = std::min(available_threads, num_cells);

                            const uint64_t chunk_size =
                                (num_cells + num_threads - 1) / num_threads;  // Ceiling division

                            std::vector<std::thread> threads{};
                            threads.reserve(num_threads);

                            for (uint64_t i = 0; i < num_threads; ++i)
                            {
                                threads.emplace_back(
                                    [&, i]
                                    {
                                        const uint64_t cell_start_index = i * chunk_size;
                                        const uint64_t cell_end_index =
                                            std::min(cell_start_index + chunk_size, num_cells);

                                        tighten_gate_influence_bounds(input_index, n, cell_start_index, cell_end_index,
                                                                      gate_nums_to_prune,
                                                                      mutex_to_protect_gate_nums_to_prune, fixpoint);
                                    });
                            }

                            for (auto& thread : threads)
                            {
                                if (thread.joinable())
                                {
                                    thread.join();
                                }
                            }

                            gate_layout.foreach_node(
                                [&](const auto& other_n)
                                {
                                    if (gate_layout.is_constant(other_n))
                                    {
                                        return;
                                    }

                                    if (gate_nums_to_prune.count(other_n) == 0)
                                    {
                                        return;
                                    }

                                    const tile<GateLyt>& other_t = gate_layout.get_tile(other_n);
                                    if (gate_designs.at(other_n).size() == gate_nums_to_prune.at(other_n).size())
                                    {
                                        throw std::runtime_error{
                                            fmt::format("All gate designs pruned for tile {}", other_t)};
                                    }

                                    for (auto it = gate_nums_to_prune.at(other_n).rbegin();
                                         it != gate_nums_to_prune.at(other_n).rend(); ++it)
                                    {
                                        std::swap(gate_designs[other_n][*it], gate_designs[other_n].back());

                                        gate_designs[other_n].pop_back();
                                    }

                                    std::cout << "pruned " << gate_nums_to_prune.at(other_n).size() << " from tile "
                                              << other_t << " (remaining: " << gate_designs.at(other_n).size() << ")"
                                              << std::endl;

                                    fixpoint     = false;
                                    big_fixpoint = false;
                                });
                        });
                }
            }
        }

        uint64_t pruned_total    = 0;
        uint64_t remaining_total = 0;

        std::stringstream ss{};

        ss << "\n==================================\n";
        ss << std::left << std::setw(10) << "TILE"
           << " | " << std::right << std::setw(8) << "#PRUNED"
           << " | " << std::right << std::setw(10) << "#REMAINING"
           << "\n";
        ss << "----------------------------------\n";

        gate_layout.foreach_node(
            [&](const auto& n)
            {
                if (gate_layout.is_constant(n))
                {
                    return;
                }

                const auto& designs = gate_designs.at(n);

                const uint64_t total     = gate_design_counts.at(n);
                const uint64_t pruned    = total - designs.size();
                const uint64_t remaining = designs.size();

                ss << std::left << std::setw(10) << gate_layout.get_tile(n) << " | " << std::right << std::setw(8)
                   << pruned << " | " << std::right << std::setw(10) << remaining << "\n";

                pruned_total += pruned;
                remaining_total += remaining;
            });

        ss << "----------------------+----------- +\n";
        ss << std::right << std::setw(21) << pruned_total << " | " << std::right << std::setw(10) << remaining_total
           << "\n";
        ss << "==================================\n";

        if (pruned_total > 0)
        {
            std::cout << ss.str();
        }
    }

  private:
    // forall super circuit input index, for all nodes, for all cells at node, for all other nodes, bounds on influence
    // from possible gate designs of other node (excl. skeleton)
    std::vector<std::unordered_map<
        mockturtle::node<GateLyt>,
        phmap::flat_hash_map<cell<CellLyt>, std::unordered_map<mockturtle::node<GateLyt>, std::array<double, 2>>>>>
        gate_design_influence_bounds_per_input;

    void
    tighten_gate_influence_bounds(const uint64_t input_index, const mockturtle::node<GateLyt>& n,
                                  const uint64_t cell_start_index, const uint64_t cell_end_index,
                                  std::unordered_map<mockturtle::node<GateLyt>, std::set<uint64_t>>& gate_nums_to_prune,
                                  std::mutex& mutex_to_protect_gate_nums_to_prune, bool& fixpoint) noexcept
    {
        const auto update_bounds = [&](std::array<double, 2>& current_bounds, const std::array<double, 2>& new_bounds)
        {
            if (std::isinf(current_bounds[0]))
            {
                fixpoint = false;

                current_bounds = new_bounds;

                return;
            }

            if (current_bounds[0] - new_bounds[0] < -std::numeric_limits<double>::epsilon())
            {
                fixpoint = false;

                // std::cout << "new LB: " << new_bounds[0] << std::endl;
                current_bounds[0] = new_bounds[0];
            }

            if (current_bounds[1] - new_bounds[1] > std::numeric_limits<double>::epsilon())
            {
                fixpoint = false;

                // std::cout << "new UB: " << new_bounds[1] << std::endl;
                current_bounds[1] = new_bounds[1];
            }
        };

        charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED> bounds_cds{
            canvasses_lyt, sim_params, sidb_charge_state::NONE};
        bounds_cds.initialize_matrices_for_electrostatic_calculation();  // todo optimise

        gate_layout.foreach_node(
            [&](const auto& other_n)
            {
                if (gate_layout.is_constant(other_n))
                {
                    return;
                }

                typename charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>::
                    local_external_potential_map_t& bounded_influence_from_other_canvasses =
                        bounds_cds.get_local_external_potentials_reference();

                const auto collect_bounds =
                    [&](const cell<CellLyt>& c, const mockturtle::node<GateLyt>& containing_node)
                {
                    std::array<double, 2> bounds{0, 0};

                    gate_layout.foreach_node(
                        [&](const auto& other_other_n)
                        {
                            if (gate_layout.is_constant(other_other_n) || other_n == other_other_n)
                            {
                                return;
                            }

                            const std::array<double, 2>& influence_bounds_from_other_n =
                                get_gate_design_influence_bounds(input_index, containing_node, c, other_other_n);

                            bounds[0] += influence_bounds_from_other_n[0];
                            bounds[1] += influence_bounds_from_other_n[1];
                        });

                    bounded_influence_from_other_canvasses[c] = std::move(bounds);
                };

                bounds_cds.foreach_cell(
                    [&](const auto& c)
                    { collect_bounds(c, gate_layout.get_node(bounds_cds.template get_cell_tile<tile<GateLyt>>(c))); });

                bounds_cds.update_local_external_potential();
                bounds_cds.determine_effective_charge_transition_thresholds();

                for (uint64_t j = cell_start_index; j < cell_end_index; ++j)
                {
                    const cell<CellLyt>& c = canvasses.at(n).positions.at(j);

                    // std::cout << "\ncell: " << c.x << ',' << c.y << ',' << c.z << std::endl;

                    std::array<double, 2> bounds{std::numeric_limits<double>::infinity(),
                                                 -std::numeric_limits<double>::infinity()};

                    uint64_t gate_num = 0;
                    for (const canvas_combination& gate : gate_designs.at(other_n))
                    {
                        // std::cout << "gate " << gate_num << "\t\t";  //"cells: ";
                        CellLyt canvas_of_other_n{};

                        for (const uint64_t gate_design_cell_index : gate)
                        {
                            canvas_of_other_n.assign_cell_type(
                                canvasses.at(other_n).positions.at(gate_design_cell_index),
                                sidb_technology::cell_type::NORMAL);
                        }

                        charge_distribution_surface<CellLyt> canvas_cds{canvas_of_other_n, sim_params,
                                                                        sidb_charge_state::NEGATIVE,
                                                                        cds_configuration::CHARGE_LOCATION_ONLY};
                        canvas_cds.initialize_matrices_for_electrostatic_calculation();

                        canvas_cds.add_sidb_defect_to_potential_landscape(
                            c, sidb_defect{sidb_defect_type::DB, 0, sim_params.epsilon_r, sim_params.lambda_tf});

                        bool c_is_part_of_gate = canvas_cds.get_defects().empty();

                        const auto max_index = canvas_cds.get_max_charge_index();

                        const auto is_bit_set = [&](const mockturtle::node<GateLyt>& pi_n)
                        {
                            uint64_t input_number = 0;
                            while (gate_layout.pi_at(input_number) != pi_n)
                            {
                                input_number++;
                            }
                            return (input_index & (uint64_t{1ull} << (gate_layout.num_pis() - 1 - input_number))) !=
                                   0ull;
                        };

                        bool at_least_one_charge_index_valid = false;

                        for (uint64_t charge_index = 0; charge_index <= max_index; charge_index++)
                        {
                            // std::flush(std::cout);
                            // std::cout << "charge_index: " << charge_index;

                            canvas_cds.assign_charge_index(charge_index,
                                                           charge_distribution_mode::UPDATE_CHARGE_DISTRIBUTION);

                            // todo START EXPERIMENTAL CODE

                            // enforce negatively charged output perturber
                            if (gate_layout.is_po(other_n) && canvasses.at(other_n).min_count == 1)
                            {
                                if (charge_index > 0)
                                {
                                    continue;
                                }
                            }

                            // enforce BDL wire operation
                            else if (gate_layout.is_wire(other_n))
                            {
                                const std::vector<bdl_pair<cell<CellLyt>>>& bdl_pairs = detect_bdl_pairs(
                                    canvas_of_other_n, std::nullopt,
                                    detect_bdl_pairs_params{0, detect_bdl_pairs_params{}.maximum_distance});

                                assert(bdl_pairs.size() == 1 && "there can be at most one BDL pair per wire tile");

                                if (!((canvas_cds.get_charge_state(bdl_pairs.front().upper) ==
                                           sidb_charge_state::NEUTRAL &&
                                       canvas_cds.get_charge_state(bdl_pairs.front().lower) ==
                                           sidb_charge_state::NEGATIVE) ||
                                      (canvas_cds.get_charge_state(bdl_pairs.front().upper) ==
                                           sidb_charge_state::NEGATIVE &&
                                       canvas_cds.get_charge_state(bdl_pairs.front().lower) ==
                                           sidb_charge_state::NEUTRAL)))
                                {
                                    continue;
                                }

                                if (!(gate_layout.is_pi(other_n) || gate_layout.is_po(other_n)))
                                {
                                    bool south_east = true;

                                    const tile<GateLyt>& t = gate_layout.get_tile(other_n);

                                    for (int8_t x = -2; x < 0; ++x)  // todo make parametric?
                                    {
                                        if (gate_layout.get_node(tile<GateLyt>{t.x + x, t.y + 1}) != 0)
                                        {
                                            south_east = false;
                                            break;
                                        }
                                    }

                                    const bool flip =
                                        !south_east && bdl_pairs.front().upper.y == bdl_pairs.front().lower.y;

                                    const auto is_input_bit_set = [&](const uint64_t input_number)
                                    {
                                        return (input_index &
                                                (uint64_t{1ull} << (gate_layout.num_pis() - 1 - input_number))) != 0ull;
                                    };

                                    if (flip ^ is_input_bit_set(south_east ? 0 : 1))
                                    {
                                        if (!(canvas_cds.get_charge_state(bdl_pairs.front().upper) ==
                                                  sidb_charge_state::NEUTRAL &&
                                              canvas_cds.get_charge_state(bdl_pairs.front().lower) ==
                                                  sidb_charge_state::NEGATIVE))
                                        {
                                            continue;
                                        }
                                    }
                                    else
                                    {
                                        if (!(canvas_cds.get_charge_state(bdl_pairs.front().upper) ==
                                                  sidb_charge_state::NEGATIVE &&
                                              canvas_cds.get_charge_state(bdl_pairs.front().lower) ==
                                                  sidb_charge_state::NEUTRAL))
                                        {
                                            continue;
                                        }
                                    }
                                }
                            }

                            if (gate_layout.is_pi(other_n))
                            {
                                const std::vector<bdl_pair<cell<CellLyt>>>& bdl_pairs = detect_bdl_pairs(
                                    canvas_of_other_n, std::nullopt,
                                    detect_bdl_pairs_params{0, detect_bdl_pairs_params{}.maximum_distance});

                                assert(bdl_pairs.size() == 1 && "there can be at most one BDL pair per pi tile");

                                bool is_NW_pi = false;

                                const tile<GateLyt>& t = gate_layout.get_tile(other_n);
                                for (int8_t x = 0; x < 4; ++x)  // todo make parametric?
                                {
                                    if (gate_layout.get_node(tile<GateLyt>{t.x + x, t.y + 1}) != 0)
                                    {
                                        is_NW_pi = true;

                                        break;
                                    }
                                }

                                const bool flip = !is_NW_pi && bdl_pairs.front().upper.y == bdl_pairs.front().lower.y;

                                const auto is_input_bit_set = [&](const uint64_t input_number)
                                {
                                    return (input_index &
                                            (uint64_t{1ull} << (gate_layout.num_pis() - 1 - input_number))) != 0ull;
                                };

                                if (flip ^ is_input_bit_set(is_NW_pi ? 0 : 1))
                                {
                                    if (!(canvas_cds.get_charge_state(bdl_pairs.front().upper) ==
                                              sidb_charge_state::NEUTRAL &&
                                          canvas_cds.get_charge_state(bdl_pairs.front().lower) ==
                                              sidb_charge_state::NEGATIVE))
                                    {
                                        continue;
                                    }
                                }
                                else
                                {
                                    if (!(canvas_cds.get_charge_state(bdl_pairs.front().upper) ==
                                              sidb_charge_state::NEGATIVE &&
                                          canvas_cds.get_charge_state(bdl_pairs.front().lower) ==
                                              sidb_charge_state::NEUTRAL))
                                    {
                                        continue;
                                    }
                                }
                            }

                            if (gate_layout.is_po(other_n) && canvasses.at(other_n).min_count == 2)
                            {
                                // todo HARDCODED NAND!!

                                const std::vector<bdl_pair<cell<CellLyt>>>& bdl_pairs = detect_bdl_pairs(
                                    canvas_of_other_n, std::nullopt,
                                    detect_bdl_pairs_params{0, detect_bdl_pairs_params{}.maximum_distance});

                                assert(bdl_pairs.size() == 1 && "there can be at most one BDL pair per pi tile");

                                if (get_bit(create_and_tt(), input_index))
                                {
                                    if (!(canvas_cds.get_charge_state(bdl_pairs.front().upper) ==
                                              sidb_charge_state::NEUTRAL &&
                                          canvas_cds.get_charge_state(bdl_pairs.front().lower) ==
                                              sidb_charge_state::NEGATIVE))
                                    {
                                        continue;
                                    }
                                }
                                else
                                {
                                    if (!(canvas_cds.get_charge_state(bdl_pairs.front().upper) ==
                                              sidb_charge_state::NEGATIVE &&
                                          canvas_cds.get_charge_state(bdl_pairs.front().lower) ==
                                              sidb_charge_state::NEUTRAL))
                                    {
                                        continue;
                                    }
                                }
                            }

                            // todo END EXPERIMENTAL CODE

                            canvas_cds.foreach_cell(
                                [&](const cell<CellLyt>& canvas_c)
                                {
                                    bounds_cds.assign_charge_state(canvas_c, canvas_cds.get_charge_state(canvas_c),
                                                                   charge_index_mode::KEEP_CHARGE_INDEX);
                                });

                            bounds_cds.template update_after_charge_change<true>(
                                dependent_cell_mode::FIXED, energy_calculation::KEEP_OLD_ENERGY_VALUE);

                            if (!bounds_cds.is_physically_valid())
                            {
                                // std::cout << "notvalid\t";
                                continue;
                            }

                            at_least_one_charge_index_valid = true;

                            const auto get_pot_at_c = [&]
                            {
                                if (!c_is_part_of_gate)
                                {
                                    canvas_cds.update_local_defect_potential();

                                    return *canvas_cds.get_local_defect_potential(c);
                                }

                                double pot_sum = 0;

                                canvas_cds.foreach_cell(
                                    [&](const cell<CellLyt>& canvas_c)
                                    {
                                        if (canvas_c != c)
                                        {
                                            pot_sum += canvas_cds.get_potential_between_sidbs(c, canvas_c);
                                        }
                                    });

                                return pot_sum;
                            };

                            const double pot_at_c = get_pot_at_c();

                            // std::cout << "VALID: " << fmt::format("{:.6f}", pot_at_c) << "\t";

                            bounds[0] = std::min(bounds[0], pot_at_c);
                            bounds[1] = std::max(bounds[1], pot_at_c);
                        }

                        canvas_cds.foreach_cell(
                            [&](const cell<CellLyt>& canvas_c)
                            {
                                bounds_cds.assign_charge_state(canvas_c, sidb_charge_state::NONE,
                                                               charge_index_mode::KEEP_CHARGE_INDEX);
                            });

                        if (!at_least_one_charge_index_valid)
                        {
                            const std::lock_guard guard{mutex_to_protect_gate_nums_to_prune};

                            // std::cout << "PRUNING POSSIBLE!!!!!!!!!!!!!!!!!!!!!!" << std::endl;

                            gate_nums_to_prune[other_n].emplace(gate_num);
                        }

                        ++gate_num;
                    }

                    // std::cout << '\n';
                    //
                    // std::cout << "old bounds: " <<
                    // gate_design_influence_bounds_per_input[input_index][n][c][other_n][0]
                    //           << ',' << gate_design_influence_bounds_per_input[input_index][n][c][other_n][1]
                    //           << "\tnew bounds: " << bounds[0] << ',' << bounds[1] << std::endl;

                    update_bounds(gate_design_influence_bounds_per_input[input_index][n][c][other_n], bounds);
                }
            });
    }

    void set_initial_gate_design_influence_bounds() noexcept
    {
        for (uint64_t input_index = 0; input_index < 1 << gate_layout.num_pis(); ++input_index)
        {
            gate_design_influence_bounds_per_input.emplace_back();

            // initialize all bounds to (-inf,inf)

            gate_layout.foreach_node(
                [&](const auto& n)
                {
                    if (gate_layout.is_constant(n))
                    {
                        return;
                    }

                    foreach_node<std::array<double, 2>> bounds_per_node{};

                    gate_layout.foreach_node(
                        [&](const auto& other_n)
                        {
                            if (gate_layout.is_constant(other_n))
                            {
                                return;
                            }

                            bounds_per_node[other_n] = std::array<double, 2>{-std::numeric_limits<double>::infinity(),
                                                                             std::numeric_limits<double>::infinity()};
                        });

                    for (const cell<CellLyt>& c : canvasses.at(n).positions)
                    {
                        gate_design_influence_bounds_per_input[input_index][n][c] = bounds_per_node;
                    }
                });
        }
    }
};

template <typename CellLyt, typename GateLyt>
class sidb_bdl_sub_circuit
{
  public:
    explicit sidb_bdl_sub_circuit(const sidb_bdl_circuit<CellLyt, GateLyt>& bdl_super_circuit) noexcept :
            super_circuit{bdl_super_circuit},
            nodes{get_all_nodes(super_circuit.gate_layout)},
            canvasses_lyt{super_circuit.canvasses_lyt},
            function{derive_function(super_circuit.gate_layout, nodes)},
            consistent_super_circuit_input_indices_per_input{
                make_super_circuit_input_identity(super_circuit.gate_layout.num_pis())}
    {}

    explicit sidb_bdl_sub_circuit(const sidb_bdl_circuit<CellLyt, GateLyt>&     bdl_super_circuit,
                                  const std::vector<mockturtle::node<GateLyt>>& sub_circuit_nodes) noexcept :
            super_circuit{bdl_super_circuit},
            nodes{sub_circuit_nodes},
            canvasses_lyt{sidb_bdl_circuit<CellLyt, GateLyt>::make_canvas_layout(
                super_circuit.gate_layout, super_circuit.canvasses, true, sub_circuit_nodes)},
            function{derive_function(super_circuit.gate_layout, nodes)},
            consistent_super_circuit_input_indices_per_input{
                collect_consistent_super_circuit_input_indices(super_circuit, nodes)}
    {}

    const sidb_bdl_circuit<CellLyt, GateLyt>& super_circuit{};

    const std::vector<mockturtle::node<GateLyt>> nodes{};

    const CellLyt canvasses_lyt;

    const kitty::dynamic_truth_table function{create_id_tt()};

    // [[nodiscard]] std::optional<sidb_technology::cell_type>
    // is_not_internal_perturber(const CellLyt& lyt, const cell<CellLyt>& c) const noexcept
    // {
    //     if (const auto ct = lyt.get_cell_type(c);
    //         (ct != sidb_technology::cell_type::OUTPUT_PERTURBER ||
    //          super_circuit.skeleton.get_cell_type(c) == sidb_technology::cell_type::OUTPUT_PERTURBER) &&
    //         (ct != sidb_technology::cell_type::INPUT ||
    //          super_circuit.skeleton.get_cell_type(c) == sidb_technology::cell_type::INPUT))
    //     {
    //         return ct;
    //     }
    //
    //     return std::nullopt;
    // }

    [[nodiscard]] bool input_index_possible_in_super_circuit(const uint64_t input_index) const noexcept
    {
        return !consistent_super_circuit_input_indices_per_input.at(input_index).empty();
    }
    [[nodiscard]] std::reference_wrapper<const std::vector<uint64_t>>
    get_consistent_super_circuit_input_indices(const uint64_t input_index) const noexcept
    {
        return std::cref(consistent_super_circuit_input_indices_per_input.at(input_index));
    }
    // [[nodiscard]] double get_skeleton_influence(const uint64_t       sub_circuit_input_index,
    //                                             const uint64_t       super_circuit_input_index,
    //                                             const cell<CellLyt>& c) const noexcept
    // {
    //     return *super_circuit_simulated_bdl_wires_per_input.at(sub_circuit_input_index)
    //                 .at(super_circuit_input_index)
    //                 ->get_local_internal_potential(c);
    // }

    void
    collect_influence_bounds(const CellLyt& lyt, const uint64_t input_index, const uint64_t super_circuit_input_index,
                             typename charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>::
                                 local_external_potential_map_t& influence_bounds) const noexcept
    {
        influence_bounds.reserve(lyt.num_cells());

        lyt.foreach_cell(
            [&](const auto& c)
            {
                // if (!is_not_internal_perturber(lyt, c))
                // {
                //     return;
                // }

                const mockturtle::node<GateLyt>& n = super_circuit.gate_layout.get_node(
                    super_circuit.canvasses_lyt.template get_cell_tile<tile<GateLyt>>(c));

                // const auto collect_skeleton_influence = [&]
                // {
                //     const charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>&
                //         simulated_bdl_wires =
                //             *super_circuit_simulated_bdl_wires_per_input.at(input_index).at(super_circuit_input_index);
                //
                //     assert(simulated_bdl_wires.get_local_internal_potential(c).has_value() &&
                //            "c is not part of the layout");
                //
                //     return *simulated_bdl_wires.get_local_internal_potential(c);
                // };

                // const double skeleton_influence = collect_skeleton_influence();

                std::array<double, 2> gate_design_influence_bound_sum = {0, 0};

                super_circuit.gate_layout.foreach_node(
                    [&](const auto& super_circuit_n)
                    {
                        if (super_circuit.gate_layout.is_constant(super_circuit_n) ||
                            std::find(nodes.cbegin(), nodes.cend(), super_circuit_n) != nodes.cend())
                        {
                            return;
                        }

                        const std::array<double, 2>& influence_bounds_from_super_circuit_n =
                            super_circuit.get_gate_design_influence_bounds(super_circuit_input_index, n, c,
                                                                           super_circuit_n);

                        gate_design_influence_bound_sum[0] += influence_bounds_from_super_circuit_n[0];
                        gate_design_influence_bound_sum[1] += influence_bounds_from_super_circuit_n[1];
                    });

                influence_bounds[c][0] = gate_design_influence_bound_sum[0];
                influence_bounds[c][1] = gate_design_influence_bound_sum[1];
            });
    }

    const std::vector<std::vector<uint64_t>> consistent_super_circuit_input_indices_per_input{};
    // const std::vector<
    //     std::vector<std::optional<charge_distribution_surface<CellLyt, local_external_potential_type::BOUNDED>>>>
    //     super_circuit_simulated_bdl_wires_per_input{};

  private:
    [[nodiscard]] static std::vector<mockturtle::node<GateLyt>> get_all_nodes(const GateLyt& gate_lyt) noexcept
    {
        std::vector<mockturtle::node<GateLyt>> nodes{};

        gate_lyt.foreach_node(
            [&](const auto& n)
            {
                if (!gate_lyt.is_constant(n))
                {
                    nodes.push_back(n);
                }
            });

        return nodes;
    }

    [[nodiscard]] static kitty::dynamic_truth_table
    derive_function(const GateLyt& gate_lyt, std::vector<mockturtle::node<GateLyt>> nodes) noexcept
    {
        // todo: does this function work for identity logic clumps?

        std::sort(nodes.begin(), nodes.end(),
                  [&](const mockturtle::node<GateLyt>& n1, const mockturtle::node<GateLyt>& n2)
                  {
                      const tile<GateLyt>& t1 = gate_lyt.get_tile(n1);
                      const tile<GateLyt>& t2 = gate_lyt.get_tile(n2);

                      if (t1.y == t2.y)
                      {
                          return t1.x < t2.x;
                      }

                      return t1.y < t2.y;
                  });

        // assert(gate_lyt.is_wire(nodes.back()) &&
        //        // assert(gate_lyt.is_wire(nodes.front()) && gate_lyt.is_wire(nodes.back()) &&
        //        "BDL pairs to observe logic behavior are not present");

        std::vector<std::vector<mockturtle::node<GateLyt>>> y_eqv_class_wire_nodes{};

        for (const auto& node : nodes)
        {
            if (!gate_lyt.is_wire(node))
            {
                continue;
            }

            if (!y_eqv_class_wire_nodes.empty() &&
                gate_lyt.get_tile(y_eqv_class_wire_nodes.back().back()).y == gate_lyt.get_tile(node).y)
            {
                y_eqv_class_wire_nodes.back().push_back(node);
            }
            else
            {
                y_eqv_class_wire_nodes.push_back(std::vector<mockturtle::node<GateLyt>>{node});
            }
        }

        uint64_t last_row_num_vars = y_eqv_class_wire_nodes.front().size();

        uint32_t num_computations = 0;

        for (auto row_wise_wire_nodes_it = y_eqv_class_wire_nodes.cbegin() + 1;
             row_wise_wire_nodes_it != y_eqv_class_wire_nodes.cend(); ++row_wise_wire_nodes_it)
        {
            const uint64_t this_row_num_vars = row_wise_wire_nodes_it->size();

            if (this_row_num_vars < last_row_num_vars)  // todo support multi output & misaligned wires
            {
                num_computations += last_row_num_vars - this_row_num_vars;
            }

            last_row_num_vars = this_row_num_vars;
        }

        const uint32_t num_vars = num_computations + 1;

        using logic_clump = std::vector<mockturtle::node<GateLyt>>;

        const auto is_part_of_logic_clump = [&](const logic_clump& l, const mockturtle::node<GateLyt>& n) noexcept
        {
            const tile<GateLyt>& t = gate_lyt.get_tile(n);

            for (int8_t xi = t.x == 0 ? 0 : -1; xi <= 1; ++xi)
            {
                for (int8_t yi = t.y == 0 ? 0 : -1; yi <= 1; ++yi)
                {
                    if (const mockturtle::node<GateLyt>& other_n = gate_lyt.get_node({t.x + xi, t.y + yi});
                        other_n != 0 && other_n != n && std::find(l.cbegin(), l.cend(), other_n) != l.cend())
                    {
                        return true;
                    }
                }
            }

            return false;
        };

        std::vector<std::vector<logic_clump>> y_eqv_class_logic_clumps{};

        // todo this might fail in the following situation:
        // when two horizontally aligned logic clumps are such that the rightmost one has a tile that is seen first
        // (ie., lower y value than any of the tiles in the leftmost logic clump). the logic clumps need to be sorted in
        // their equivalence class

        for (const auto& node : nodes)
        {
            if (gate_lyt.is_wire(node))
            {
                continue;
            }

            if (y_eqv_class_logic_clumps.empty())
            {
                y_eqv_class_logic_clumps.push_back(std::vector<std::vector<mockturtle::node<GateLyt>>>{{node}});

                continue;
            }

            bool is_part_of_existing_logic_clump = false;

            for (logic_clump& l : y_eqv_class_logic_clumps.back())
            {
                if (is_part_of_logic_clump(l, node))
                {
                    l.emplace_back(node);  // add to logic clump

                    is_part_of_existing_logic_clump = true;

                    break;
                }
            }

            if (is_part_of_existing_logic_clump)
            {
                continue;
            }

            const uint64_t min_y  = gate_lyt.get_tile(y_eqv_class_logic_clumps.back().front().front()).y;  // todo check
            const uint64_t this_y = gate_lyt.get_tile(node).y;

            bool is_in_equivalent_y_class = true;

            if (min_y != this_y)
            {
                for (const std::vector<mockturtle::node<GateLyt>>& wire_eqv_class : y_eqv_class_wire_nodes)
                {
                    if (gate_lyt.get_tile(wire_eqv_class.front()).y <= min_y)
                    {
                        continue;
                    }

                    if (gate_lyt.get_tile(wire_eqv_class.front()).y > this_y)
                    {
                        break;
                    }

                    // BDL pair exists in between last logic clump and this node -- not in equivalent y class

                    is_in_equivalent_y_class = false;

                    break;
                }
            }

            if (is_in_equivalent_y_class)
            {
                y_eqv_class_logic_clumps.back().push_back(std::vector<mockturtle::node<GateLyt>>{node});
            }
            else
            {
                y_eqv_class_logic_clumps.push_back(std::vector<std::vector<mockturtle::node<GateLyt>>>{{node}});
            }
        }

        assert(std::all_of(y_eqv_class_logic_clumps.cbegin(), y_eqv_class_logic_clumps.cend(),
                           [&](const std::vector<logic_clump>& lv)
                           {
                               return std::all_of(
                                   lv.cbegin(), lv.cend(),
                                   [&](const logic_clump& l)
                                   {
                                       return std::all_of(
                                           l.cbegin(), l.cend(), [&](const mockturtle::node<GateLyt>& n)
                                           { return gate_lyt.node_function(l.front()) == gate_lyt.node_function(n); });
                                   });
                           }) &&
               "there exists a logic clump that is heterogeneous");

        using tt = kitty::dynamic_truth_table;

        // kitty::dynamic_truth_table tt{num_computations + 1};

        uint32_t logic_fan_in = y_eqv_class_logic_clumps.empty() ? num_vars : 0;

        if (!y_eqv_class_logic_clumps.empty())
        {
            for (const logic_clump& l : y_eqv_class_logic_clumps.front())
            {
                logic_fan_in += gate_lyt.node_function(l.front()).num_vars();
            }
        }

        std::vector<tt> vars{};

        for (uint8_t i = 0; i < logic_fan_in; ++i)
        {
            tt var{logic_fan_in};
            kitty::create_nth_var(var, i);
            // std::cout << "var: ";
            // kitty::print_binary(var);
            // std::cout << std::endl;
            vars.push_back(std::move(var));
        }

        std::reverse(vars.begin(), vars.end());

        const auto print_vars = [&]
        {
            // std::cout << "VARS:" << std::endl;
            // for (const auto& var : vars)
            // {
            //     kitty::print_binary(var);
            //     std::cout << '\t';
            // }
            // std::cout << std::endl;
        };

        for (const std::vector<logic_clump>& lv : y_eqv_class_logic_clumps)
        {
            print_vars();

            uint32_t var_offset = 0;

            for (const logic_clump& l : lv)
            {
                const tt& this_tt = gate_lyt.node_function(l.front());

                // std::cout << "this tt: ";
                // kitty::print_binary(this_tt);
                // std::cout << std::endl;

                std::vector<tt> input_vars{};
                for (uint32_t i = 0; i < this_tt.num_vars(); ++i)
                {
                    input_vars.emplace_back(vars.at(var_offset + i));
                }

                vars[var_offset] = kitty::compose_truth_table<tt, tt>(this_tt, input_vars);

                for (uint32_t i = 0; i < this_tt.num_vars() - 1; ++i)
                {
                    vars.erase(vars.begin() + var_offset + 1 + i);  // erase consumed variable(s)
                }

                var_offset++;
            }
        }

        print_vars();

        assert(vars.size() == 1 && "only single output is supported");

        return vars.front();
    }

    [[nodiscard]] static std::vector<std::vector<uint64_t>>
    make_super_circuit_input_identity(const uint64_t num_inputs) noexcept
    {
        std::vector<std::vector<uint64_t>> tiles{};
        tiles.reserve(num_inputs);

        for (uint64_t i = 0; i < 1 << num_inputs; ++i)
        {
            tiles.push_back({i});
        }

        return tiles;
    }

    [[nodiscard]] static std::vector<std::vector<uint64_t>>
    collect_consistent_super_circuit_input_indices(const sidb_bdl_circuit<CellLyt, GateLyt>&     super_circuit,
                                                   const std::vector<mockturtle::node<GateLyt>>& nodes) noexcept
    {
        uint64_t num_sub_circuit_inputs = 0;
        for (const auto& n : nodes)
        {
            if (super_circuit.gate_layout.is_pi(n))
            {
                num_sub_circuit_inputs++;

                // continue;
            }

            // // count also wire segment that is connected to a pi that is not included in the nodes
            // const tile<GateLyt>& t = super_circuit.gate_layout.get_tile(n);
            //
            // for (int8_t x = -3; x < 3; ++x)  // todo make parametric?
            // {
            //     // todo works only with at most one non-PI wire segment
            //     if (const tile<GateLyt> pi_t{t.x + x, t.y - 1}; super_circuit.gate_layout.is_pi_tile(pi_t))
            //     {
            //         if (std::find(nodes.cbegin(), nodes.cend(), super_circuit.gate_layout.get_node(pi_t)) ==
            //         nodes.cend())
            //         {
            //             num_sub_circuit_inputs++;
            //         }
            //
            //         break;
            //     }
            // }
        }

        std::vector<std::vector<uint64_t>> consistent_super_circuit_input_indices_per_sub_circuit_index{};
        consistent_super_circuit_input_indices_per_sub_circuit_index.reserve(1 << num_sub_circuit_inputs);

        for (uint64_t sub_circuit_iix = 0; sub_circuit_iix < 1 << num_sub_circuit_inputs; ++sub_circuit_iix)
        {
            std::vector<uint64_t> consistent_super_circuit_input_indices{};
            consistent_super_circuit_input_indices.reserve(1 << super_circuit.gate_layout.num_pis());

            for (uint64_t super_circuit_iix = 0; super_circuit_iix < 1 << super_circuit.gate_layout.num_pis();
                 ++super_circuit_iix)
            {
                // todo this is a very crude implementation, works (?) only for a single gate

                if (num_sub_circuit_inputs == 0)
                {
                    consistent_super_circuit_input_indices.push_back(super_circuit_iix);

                    continue;
                }

                if (num_sub_circuit_inputs == 2)
                {
                    if (sub_circuit_iix == super_circuit_iix)
                    {
                        consistent_super_circuit_input_indices.push_back(super_circuit_iix);
                    }

                    continue;
                }

                const auto is_input_bit_set = [&](const uint64_t input_number)
                {
                    return (super_circuit_iix &
                            (uint64_t{1ull} << (super_circuit.gate_layout.num_pis() - 1 - input_number))) != 0ull;
                };

                bool has_NW_pi = false;

                for (const auto& n : nodes)
                {
                    if (super_circuit.gate_layout.is_pi(n))
                    {
                        const tile<GateLyt>& t = super_circuit.gate_layout.get_tile(n);
                        for (int8_t x = 0; x < 4; ++x)  // todo make parametric?
                        {
                            if (super_circuit.gate_layout.get_node(tile<GateLyt>{t.x + x, t.y + 1}) != 0)
                            {
                                has_NW_pi = true;

                                break;
                            }
                        }
                    }
                }

                if (is_input_bit_set(has_NW_pi ? 0 : 1) == sub_circuit_iix)
                {
                    consistent_super_circuit_input_indices.push_back(super_circuit_iix);
                }
            }

            consistent_super_circuit_input_indices_per_sub_circuit_index.push_back(
                std::move(consistent_super_circuit_input_indices));
        }

        return consistent_super_circuit_input_indices_per_sub_circuit_index;
    }
};

template <typename CellLyt, typename GateLyt>
struct sidb_cell_level_bdl_circuit
{
    /**
     * SiDB cell-level layout.
     */
    const CellLyt&                                cell_layout;
    const sidb_bdl_sub_circuit<CellLyt, GateLyt>& circuit{};

    explicit sidb_cell_level_bdl_circuit(const CellLyt&                                lyt,
                                         const sidb_bdl_sub_circuit<CellLyt, GateLyt>& bdl_circuit) noexcept :
            cell_layout{lyt},
            circuit{bdl_circuit}
    {}
};

}  // namespace fiction

#endif  // SIDB_BDL_CIRCUIT_HPP

//
// Created by marcel on 05.04.23.
//

#ifndef FICTION_SIDB_SIMULATION_RESULT_HPP
#define FICTION_SIDB_SIMULATION_RESULT_HPP

#include "fiction/algorithms/simulation/sidb/minimum_energy.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/technology/constants.hpp"

#include <any>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fiction
{

/**
 * This struct defines a unified return type for all SiDB simulation algorithms. It contains the name of the algorithm,
 * the total simulation runtime, the charge distributions determined by the algorithm, the physical parameters used in
 * the simulation, and (optional) algorithm-specific named simulation parameters.
 *
 * @tparam Lyt SiDB cell-level layout type.
 * @tparam ExtPotType Type of external local potential (single-valued / bounded).
 */
template <typename Lyt, local_external_potential_type ExtPotType = local_external_potential_type::SINGLE_VALUED>
struct sidb_simulation_result
{
    /**
     * Default constructor. It only exists to allow for the use of `static_assert` statements that restrict the type of
     * `Lyt`.
     */
    sidb_simulation_result() noexcept
    {
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
    }
    /**
     * Name of the algorithm used to determine the charge distributions.
     */
    std::string algorithm_name{};
    /**
     * Total simulation runtime in seconds.
     */
    std::chrono::duration<double> simulation_runtime{};
    /**
     * Charge distributions determined by the algorithm.
     */
    std::vector<charge_distribution_surface<Lyt, ExtPotType>> charge_distributions{};
    /**
     * Physical parameters used in the simulation.
     */
    sidb_simulation_parameters simulation_parameters{};
    /**
     * Additional named simulation parameters. This is used to store algorithm-dependent parameters that are not part of
     * the `sidb_simulation_parameters` struct.
     *
     * The key of the map is the name of the parameter, the element is the value of the parameter.
     */
    std::unordered_map<std::string, std::any> additional_simulation_parameters{};
    /**
     * This function computes the ground state of the charge distributions.
     *
     * @note If degenerate states exist in the simulation result, this function will return multiple ground states that
     * all possess the same system energy. todo
     *
     * @return A vector of charge distributions with the minimal energy.
     */
    [[nodiscard]] std::vector<charge_distribution_surface<Lyt, ExtPotType>> groundstates() const noexcept
    {
        static_assert(ExtPotType != local_external_potential_type::BOUNDED, "ExtPotType cannot be BOUNDED.");

        std::vector<charge_distribution_surface<Lyt, ExtPotType>> groundstate_charge_distributions{};

        std::set<uint64_t> charge_indices{};

        // Find all unique charge indices. This is done because simulation results can have multiple identical
        // charge distributions.
        for (auto& cds : charge_distributions)
        {
            cds.charge_distribution_to_index_general();
            charge_indices.insert(cds.get_charge_index_and_base().first);
        }

        // Find the minimum energy
        double min_energy = minimum_energy(charge_distributions.cbegin(), charge_distributions.cend());

        for (const auto charge_index : charge_indices)
        {
            const auto cds_it = std::find_if(charge_distributions.cbegin(), charge_distributions.cend(),
                                             [&](const auto& cds)
                                             {
                                                 return cds.get_charge_index_and_base().first == charge_index &&
                                                        std::abs(cds.get_electrostatic_potential_energy() -
                                                                 min_energy) < constants::ERROR_MARGIN;
                                             });

            if (cds_it != charge_distributions.cend())
            {
                groundstate_charge_distributions.push_back(std::move(*cds_it));
            }
        }

        return groundstate_charge_distributions;
    }

    void reduce_to_groundstates_under_bounded_energy() noexcept
    {
        static_assert(ExtPotType == local_external_potential_type::BOUNDED, "ExtPotType needs to be BOUNDED.");

        // Sort by energy[0] ascending, then energy[1] descending
        std::sort(charge_distributions.begin(), charge_distributions.end(),
                  [](const auto& a, const auto& b)
                  {
                      const auto e1 = a.get_electrostatic_potential_energy();
                      const auto e2 = b.get_electrostatic_potential_energy();

                      if (std::abs(e1[0] - e2[0]) < constants::ERROR_MARGIN)
                      {
                          return e1[1] > e2[1];
                      }

                      return e1[0] < e2[0];
                  });

        double max_energy = charge_distributions.front().get_electrostatic_potential_energy()[1];

        for (uint64_t ix = 1; ix < charge_distributions.size(); ++ix)
        {
            const auto& energies = charge_distributions[ix].get_electrostatic_potential_energy();

            if (energies[0] > max_energy - constants::ERROR_MARGIN)
            {
                charge_distributions.resize(ix);

                break;  // Done: all future entries start after max_energy
            }

            max_energy = std::max(max_energy, energies[1]);
        }
    }

    std::vector<double> get_ordered_weights_under_bounded_energy(const double alpha = 1.0) noexcept
    {
        static_assert(ExtPotType == local_external_potential_type::BOUNDED, "ExtPotType needs to be BOUNDED.");

        if (charge_distributions.empty())
        {
            return {};
        }

        double min_E_low  = charge_distributions.front().get_electrostatic_potential_energy()[0];
        double min_E_high = charge_distributions.front().get_electrostatic_potential_energy()[1];
        double max_E_high = min_E_high;

        // Determine min and max high energies
        for (const auto& cds : charge_distributions)
        {
            const auto energies = cds.get_electrostatic_potential_energy();
            min_E_low           = std::min(min_E_low, energies[0]);
            min_E_high          = std::min(min_E_high, energies[1]);
            max_E_high          = std::max(max_E_high, energies[1]);
        }

        constexpr double    epsilon = 1e-12;
        std::vector<double> weights;

        for (const auto& cds : charge_distributions)
        {
            const auto energies = cds.get_electrostatic_potential_energy();

            const double decay     = std::exp(-alpha * (energies[0] - min_E_low));
            const double certainty = 1.0 - (max_E_high - energies[1]) / (max_E_high - min_E_high + epsilon);

            weights.push_back(decay * certainty);
        }

        // Normalize weights
        const double total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);
        if (total_weight > epsilon)
        {
            for (auto& w : weights)
            {
                w /= total_weight;
            }
        }

        return weights;
    }

    // std::vector<double> get_ordered_weights_under_bounded_energy() noexcept
    // {
    //     static_assert(ExtPotType == local_external_potential_type::BOUNDED, "ExtPotType needs to be BOUNDED.");
    //
    //     if (charge_distributions.empty())
    //     {
    //         return {};
    //     }
    //
    //     // // Ensure sorting order matches reduce_to_groundstates_under_bounded_energy
    //     // std::sort(charge_distributions.begin(), charge_distributions.end(),
    //     //           [](const auto& a, const auto& b)
    //     //           {
    //     //               const auto e1 = a.get_electrostatic_potential_energy();
    //     //               const auto e2 = b.get_electrostatic_potential_energy();
    //     //
    //     //               if (std::abs(e1[0] - e2[0]) < constants::ERROR_MARGIN)
    //     //               {
    //     //                   return e1[1] > e2[1];
    //     //               }
    //     //
    //     //               return e1[0] < e2[0];
    //     //           });
    //
    //     const std::size_t   n = charge_distributions.size();
    //     std::vector<double> raw_weights(n, 0.0);
    //
    //     const auto   base_energy = charge_distributions.front().get_electrostatic_potential_energy();
    //     const double lower_bound = base_energy[0];
    //     const double upper_bound = base_energy[1];
    //
    //     const double bound_width = std::max(upper_bound - lower_bound, constants::ERROR_MARGIN);  // avoid div by 0
    //
    //     for (std::size_t i = 0; i < n; ++i)
    //     {
    //         const auto energy = charge_distributions[i].get_electrostatic_potential_energy();
    //
    //         // Distance from lower bound, normalized by bound width
    //         const double relative_energy_distance = (energy[0] - lower_bound) / bound_width;
    //
    //         // Weight: inverse exponential based on distance to lower bound (more energy = less weight)
    //         // Modified by number of configurations: more configurations flatten the curve
    //         raw_weights[i] = std::exp(-relative_energy_distance * static_cast<double>(n));
    //     }
    //
    //     // Normalize weights to sum to 1
    //     const double weight_sum = std::accumulate(raw_weights.begin(), raw_weights.end(), 0.0);
    //     // std::cout <<"\nweights:" << weight_sum << std::endl;
    //     if (weight_sum > 0.0)
    //     {
    //         for (auto& w : raw_weights)
    //         {
    //             w /= weight_sum;
    //         }
    //     }
    //
    //     return raw_weights;
    // }
};

}  // namespace fiction

#endif  // FICTION_SIDB_SIMULATION_RESULT_HPP

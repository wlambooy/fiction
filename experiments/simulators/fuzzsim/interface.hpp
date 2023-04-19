//
// Created by Jan Drewniok on 31.01.23.
//

#ifndef FUZZSIM_SIQAD_PLUGIN_INTERFACE_HPP
#define FUZZSIM_SIQAD_PLUGIN_INTERFACE_HPP

#include "simulators/logger.hpp"
#include "simulators/siqadconn.cc"
#include "simulators/siqadconn.h"

#include "fiction/algorithms/simulation/sidb/fuzzsim.hpp"
#include "fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp"
#include "fiction/layouts/cell_level_layout.hpp"
#include "fiction/technology/charge_distribution_surface.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"

#include <fmt/format.h>

#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

class fuzzsim_interface
{
  public:
    //! Constructor for QuickSimInterface
    fuzzsim_interface(const std::string_view& t_in_path, const std::string_view& t_out_path, const bool verbose = true,
                      const int log_l = logger::MSG) :
            log_level{log_l},
            in_path{t_in_path},
            out_path{t_out_path}
    {
        sqconn = std::make_unique<SiQADConnector>("FuzzSim", static_cast<std::string>(in_path),
                                                  static_cast<std::string>(out_path), verbose);
        initialize_fiction_layout();
    }

    int run_simulation()
    {
        while (sim_results.valid_lyts.empty())
        {
            fiction::fuzzsim<fiction::sidb_cell_clk_lyt_siqad>(layout, sim_par, &sim_results);
        }

        return EXIT_SUCCESS;
    }

    void write_sim_results()
    {
        sim_results.report();

        // create the vector of strings for the db locations
        const auto data = sim_results.valid_lyts.cbegin()->get_all_sidb_location_in_nm();

        std::vector<std::pair<std::string, std::string>> dbl_data{};
        dbl_data.reserve(data.size());

        for (const auto& d : data)
        {
            dbl_data.emplace_back(std::to_string(d.first * 10E9), std::to_string(d.second * 10E9));
        }

        sqconn->setExport("db_loc", dbl_data);

        std::vector<std::vector<std::string>> db_dist_data{};
        db_dist_data.reserve(sim_results.valid_lyts.size());

        for (const auto& lyt : sim_results.valid_lyts)
        {
             db_dist_data.push_back({{
                fiction::charge_configuration_to_string(lyt.get_all_sidb_charges()),  // config
                std::to_string(lyt.get_system_energy()),                              // energy
                std::to_string(1),                                                    // occurrence freq
                "1",                                                                  // metastability
                "3"  // 3-state (GUI still does not work for 2)
            }});
        }

        sqconn->setExport("db_charge", db_dist_data);

        sqconn->writeResultsXml();
    }

    [[nodiscard]] fiction::fuzzsim_params& get_fuzzsim_params() noexcept
    {
        return sim_par;
    }

    [[nodiscard]] fiction::fuzzsim_stats<fiction::sidb_cell_clk_lyt_siqad> get_simulation_results() const noexcept
    {
        return sim_results;
    }

    [[nodiscard]] uint64_t get_auto_fail() const
    {
        return auto_fail;
    }

    [[nodiscard]] uint64_t get_cell_num() const
    {
        return layout.num_cells();
    }

  private:
    void initialize_fiction_layout()
    {
        logger log(log_level);

        // grab all physical locations
        log.debug() << "Grab all physical locations..." << std::endl;

        auto* const db_collection = sqconn->dbCollection();

        std::vector<std::pair<double, double>> db_locs{};
        db_locs.reserve(db_collection->db_tree_inner->size());

        for (const auto& db : *db_collection)
        {
            db_locs.push_back(fiction::sidb_nm_position<fiction::sidb_cell_clk_lyt_siqad>(
                fiction::sidb_simulation_parameters{}, {db->n, db->m, db->l}));

            layout.assign_cell_type({db->n, db->m, db->l}, fiction::sidb_cell_clk_lyt_siqad::cell_type::NORMAL);

            log.debug() << fmt::format("DB loc: x={}, y={}", db_locs.back().first, db_locs.back().second) << std::endl;
        }

        fiction::sidb_simulation_parameters params{2};

        try
        {
            // variables: physical
            params.mu = std::stod(sqconn->getParameter("muzm"));

            params.epsilon_r = std::round(std::stod(sqconn->getParameter("eps_r")) * 100) / 100;  // round to two digits
            params.lambda_tf = std::stod(sqconn->getParameter("debye")) * 1E-9;

            auto_fail = std::stoi(sqconn->getParameter("autofail"));

            // variables: QuickSim
            const auto iteration_steps = static_cast<uint64_t>(std::stoi(sqconn->getParameter("iteration_steps")));
            const auto alpha           = std::stod(sqconn->getParameter("alpha"));

            // variables: FuzzSim
            const auto trials = static_cast<uint64_t>(std::stoi(sqconn->getParameter("quicksim_trials")));
            const auto max_perm_log = static_cast<uint64_t>(std::stoi(sqconn->getParameter("maximum_permutations_logarithm")));
            const auto max_charge_layouts = static_cast<uint64_t>(std::stoi(sqconn->getParameter("maximum_charge_layouts")));

            sim_par = fiction::fuzzsim_params{fiction::quicksim_params{params, iteration_steps, alpha},
                                              trials, max_perm_log, max_charge_layouts};

            // prevent number of threads to be negative
            if (const auto number_threads = static_cast<int64_t>(std::stod(sqconn->getParameter("num_threads")));
                number_threads >= 0)
            {
                // Update sim_par with number_threads
                sim_par.qs_params.number_threads = static_cast<uint64_t>(number_threads);
            }

            log.echo() << "Retrieval from SiQADConn complete." << std::endl;
        }
        catch (...)
        {
            log.critical() << "Could not retrieve all parameters from SiQADConn." << std::endl;
            throw;
        }
    }

    // Instances
    std::unique_ptr<SiQADConnector> sqconn = nullptr;

    // variables
    uint64_t                                                  auto_fail;
    const int                                                 log_level;
    const std::string_view                                    in_path;
    const std::string_view                                    out_path;
    fiction::sidb_cell_clk_lyt_siqad                          layout{};
    fiction::fuzzsim_params                                   sim_par{};
    fiction::fuzzsim_stats<fiction::sidb_cell_clk_lyt_siqad> sim_results{};
};

#endif  // FUZZSIM_SIQAD_PLUGIN_INTERFACE_HPP
//
// Created by Jan Drewniok on 31.01.23.
//

#include "simulators/fuzzsim/interface.hpp"
#include "simulators/logger.hpp"
#include "simulators/timer.hpp"

#include <fmt/format.h>

#include <cstdlib>
#include <iostream>
#include <string_view>
#include <vector>

int main(int argc, char* argv[])
{
    std::cout << "Physeng invoked" << std::endl;

    const std::vector<std::string_view> cml_args(argv, argv + argc);

    std::cout << "*** Argument Parsing ***" << std::endl;

    if (argc < 3)
    {
        std::cerr << "Error: Too few arguments. Expected input and output file names.\n";

        return EXIT_FAILURE;
    }

    const std::string_view if_name = cml_args[1];
    const std::string_view of_name = cml_args[2];

    bool verbose   = false;
    auto log_level = logger::DBG;

    for (auto it = cml_args.cbegin() + 3; it != cml_args.cend(); ++it)
    {
        if (*it == "--debug")
        {
            // show additional debug information
            std::cout << "--debug: Showing additional outputs." << std::endl;
            verbose   = true;
            log_level = logger::DBG;
        }
        else
        {
            std::cout << fmt::format("Unrecognized command-line argument: {}", *it) << std::endl;
        }
    }

    logger log{log_level};

    simglobal::stopwatch stopwatch;

    log.echo() << "In File: " << if_name << std::endl;
    log.echo() << "Out File: " << of_name << std::endl;

    log.echo() << "\n*** Initiate FuzzSim interface ***" << std::endl;
    log.echo() << "\n*** Read Simulation parameters ***" << std::endl;
    fuzzsim_interface fuzzsim_interface{if_name, of_name, verbose};

    if (fuzzsim_interface.get_auto_fail() < fuzzsim_interface.get_cell_num())
    {
        log.warning() << "Problem size > autofail threshold, exiting." << std::endl;
        return EXIT_FAILURE;
    }

    log.echo() << "\n*** Invoke simulation ***" << std::endl;
    stopwatch.start();
    fuzzsim_interface.run_simulation();
    stopwatch.end();

    log.echo() << "\n*** Write simulation results ***" << std::endl;
    fuzzsim_interface.write_sim_results();

    log.echo() << "\n*** FuzzSim Complete ***" << std::endl;

    stopwatch.print_stopwatch(log_level);

    return EXIT_SUCCESS;
}
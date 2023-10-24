//
// Created by marcel on 22.05.23.
//

#ifndef PYFICTION_CRITICAL_TEMPERATURE_HPP
#define PYFICTION_CRITICAL_TEMPERATURE_HPP

#include "pyfiction/documentation.hpp"
#include "pyfiction/types.hpp"

#include <fiction/algorithms/simulation/sidb/critical_temperature.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace pyfiction
{

namespace detail
{

template <typename Lyt>
void critical_temperature(pybind11::module& m)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    /**
     * Critical temperature statistics.
     */
    py::class_<fiction::critical_temperature_stats<Lyt>>(m, "critical_temperature_stats",
                                                         DOC(fiction_critical_temperature_stats))
        .def(py::init<>())
        .def_readwrite("algorithm_name", &fiction::critical_temperature_stats<Lyt>::algorithm_name,
                       DOC(fiction_critical_temperature_stats_algorithm_name))
        .def_readwrite("critical_temperature", &fiction::critical_temperature_stats<Lyt>::critical_temperature,
                       DOC(fiction_critical_temperature_stats_critical_temperature))
        .def_readwrite("num_valid_lyt", &fiction::critical_temperature_stats<Lyt>::num_valid_lyt,
                       DOC(fiction_critical_temperature_stats_num_valid_lyt))
        .def_readwrite("energy_between_ground_state_and_first_erroneous",
                       &fiction::critical_temperature_stats<Lyt>::energy_between_ground_state_and_first_erroneous,
                       DOC(fiction_critical_temperature_stats_energy_between_ground_state_and_first_erroneous))
        .def("report", &fiction::critical_temperature_stats<Lyt>::report,
             DOC(fiction_critical_temperature_stats_report))

        ;

    m.def("critical_temperature_gate_based", &fiction::critical_temperature_gate_based<Lyt, py_tt>, "lyt"_a, "spec"_a,
          "params"_a = fiction::critical_temperature_params{}, "stats"_a = nullptr, DOC(fiction_critical_temperature));
}

}  // namespace detail

inline void critical_temperature(pybind11::module& m)
{
    namespace py = pybind11;

    /**
     * Critical temperature mode.
     */
    py::enum_<fiction::critical_temperature_params::critical_temperature_mode>(
        m, "critical_temperature_mode", DOC(fiction_critical_temperature_params_critical_temperature_mode))
        .value("GATE_BASED_SIMULATION",
               fiction::critical_temperature_params::critical_temperature_mode::GATE_BASED_SIMULATION,
               DOC(fiction_critical_temperature_params_critical_temperature_mode_GATE_BASED_SIMULATION))
        .value("NON_GATE_BASED_SIMULATION",
               fiction::critical_temperature_params::critical_temperature_mode::NON_GATE_BASED_SIMULATION,
               DOC(fiction_critical_temperature_params_critical_temperature_mode_NON_GATE_BASED_SIMULATION))

        ;

    /**
     * Simulation engine.
     */
    py::enum_<fiction::critical_temperature_params::simulation_engine>(
        m, "simulation_engine", DOC(fiction_critical_temperature_params_simulation_engine))
        .value("EXACT", fiction::critical_temperature_params::simulation_engine::EXACT,
               DOC(fiction_critical_temperature_params_simulation_engine_EXACT))
        .value("APPROXIMATE", fiction::critical_temperature_params::simulation_engine::APPROXIMATE,
               DOC(fiction_critical_temperature_params_simulation_engine_APPROXIMATE))

        ;

    /**
     * Critical temperature parameters.
     */
    py::class_<fiction::critical_temperature_params>(m, "critical_temperature_params",
                                                     DOC(fiction_critical_temperature_params))
        .def(py::init<>())
        .def_readwrite("simulation_params", &fiction::critical_temperature_params::simulation_params,
                       DOC(fiction_critical_temperature_params_simulation_params))
        .def_readwrite("engine", &fiction::critical_temperature_params::engine,
                       DOC(fiction_critical_temperature_params_engine))
        .def_readwrite("confidence_level", &fiction::critical_temperature_params::confidence_level,
                       DOC(fiction_critical_temperature_params_confidence_level))
        .def_readwrite("max_temperature", &fiction::critical_temperature_params::max_temperature,
                       DOC(fiction_critical_temperature_params_max_temperature))
        .def_readwrite("bdl_params", &fiction::critical_temperature_params::bdl_params,
                       DOC(fiction_critical_temperature_params_bdl_params))

        ;

    // NOTE be careful with the order of the following calls! Python will resolve the first matching overload!

    detail::critical_temperature<py_charge_distribution_surface>(m);
}

}  // namespace pyfiction

#endif  // PYFICTION_CRITICAL_TEMPERATURE_HPP

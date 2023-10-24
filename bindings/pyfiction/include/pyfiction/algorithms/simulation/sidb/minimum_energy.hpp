//
// Created by marcel on 22.05.23.
//

#ifndef PYFICTION_MINIMUM_ENERGY_HPP
#define PYFICTION_MINIMUM_ENERGY_HPP

#include "pyfiction/documentation.hpp"
#include "pyfiction/types.hpp"

#include <fiction/algorithms/simulation/sidb/minimum_energy.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace pyfiction
{

namespace detail
{

template <typename Lyt>
void minimum_energy(pybind11::module& m)
{
    using namespace pybind11::literals;

    m.def("minimum_energy", &fiction::minimum_energy<Lyt>, "charge_lyts"_a, DOC(fiction_minimum_energy));
}

}  // namespace detail

inline void minimum_energy(pybind11::module& m)
{
    // NOTE be careful with the order of the following calls! Python will resolve the first matching overload!

    detail::minimum_energy<py_charge_distribution_surface>(m);
}

}  // namespace pyfiction

#endif  // PYFICTION_MINIMUM_ENERGY_HPP

//
// Created by marcel on 25.04.23.
//

#ifndef PYFICTION_READ_SQD_LAYOUT_HPP
#define PYFICTION_READ_SQD_LAYOUT_HPP

#include "pyfiction/documentation.hpp"
#include "pyfiction/types.hpp"

#include <fiction/io/read_sqd_layout.hpp>

#include <pybind11/pybind11.h>

#include <string_view>

namespace pyfiction
{

namespace detail
{

template <typename Lyt>
void read_sqd_layout(pybind11::module& m)
{
    using namespace pybind11::literals;

    Lyt (*read_sqd_layout_function_pointer)(const std::string_view&, const std::string_view&) =
        &fiction::read_sqd_layout<Lyt>;

    m.def("read_sqd_layout", read_sqd_layout_function_pointer, "filename"_a, "layout_name"_a = "",
          DOC(fiction_read_sqd_layout_3));
}

}  // namespace detail

inline void read_sqd_layout(pybind11::module& m)
{
    detail::read_sqd_layout<py_sidb_layout>(m);
}

}  // namespace pyfiction

#endif  // PYFICTION_READ_SQD_LAYOUT_HPP

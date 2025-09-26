//
// Created by Willem Lambooy on 13/08/2025.
//

#ifndef SIDB_BDL_SKELETONS_HPP
#define SIDB_BDL_SKELETONS_HPP

#include <fiction/io/read_sqd_layout.hpp>
#include <fiction/technology/cell_technologies.hpp>
#include <fiction/technology/sidb_skeleton_gate_library.hpp>
#include <fiction/technology/sidb_skeletons/sidb_skeleton.hpp>

#include <string_view>

using namespace fiction;

class sidb_bdl_skeleton_1 : public sidb_skeleton_gate_library<27, 30>  // width and height of a hexagon
{
  public:
    sidb_bdl_skeleton_1() noexcept :
            sidb_skeleton_gate_library{read_sqd_layout<sidb_skeleton_t>(
                "/home/willem/fiction/include/fiction/technology/sidb_skeletons/new_mini_bestagon.sqd", "skeleton")}
    {}
};

class sidb_bdl_skeleton_august : public sidb_skeleton_gate_library<42, 27>  // width and height of a hexagon
{
  public:
    sidb_bdl_skeleton_august() noexcept :
            sidb_skeleton_gate_library{read_sqd_layout<sidb_skeleton_t>(
                "/home/willem/fiction/include/fiction/technology/sidb_skeletons/new_mini_bestagon_august.sqd",
                "skeleton")}
    {}
};

class sidb_bdl_skeleton_original_bestagon : public sidb_skeleton_gate_library<60, 46>  // width and height of a hexagon
{
  public:
    sidb_bdl_skeleton_original_bestagon() noexcept :
            sidb_skeleton_gate_library{read_sqd_layout<sidb_skeleton_t>(
                "/home/willem/fiction/include/fiction/technology/sidb_skeletons/original_bestagon.sqd", "skeleton")}
    {}
};

#endif  // SIDB_BDL_SKELETONS_HPP

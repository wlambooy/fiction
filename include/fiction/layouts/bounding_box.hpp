//
// Created by marcel on 13.01.22.
//

#ifndef FICTION_BOUNDING_BOX_HPP
#define FICTION_BOUNDING_BOX_HPP

#include "fiction/layouts/cell_level_layout.hpp"
#include "fiction/layouts/coordinates.hpp"
#include "fiction/technology/cell_ports.hpp"
#include "fiction/traits.hpp"
#include "fiction/types.hpp"
#include "fiction/utils/layout_utils.hpp"

// data types cannot properly be converted to bit field types
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"

namespace fiction
{

/**
 * A 2D bounding box object that computes a minimum-sized box around all non-empty coordinates in a given layout.
 * Layouts can be of arbitrary size and, thus, may be larger than their contained elements. Sometimes, it might be
 * necessary to know exactly which space the associated layout internals occupy. A bounding box computes coordinates
 * that span a minimum-sized rectangle that encloses all non-empty layout coordinates.
 *
 * The bounding box does not automatically updated when the layout changes. Call `update_bounding_box()` to recompute
 * it.
 *
 * @tparam Lyt Gate-level or cell-level layout type.
 */
template <typename Lyt>
class bounding_box_2d
{
  public:
    /**
     * Standard constructor that computes an initial bounding box.
     *
     * @param lyt Gate-level or cell-level layout whose bounding box is desired.
     */
    explicit bounding_box_2d(const Lyt& lyt) noexcept : layout{lyt}
    {
        static_assert(is_coordinate_layout_v<Lyt>, "Lyt is not a coordinate layout");
        update_bounding_box();
    }
    /**
     * The bounding box is not automatically updated when the layout changes. This function recomputes the bounding box.
     */
    void update_bounding_box()
    {
        min = {0, 0, 0};
        max = {0, 0, 0};

        // empty layouts don't need further computation
        if (layout.is_empty())
        {
            return;
        }

        // the layout is based on SiQAD coordinates
        if constexpr (has_siqad_coord_v<Lyt>)
        {
            int32_t min_x = std::numeric_limits<int32_t>::max();
            int32_t max_x = std::numeric_limits<int32_t>::min();

            int32_t min_y = std::numeric_limits<int32_t>::max();
            int32_t max_y = std::numeric_limits<int32_t>::min();

            uint8_t min_z = 1;
            uint8_t max_z = 0;

            layout.foreach_cell(
                [&min_x, &max_x, &min_y, &max_y, &min_z, &max_z](const auto& c)
                {
                    if (c.x < min_x)
                    {
                        min_x = c.x;
                    }
                    if (c.x > max_x)
                    {
                        max_x = c.x;
                    }

                    if (c.y == min_y && c.z < min_z)
                    {
                        min_z = c.z;
                    }
                    if (c.y < min_y)
                    {
                        min_y = c.y;
                        min_z = c.z;
                    }

                    if (c.y == max_y && c.z > max_z)
                    {
                        max_z = c.z;
                    }
                    if (c.y > max_y)
                    {
                        max_y = c.y;
                        max_z = c.z;
                    }
                });
            min = {min_x, min_y, min_z};
            max = {max_x, max_y, max_z};
        }
        else
        {
            // set min to max coordinate in the layout
            min = {layout.x(), layout.y()};

            layout.foreach_coordinate(
                [this](const auto& c)
                {
                    if (!is_empty_coordinate(c))
                    {
                        if (c.x < min.x)
                        {
                            min.x = c.x;
                        }
                        if (c.y < min.y)
                        {
                            min.y = c.y;
                        }
                        if (c.x > max.x)
                        {
                            max.x = c.x;
                        }
                        if (c.y > max.y)
                        {
                            max.y = c.y;
                        }
                    }
                });
        }

        x_size = max.x - min.x;
        y_size = max.y - min.y;
    }
    /**
     * Returns the minimum corner of the bounding box.
     *
     * In a `cartesian_layout<offset::ucoord_t>` object, this location represents the most north-western coordinate of
     * the bounding box enclosing every non-empty coordinate.
     *
     * @return The minimum enclosing coordinate in the associated layout.
     */
    [[nodiscard]] coordinate<Lyt> get_min() const noexcept
    {
        return min;
    }
    /**
     * Returns the maximum corner of the bounding box.
     *
     * In a `cartesian_layout<offset::ucoord_t>` object, this location represents the most south-eastern coordinate of
     * the bounding box enclosing every non-empty coordinate.
     *
     * @return The maximum enclosing coordinate in the associated layout.
     */
    [[nodiscard]] coordinate<Lyt> get_max() const noexcept
    {
        return max;
    }
    /**
     * Returns the horizontal size of the bounding box in layout coordinates.
     *
     * @return Bounding box size along the x-axis.
     */
    [[nodiscard]] auto get_x_size() const noexcept
    {
        return x_size;
    }
    /**
     * Returns the vertical size of the bounding box in layout coordinates.
     *
     * @return Bounding box size along the y-axis.
     */
    [[nodiscard]] auto get_y_size() const noexcept
    {
        return y_size;
    }

  private:
    const Lyt& layout;

    coordinate<Lyt> min{0, 0, 0}, max{0, 0, 0};

    decltype(min.x) x_size{};
    decltype(min.y) y_size{};

    [[nodiscard]] bool is_empty_coordinate(const coordinate<Lyt>& c) const noexcept
    {
        static_assert(has_is_empty_tile_v<Lyt> || has_is_empty_cell_v<Lyt>,
                      "Lyt does neither implement the is_empty_tile nor the is_empty_cell function");

        if constexpr (has_is_empty_tile_v<Lyt>)
        {
            return layout.is_empty_tile(c);
        }
        else if constexpr (has_is_empty_cell_v<Lyt>)
        {
            return layout.is_empty_cell(c);
        }

        return false;
    }
};

}  // namespace fiction

#pragma GCC diagnostic pop

#endif  // FICTION_BOUNDING_BOX_HPP

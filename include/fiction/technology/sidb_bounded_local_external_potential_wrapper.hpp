//
// Created by Willem Lambooy on 23/04/2025.
//

#ifndef SIDB_BOUNDED_LOCAL_EXTERNAL_POTENTIAL_WRAPPER_HPP
#define SIDB_BOUNDED_LOCAL_EXTERNAL_POTENTIAL_WRAPPER_HPP

#include <fiction/traits.hpp>

#include <type_traits>
#include <unordered_map>

namespace fiction
{

/**
 * This struct stores parameters for the `sidb_bounded_local_external_potential_wrapper`
 */
struct sidb_bounded_local_external_potential_wrapper_params
{
};

/**
 * A layout type to layer on top of any SiDB cell-level layout. It implements an interface to store and access
 * fabrication defects on the H-Si(100) 2x1 surface.
 *
 * @tparam Lyt SiDB cell-level layout type.
 * @tparam has_sidb_defect_surface Automatically determines whether a defect interface is already present.
 */
template <typename Lyt,
          bool has_sidb_defect_surface = std::conjunction_v<has_assign_sidb_defect<Lyt>, has_get_sidb_defect<Lyt>>>
class sidb_bounded_local_external_potential_wrapper : public Lyt
{};

template <typename Lyt>
class sidb_bounded_local_external_potential_wrapper<Lyt, true> : public Lyt
{
public:
    explicit sidb_bounded_local_external_potential_wrapper(const Lyt& lyt, [[maybe_unused]] const sidb_bounded_local_external_potential_wrapper_params& ps = {}) : Lyt(lyt)
    {}
};

template <typename Lyt>
class sidb_bounded_local_external_potential_wrapper<Lyt, false> : public Lyt
{
  public:
    struct sidb_bounded_local_external_potential_wrapper_storage
    {
        explicit sidb_bounded_local_external_potential_wrapper_storage(sidb_bounded_local_external_potential_wrapper_params ps = {}) : params(std::move(ps)) {}

        sidb_bounded_local_external_potential_wrapper_params params{};

        std::unordered_map<typename Lyt::coordinate, sidb_defect> defective_coordinates{};
    };

    using storage = std::shared_ptr<sidb_bounded_local_external_potential_wrapper_storage>;

    /**
     * Standard constructor for empty layouts.
     *
     * @param ps SiDB defect surface parameters.
     */
    explicit sidb_bounded_local_external_potential_wrapper(const sidb_bounded_local_external_potential_wrapper_params& ps = {}) :
            Lyt(),
            strg{std::make_shared<sidb_bounded_local_external_potential_wrapper_storage>(ps)}
    {
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
        static_assert(is_charge_distribution_surface_v<Lyt>, "Lyt is not a charge distribution surface");
    }

    /**
     * Standard constructor that layers the SiDB defect interface onto a layout with aspect ratio ar as input.
     *
     * @param ar aspect ratio of the layout.
     * @param ps SiDB defect surface parameters.
     */
    explicit sidb_bounded_local_external_potential_wrapper(const typename Lyt::aspect_ratio& ar, const sidb_bounded_local_external_potential_wrapper_params& ps = {}) :
            Lyt(ar),
            strg{std::make_shared<sidb_bounded_local_external_potential_wrapper_storage>(ps)}
    {
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
        static_assert(is_charge_distribution_surface_v<Lyt>, "Lyt is not a charge distribution surface");
    }

    /**
     * Standard constructor that layers the SiDB defect interface onto an existing layout.
     *
     * @param lyt Existing layout that is to be extended by an SiDB defect interface.
     * @param ps SiDB defect surface parameters.
     */
    explicit sidb_bounded_local_external_potential_wrapper(const Lyt& lyt, const sidb_bounded_local_external_potential_wrapper_params& ps = {}) :
            Lyt(lyt),
            strg{std::make_shared<sidb_bounded_local_external_potential_wrapper_storage>(ps)}
    {
        static_assert(is_cell_level_layout_v<Lyt>, "Lyt is not a cell-level layout");
        static_assert(has_sidb_technology_v<Lyt>, "Lyt is not an SiDB layout");
        static_assert(is_charge_distribution_surface_v<Lyt>, "Lyt is not a charge distribution surface");
    }
    /**
     * Clones the layout returning a deep copy.
     *
     * @return Deep copy of the layout.
     */
    [[nodiscard]] sidb_bounded_local_external_potential_wrapper clone() const noexcept
    {
        sidb_bounded_local_external_potential_wrapper copy{Lyt::clone()};
        copy.strg = std::make_shared<sidb_bounded_local_external_potential_wrapper_storage>(*strg);

        return copy;
    }

  private:
    storage strg;
};

template <class T>
sidb_bounded_local_external_potential_wrapper(const T&) -> sidb_bounded_local_external_potential_wrapper<T>;

}  // namespace fiction

#endif //SIDB_BOUNDED_LOCAL_EXTERNAL_POTENTIAL_WRAPPER_HPP

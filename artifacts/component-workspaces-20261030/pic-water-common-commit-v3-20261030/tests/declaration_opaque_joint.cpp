#include "closed_parent_joint_step.hpp"
#include <type_traits>
using Trial=raikou::native_water::ClosedParentTrial;
static_assert(!std::is_default_constructible_v<Trial>);
static_assert(!std::is_copy_constructible_v<Trial>);
static_assert(!std::is_copy_assignable_v<Trial>);
static_assert(std::is_nothrow_move_constructible_v<Trial>);
static_assert(std::is_nothrow_move_assignable_v<Trial>);
template<class T>concept ExposesNative=requires(T&t){t.native_;};
template<class T>concept ExposesField=requires(T&t){t.field_heat_;};
template<class T>concept ExposesRng=requires(T&t){t.random_;};
static_assert(!ExposesNative<Trial> && !ExposesField<Trial> && !ExposesRng<Trial>);
static_assert(std::is_same_v<decltype(std::declval<const Trial&>().audit()),
  const raikou::native_water::Json&>);
// Declaration-only qualification. No accepted owner or numerical operator is
// constructed, prepared, published or advanced by this translation unit.

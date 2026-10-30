#include "closed_parent_joint_step.hpp"
#include <type_traits>
using Trial=raikou::native_water::ClosedParentTrial;
template<class T> concept PublicField = requires(T& t){t.field_heat;};
template<class T> concept PublicNative = requires(T& t){t.native;};
template<class T> concept PublicRandom = requires(T& t){t.random;};
template<class T> concept PublicCensus = requires(T& t){t.post_first_census;};
template<class T> concept PublicStep = requires(T& t){t.private_step;};
template<class T> concept PublicExpected = requires(T& t){t.expected;};
static_assert(!PublicField<Trial> && !PublicNative<Trial> && !PublicRandom<Trial>);
static_assert(!PublicCensus<Trial> && !PublicStep<Trial> && !PublicExpected<Trial>);
static_assert(!std::is_default_constructible_v<Trial>);
static_assert(!std::is_copy_constructible_v<Trial> && !std::is_copy_assignable_v<Trial>);
static_assert(std::is_nothrow_move_constructible_v<Trial> && std::is_nothrow_move_assignable_v<Trial>);
static_assert(std::is_same_v<decltype(std::declval<const Trial&>().audit()),
                             const raikou::native_water::Json&>);

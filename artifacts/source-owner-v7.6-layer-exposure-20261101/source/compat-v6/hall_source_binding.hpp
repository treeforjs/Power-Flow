#pragma once
// Compatibility include used only to compile the sealed source v6 qualification test against
// owner v7 (benchmark settings): the sealed binding name resolves to the v7 binding.
#include "pf_arrival_v11/arrival_owner.hpp"
#include "hall_source_binding_v7.hpp"
namespace hall_source_binding = hall_source_binding_v7;

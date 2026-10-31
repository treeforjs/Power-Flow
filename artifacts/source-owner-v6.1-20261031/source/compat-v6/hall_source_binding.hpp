#pragma once
// Compatibility include used only to compile the sealed source v6 qualification test against
// owner v6.1: the sealed binding name resolves to the v6.1 binding (same bodies, v6.1 types).
#include "pf_arrival_v11/arrival_owner.hpp"
#include "hall_source_binding_v6_1.hpp"
namespace hall_source_binding = hall_source_binding_v6_1;

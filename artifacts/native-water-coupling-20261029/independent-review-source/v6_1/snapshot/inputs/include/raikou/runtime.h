#pragma once
#include "raikou/json.h"
#include <string_view>
namespace raikou {
Json execute(std::string_view operation, const Json& request = Json::object());
}

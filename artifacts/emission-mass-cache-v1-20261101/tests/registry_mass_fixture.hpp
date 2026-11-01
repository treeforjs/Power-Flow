#pragma once
#include <bit>
#include <cstdint>
#include "surface_emission.hpp"
inline surface_emission::Config registry_mass_fixture(surface_emission::Config c){c.electron_mass_kg=std::bit_cast<double>(std::uint64_t{4157519395190634280});c.electron_mass_identity="source-v7.1/blob-electron-mass-v1/e5ddf4f6a2da40cfff27cf41f07160a2d53630a195c50754f081444b44a8776b";return c;}

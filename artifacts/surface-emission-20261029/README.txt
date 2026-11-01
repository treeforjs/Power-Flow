Read MODEL.txt for equations, primary sources, domain/approximation, phase space,
energy/current ownership and strong-B subgrid integration. Public API is
include/surface_emission.hpp. The preferred conditional source interface is
prepare -> begin_layer -> advance_layer, followed by parent atomic commit of
electron inventory, current paths, complete wall heat and field/reference work.

Build and reproduce on one CPU7, nice15, one job, frozen userspace, holding
/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core7.lock.
The reproduce.sh helper does not acquire the lock itself. No dependencies beyond
the frozen C++20 compiler/standard library are required. No GPU/PIC state is
accessed by this library or its numerical tests.

Reports and interface JSON live outside projects at
/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-thermal-field-emission-20261029.
The parent owns integration/launch. All source and reports are candidates for
joint admission; physical_foil_model_admitted remains false.

// Read-only reconstruction of the sealed heat owner's Cartesian scalar Q1 map.
#include "mfem.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <sched.h>
#include <sys/resource.h>
#include <unistd.h>

static void need(bool value, const char* message) {
    if (!value) { throw std::invalid_argument(message); }
}

int main(int argc, char** argv) {
    try {
        need(argc == 8, "usage: pilot-heat-coordinates-v4 output.json nx ny x0 y0 width thickness");
        cpu_set_t mask; CPU_ZERO(&mask);
        need(sched_getaffinity(0, sizeof(mask), &mask) == 0 && CPU_COUNT(&mask) == 1
             && CPU_ISSET(7, &mask) && getpriority(PRIO_PROCESS, 0) >= 15,
             "requires CPU7/nice15 bounded execution under the external test lease");
        const std::filesystem::path output(argv[1]);
        need(!std::filesystem::exists(output), "coordinate certificate destination exists");
        const int nx = std::stoi(argv[2]), ny = std::stoi(argv[3]);
        const double x0 = std::stod(argv[4]), y0 = std::stod(argv[5]),
                     width = std::stod(argv[6]), thickness = std::stod(argv[7]);
        need(nx > 0 && ny > 0 && nx <= 4096 && ny <= 4096 && nx * ny <= 1000000,
             "invalid bounded Cartesian cell counts");
        need(std::isfinite(x0) && std::isfinite(y0) && std::isfinite(width)
             && std::isfinite(thickness) && width > 0 && thickness > 0,
             "invalid finite positive geometry");
        // Match the owner's constructor and mutation order exactly.
        mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL,
                                                     true, width, thickness);
        mfem::H1_FECollection collection(1, 2);
        mfem::FiniteElementSpace space(&mesh, &collection);
        for (int v = 0; v < mesh.GetNV(); ++v) {
            auto* x = mesh.GetVertex(v); x[0] += x0; x[1] += y0;
        }
        std::set<int> dofs;
        std::vector<double> xs, ys;
        for (int v = 0; v < mesh.GetNV(); ++v) {
            auto* x = mesh.GetVertex(v); xs.push_back(x[0]); ys.push_back(x[1]);
            mfem::Array<int> d; space.GetVertexDofs(v, d);
            need(d.Size() == 1 && d[0] >= 0 && d[0] < space.GetVSize()
                 && dofs.insert(d[0]).second, "scalar Q1 vertex/DOF map is not bijective");
        }
        need(int(dofs.size()) == space.GetVSize() && space.GetVSize() == (nx + 1) * (ny + 1),
             "native scalar Q1 checkpoint dimensions differ");
        for (auto* axis : {&xs, &ys}) {
            std::sort(axis->begin(), axis->end());
            axis->erase(std::unique(axis->begin(), axis->end()), axis->end());
        }
        need(xs.size() == std::size_t(nx + 1) && ys.size() == std::size_t(ny + 1),
             "native Cartesian axes are unresolved");
        std::ofstream file(output);
        file << std::setprecision(17) << "{\n\"schema\":\"native-pilot-Q1-vertex-dof-v4\",\n"
             << "\"mfem_version\":" << std::quoted(mfem::GetVersionStr()) << ",\n"
             << "\"executable_path\":" << std::quoted(std::filesystem::read_symlink("/proc/self/exe").string()) << ",\n"
             << "\"nx\":" << nx << ",\"ny\":" << ny << ",\"x0_m\":" << x0 << ",\"y0_m\":" << y0
             << ",\"width_m\":" << width << ",\"thickness_m\":" << thickness << ",\n"
             << "\"node_count\":" << space.GetVSize() << ",\"cell_count\":" << mesh.GetNE() << ",\n";
        auto array = [&](const char* name, const std::vector<double>& values) {
            file << std::quoted(name) << ":[";
            for (std::size_t i = 0; i < values.size(); ++i) { file << (i ? "," : "") << values[i]; }
            file << "],\n";
        };
        array("x_edges_m", xs); array("y_edges_m", ys);
        file << "\"vertices\":[\n";
        for (int v = 0; v < mesh.GetNV(); ++v) {
            auto* x = mesh.GetVertex(v); mfem::Array<int> d; space.GetVertexDofs(v, d);
            const int i = std::lower_bound(xs.begin(), xs.end(), x[0]) - xs.begin();
            const int j = std::lower_bound(ys.begin(), ys.end(), x[1]) - ys.begin();
            file << (v ? ",\n" : "") << "{\"vertex\":" << v << ",\"dof\":" << d[0]
                 << ",\"i\":" << i << ",\"j\":" << j << ",\"x_m\":" << x[0] << ",\"y_m\":" << x[1] << "}";
        }
        file << "\n],\"cpu\":7,\"nice\":" << getpriority(PRIO_PROCESS, 0)
             << ",\"owner_state_read_or_mutated\":false,\"numerical_solution_performed\":false}\n";
        file.close(); need(bool(file), "coordinate certificate write failed");
        std::cout << output << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; return 2;
    }
}

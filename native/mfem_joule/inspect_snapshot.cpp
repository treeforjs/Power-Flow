// Inspect native MFEM VisIt snapshots without assigning physical units.
#include "mfem.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

int main(int argc, char **argv)
{
   mfem::Mpi::Init(argc, argv);
   mfem::Hypre::Init();
   const char *name = "Joule";
   const char *prefix = "./";
   int cycle = 0;
   mfem::OptionsParser options(argc, argv);
   options.AddOption(&name, "-n", "--name", "VisIt collection basename");
   options.AddOption(&prefix, "-d", "--directory", "Collection parent directory");
   options.AddOption(&cycle, "-c", "--cycle", "Saved cycle");
   options.Parse();
   if (!options.Good()) { if (mfem::Mpi::Root()) options.PrintUsage(std::cerr); return 2; }
   mfem::VisItDataCollection data(MPI_COMM_WORLD, name);
   data.SetPrefixPath(prefix);
   data.Load(cycle);
   if (data.Error() || !data.GetMesh()) { return 3; }
   mfem::Mesh &mesh = *data.GetMesh();
   const char *fields[] = {"T", "w", "E", "B", "Phi", "F"};
   bool finite = true;
   double integrals[2] = {0., 0.};
   double extrema[4] = {std::numeric_limits<double>::infinity(),
                       -std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity(),
                       -std::numeric_limits<double>::infinity()};
   for (int f = 0; f < 6; ++f)
   {
      mfem::GridFunction *field = data.GetField(fields[f]);
      if (!field) { return 4; }
      for (int i = 0; i < field->Size(); ++i) { finite = finite && std::isfinite((*field)[i]); }
      if (f >= 2) { continue; }
      for (int e = 0; e < mesh.GetNE(); ++e)
      {
         auto *transform = mesh.GetElementTransformation(e);
         const auto &rule = mfem::IntRules.Get(mesh.GetElementBaseGeometry(e),
                                              2 * field->FESpace()->GetOrder(e) + 3);
         for (int q = 0; q < rule.GetNPoints(); ++q)
         {
            const auto &ip = rule.IntPoint(q);
            transform->SetIntPoint(&ip);
            const double value = field->GetValue(e, ip);
            finite = finite && std::isfinite(value);
            integrals[f] += value * ip.weight * transform->Weight();
            extrema[2*f] = std::min(extrema[2*f], value);
            extrema[2*f+1] = std::max(extrema[2*f+1], value);
         }
      }
   }
   int local_finite = finite ? 1 : 0, all_finite = 0;
   MPI_Allreduce(&local_finite, &all_finite, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
   double global_integrals[2];
   MPI_Allreduce(integrals, global_integrals, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   double global_extrema[4];
   for (int i = 0; i < 4; ++i)
   {
      MPI_Allreduce(extrema + i, global_extrema + i, 1, MPI_DOUBLE,
                    i % 2 == 0 ? MPI_MIN : MPI_MAX, MPI_COMM_WORLD);
   }
   if (mfem::Mpi::Root())
   {
      std::cout << std::setprecision(17)
         << "{\"schema\":\"mhdlab-mfem-joule-native-inspection-v1\","
         << "\"physical_units_assigned\":false,\"all_fields_finite\":" << (all_finite ? "true" : "false")
         << ",\"cycle\":" << data.GetCycle() << ",\"time_native\":" << data.GetTime()
         << ",\"temperature_volume_integral_native\":" << global_integrals[0]
         << ",\"joule_volume_integral_native\":" << global_integrals[1]
         << ",\"temperature_min_native\":" << global_extrema[0]
         << ",\"temperature_max_native\":" << global_extrema[1]
         << ",\"joule_min_native\":" << global_extrema[2]
         << ",\"joule_max_native\":" << global_extrema[3] << "}\n";
   }
   return all_finite && std::isfinite(global_integrals[0]) && std::isfinite(global_integrals[1]) ? 0 : 5;
}

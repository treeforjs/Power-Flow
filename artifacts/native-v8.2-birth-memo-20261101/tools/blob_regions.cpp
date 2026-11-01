// Lists the coverage regions of a chemistry v3 blob with their bound status (development aid).
#include "chem_host.hpp"
#include <cstdio>
int main(int argc, char **argv) {
  const auto b = pfc::load_blob(argv[1]);
  std::vector<int> any(b.n_regions, 0), bounded(b.n_regions, 0), unbounded(b.n_regions, 0);
  for (pfc_u32 g = 0; g < b.n_groups; ++g) {
    const auto &gr = b.view.groups[g];
    for (pfc_u32 i = 0; i <= gr.grid_n; ++i) {
      const auto reg = b.cover.bin_region[b.cover.bin_off[g] + i];
      if (reg == 0xFFFFFFFFu) continue;
      any[reg]++;
      if (b.cover.bin_row[b.cover.bin_off[g] + i] != 0xFFFFFFFFu) bounded[reg]++; else unbounded[reg]++;
    }
  }
  for (pfc_u32 r = 0; r < b.n_regions; ++r)
    std::printf("%3u bins=%d bounded=%d unbounded=%d %s\n", r, any[r], bounded[r], unbounded[r], b.region_names[r].c_str());
  for (pfc_u32 q = 0; q < b.n_proc; ++q) {
    const auto &pr = b.view.procs[q];
    const auto proj = b.view.groups[pr.group].projectile;
    std::printf("proc %u kind %u [%s]:", q, pr.kind, b.process_names[q].c_str());
    for (pfc_u32 i = 0; i < pr.n_out; ++i)
      std::printf(" %s/%u%s", b.species_names[pr.out_species[i]].c_str(), pr.out_dest[i],
                  (pr.out_dest[i] == PFC_DEST_SURVIVOR && pr.out_species[i] != proj) ? "(SPECIES CHANGE)" : "");
    std::printf("\n");
  }
  for (pfc_u32 g = 0; g < b.n_groups; ++g)
    std::printf("group %u %s + %s procs %u\n", g, b.species_names[b.view.groups[g].projectile].c_str(),
                b.species_names[b.view.groups[g].target].c_str(), b.view.groups[g].proc_n);
}
// (appended) survivor species changes and destinations per process
struct ProcDump {
  ProcDump() {}
};

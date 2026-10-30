Matched solid magnetic/heating quick runs and frozen electron-orbit checks.

128x16 and256x32 cells per foil; both run continuously from0.749788055492ns through60ns with5mmwidth,0.2mmthickness,5mmgap,25mm current length and±40mm exterior. Both temperature histories and native Az/rate fields originate in each same run. Two real MPI ranks are bound to OS CPUs6/7 via translated hwloc PU indices12/14. Numerical dependencies remain frozen.

Both direct electrode surface histories are retained; the solid cell-average bulk export covers the upper foil only. Hot corner/patch values remain mesh sensitive, and material properties are constant with no phase change. These are quick conditional studies, not material/grid-converged predictions.

Raw native runs are preserved across48MiB gzip-tar parts with checksums; analysis-source-and-previews.tar.gz contains replay/analysis/plotting code, hashed evidence, preview images and frozen pic_boris unit benchmarks. No third-party binary libraries are included. The orbit checks do not run PIC or model the full pulsed finite target. The primary long PIC run remains paused and its original inputs/checkpoints are preserved separately. New heating inputs require a fresh consistent PIC history if adopted.

Path and wording normalization for merged branches — 2026-10-31

The side branches merged into master on 2026-10-31 (electrode-heat-source-v4
with its predecessors, surface-emission, ss304-thermal-expansion and the five
ss304 archive branches) predate the 2026-10-30 normalization. The same rules
were applied to their files:

  run workspace       -> /home/treefor/sim-workspace/powerflow
  other local sources -> /home/treefor/sim-workspace/source-backups/<name>
  Windows worktrees   -> <windows-worktree>

One third-person phrase in ss304_ethos_eos.md was rewritten in a neutral
voice. No numerical values, inputs or results changed. Hash values recorded
inside earlier receipts and manifests refer to the pre-normalization file
contents, which remain available in each branch's last commit (listed in
docs/coupled_pic_handoff.md). files.json lists the sha256 of each changed
file before and after normalization.

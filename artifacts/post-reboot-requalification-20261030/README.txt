Post-reboot re-qualification — 2026-10-30

The workstation rebooted at 12:09 MDT after the corrective component work. Host
identity (kernel 7.2.11-arch1-1, NVIDIA 615.78.08 module srcversion
4200B1C66D4D12038311105) matches the pre-integration audit observation, and the
frozen userspace root matches its sealed manifest (246,673 entries, none
missing, extra or changed).

The reboot ended the stopped electron-only GPU process. Its accepted checkpoint
checkpoint_437462.bin at 26.4787 ns is intact (sha256 matches the recorded
value in both copies); 574 uncheckpointed steps after it were lost. Resume is
not authorized. paused-production-state.json supersedes the audit's
preserved-production.json.

Every re-runnable corrective suite passes again on the current host with
identical scientific outputs (results.json). The SS304 melt runner was not run
because it asserts the pre-reboot process state. reruns/ holds the logs and
re-run evidence; reruns/run-cpu7-serial.sh is the serial CPU7 driver, shown
with workspace paths normalized (inside the frozen sandbox the physical
workspace path was used). This confirms existing evidence only; no further
audit finding is closed by it.

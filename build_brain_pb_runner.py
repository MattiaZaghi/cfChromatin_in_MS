#!/usr/bin/env python3
"""
Run the shipped marker builder with the brain_pb atlas swapped in for bulk brain.

Drop-in for build_tissue_markers_sadeh.py: same CLI, same output, same flags.

    python3 build_brain_pb_runner.py --top-n 400 --out <dir> ...

Why a runner instead of editing the builder, or copying it
---------------------------------------------------------
build_tissue_markers_sadeh.py ships with brain enabled and the brain_pb entry
commented out (FRAGMENT_TISSUES, and ROADMAP_EIDS["brain"]). Editing it in place
is not an option: marker_n_sweep_brain.sh re-runs that same file at every level,
so changing the tissue set there would change which panels exist mid-sweep.
Copying the 1500-line builder is worse -- two copies drift, and the bulk and CNS
copies have already drifted ~245 substantive lines apart exactly this way.

So: import the real builder under a private name, mutate its tissue maps, call
its main(). One place, no copy, and it tracks upstream changes automatically.

The mutation must be IN PLACE. build_tissue_replicates() reads ROADMAP_EIDS /
BAM_TISSUES / FRAGMENT_TISSUES out of its OWN module globals when it is called,
so rebinding those names here would have no effect on it. .pop()/assignment on
the real module's dicts does.

This is the same swap that
tissue_markers_sadeh_brain_pb/decide_n_markers/bld_override/ performs for the
empirical-N scripts; that one is an import shim, this one is an entry point.
Keep the two consistent if the tissue maps ever change.
"""

import importlib.util
import sys
from pathlib import Path

_REAL = Path(__file__).resolve().parent / "build_tissue_markers_sadeh.py"
if not _REAL.exists():                                       # pragma: no cover
    raise SystemExit(f"[pb-runner] shipped builder not found: {_REAL}")

# Load under a private name. Reusing "build_tissue_markers_sadeh" would be fine
# here (this file has a different name) but the private name makes it obvious in
# a traceback which module is the real builder and which is the wrapper.
_spec = importlib.util.spec_from_file_location("_bld_shipped", _REAL)
_real = importlib.util.module_from_spec(_spec)
sys.modules["_bld_shipped"] = _real
_spec.loader.exec_module(_real)

# --- the swap ---------------------------------------------------------------
# Drop bulk brain's Roadmap epigenomes, add the brain_pb fragment atlas. After
# this the reference emits a brain_pb panel and no brain panel, which is what
# makes gt_brain_pb_* samples scoreable on target.
_dropped = _real.ROADMAP_EIDS.pop("brain", None)
if _dropped is None:                                         # pragma: no cover
    print("[pb-runner] WARNING: 'brain' was not in ROADMAP_EIDS -- the shipped "
          "builder changed; verify this runner still does what it claims.",
          file=sys.stderr)
_real.FRAGMENT_TISSUES["brain_pb"] = _real.BRAIN_PB_FILES

print(f"[pb-runner] brain -> brain_pb  (dropped Roadmap EIDs: {_dropped})",
      file=sys.stderr)
print(f"[pb-runner] {len(_real.ROADMAP_EIDS)} roadmap + "
      f"{len(_real.FRAGMENT_TISSUES)} fragment tissues", file=sys.stderr)

# main() parses sys.argv[1:] itself, and argv is already this script's argv, so
# every flag the caller passed reaches the builder untouched.
if __name__ == "__main__":
    _real.main()

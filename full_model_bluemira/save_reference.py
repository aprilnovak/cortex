#!/usr/bin/env python3
"""
save_reference.py
=================
Snapshot current post-processing results as a named gold run.

By default reads SIM_TYPE and BREEDER_TYPE directly from inputs.py, so it always
snapshots whichever simulation is currently configured — no manual editing needed.

Expected source layout (written by neutronics_post.py + depletion_post.py):

    <SIM_TYPE>/<BREEDER_TYPE>/
        neutronics_results/
            <chunk_key>/
                profiles/          
                plots/
                spectra/
                ...
        depletion_results/
            <chunk_key>/
                activity/          
                decay_heat/       

The reference copy mirrors this layout exactly under:

    reference_runs/<label>/<SIM_TYPE>/<BREEDER_TYPE>/

"""

import json
import shutil
from datetime import datetime
from pathlib import Path

# ───────────────────────────────────────────────────
# Inputs
# ───────────────────────────────────────────────────
LABEL        = "W_armor_eurofer_structural"
ARMOR        = "W"
STRUCTURAL   = "eurofer97"

# Leave as None to use whatever inputs.py currently says.
SIM_TYPE_OVERRIDE     = None   # e.g. "tokamak" | "slab" | "fast_slab"
BREEDER_TYPE_OVERRIDE = None   # e.g. "WCLL"    | "HCLL" | "HCPB"

# ───────────────────────────────────────────────────
# Directory allocation
# ───────────────────────────────────────────────────
BASE_DIR = Path(__file__).parent

import sys
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))
import inputs as cfg

SIM_TYPE     = SIM_TYPE_OVERRIDE     or cfg.SIM_TYPE
BREEDER_TYPE = BREEDER_TYPE_OVERRIDE or cfg.BREEDER_TYPE

print(f"[save_reference] SIM_TYPE={SIM_TYPE}  BREEDER_TYPE={BREEDER_TYPE}")

# Source and destination paths
SIM_DIR        = BASE_DIR / SIM_TYPE / BREEDER_TYPE
src_neutronics = SIM_DIR / "neutronics_results"
src_depletion  = SIM_DIR / "depletion_results"

dst_root = BASE_DIR / "reference_runs" / LABEL
dst_sim  = dst_root / SIM_TYPE / BREEDER_TYPE

# Guard 
if dst_sim.exists():
    ans = input(
        f"[warn] '{dst_sim.relative_to(BASE_DIR)}' already exists. "
        f"Overwrite? [y/N] "
    ).strip().lower()
    if ans != "y":
        print("Aborted.")
        raise SystemExit(0)
    shutil.rmtree(dst_sim)

dst_sim.mkdir(parents=True)

# ───────────────────────────────────────────────────
# Copy result trees 
# ───────────────────────────────────────────────────
for src, name in [
    (src_neutronics, "neutronics_results"),
    (src_depletion,  "depletion_results"),
]:
    if src.exists():
        shutil.copytree(src, dst_sim / name)
        print(
            f"[ok] copied  {src.relative_to(BASE_DIR)}"
            f"  →  {(dst_sim / name).relative_to(BASE_DIR)}"
        )
    else:
        print(f"[skip] {src.relative_to(BASE_DIR)} not found — skipping")

# ───────────────────────────────────────────────────
# Validate data
# ───────────────────────────────────────────────────
# Neutronics: profile CSVs under profiles/ 
dst_neutronics = dst_sim / "neutronics_results"
if dst_neutronics.exists():
    n_profiles = len(list(dst_neutronics.rglob("profiles/profile_*.csv")))
    if n_profiles:
        print(f"[ok] neutronics: {n_profiles} profile CSV(s) found under profiles/")
    else:
        print("[warn] neutronics: no profiles/profile_*.csv found "
              "— did neutronics_post.py complete?")

# Depletion: CSVs under activity/ and decay_heat/
dst_depletion = dst_sim / "depletion_results"
if dst_depletion.exists():
    n_act = len(list(dst_depletion.rglob("activity/activity_all_cells.csv")))
    n_dh  = len(list(dst_depletion.rglob("decay_heat/decayheat_all_cells.csv")))
    if n_act:
        print(f"[ok] depletion: {n_act} activity_all_cells.csv found under activity/")
    else:
        print("[warn] depletion: no activity/activity_all_cells.csv found "
              "— did depletion_post.py complete?")
    if n_dh:
        print(f"[ok] depletion: {n_dh} decayheat_all_cells.csv found under decay_heat/")
    else:
        print("[warn] depletion: no decay_heat/decayheat_all_cells.csv found "
              "— did depletion_post.py complete?")

# ───────────────────────────────────────────────────
# Write per-slot metadata 
# ───────────────────────────────────────────────────
meta = {
    "label":        LABEL,
    "saved_at":     datetime.now().isoformat(),
    "sim_type":     SIM_TYPE,
    "breeder_type": BREEDER_TYPE,
    "armor":        ARMOR,
    "structural":   STRUCTURAL,
}
meta_path = dst_root / f"meta_{SIM_TYPE}_{BREEDER_TYPE}.json"
meta_path.write_text(json.dumps(meta, indent=2))
print(f"[ok] metadata → {meta_path.relative_to(BASE_DIR)}")

print(f"\nReference slot ready → {dst_sim.relative_to(BASE_DIR)}")
print("Re-run (after changing inputs.py) to add more slots to this label.")
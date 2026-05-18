#!/usr/bin/env python3
"""
reference_manager.py
====================
Load, list, and inspect named gold runs produced by save_reference.py.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import inputs as cfg

REFERENCE_ROOT: Path = cfg.BASE_DIR / "reference_runs"

# ──────────────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────────────
def load_reference(
    label: str,
    *,
    sim_type: Optional[str]     = None,
    breeder_type: Optional[str] = None,
) -> Optional[dict]:
    """
    Load a named gold run for the given sim_type / breeder_type slot.

    Falls back to cfg.SIM_TYPE and cfg.BREEDER_TYPE when not supplied.

    Returns
    -------
    dict with keys:
        "meta"           : metadata dict
        "neutronics_dir" : Path to neutronics_results/  (or None if absent)
        "depletion_dir"  : Path to depletion_results/   (or None if absent)
    or None if the slot directory does not exist.
    """
    sim_type     = sim_type     or cfg.SIM_TYPE
    breeder_type = breeder_type or cfg.BREEDER_TYPE

    slot = REFERENCE_ROOT / label / sim_type / breeder_type

    # Will always flag before reference run is stored
    if not slot.exists():
        print(
            f"[reference] WARNING: slot '{label}/{sim_type}/{breeder_type}' not found "
            f"under {REFERENCE_ROOT} — overlay disabled."
        )
        return None

    meta_path = REFERENCE_ROOT / label / f"meta_{sim_type}_{breeder_type}.json"
    meta: dict = {}
    if meta_path.is_file():
        meta = json.loads(meta_path.read_text())
    else:
        print(f"[reference] WARNING: metadata file not found: {meta_path.name}")

    neutronics_dir = slot / "neutronics_results"
    depletion_dir  = slot / "depletion_results"

    return {
        "meta":           meta,
        "neutronics_dir": neutronics_dir if neutronics_dir.exists() else None,
        "depletion_dir":  depletion_dir  if depletion_dir.exists()  else None,
    }

def _count_slot_files(slot: Path) -> dict:
    """
    Count the key output files in a reference slot so list_references()
    can report slot completeness.

    Returns a dict with integer counts for each category.
    """
    counts = {
        "n_neutronics_profiles": 0,
        "n_activity_csvs":       0,
        "n_decayheat_csvs":      0,
    }

    neutronics_dir = slot / "neutronics_results"
    if neutronics_dir.exists():
        counts["n_neutronics_profiles"] = len(
            list(neutronics_dir.rglob("profiles/profile_*.csv"))
        )

    depletion_dir = slot / "depletion_results"
    if depletion_dir.exists():
        counts["n_activity_csvs"] = len(
            list(depletion_dir.rglob("activity/activity_all_cells.csv"))
        )
        counts["n_decayheat_csvs"] = len(
            list(depletion_dir.rglob("decay_heat/decayheat_all_cells.csv"))
        )

    return counts

def list_references() -> list[dict]:
    """Return a list of dicts describing every saved slot across all labels."""
    if not REFERENCE_ROOT.exists():
        return []

    rows = []
    for label_dir in sorted(REFERENCE_ROOT.iterdir()):
        if not label_dir.is_dir():
            continue
        label = label_dir.name
        for sim_dir in sorted(label_dir.iterdir()):
            if not sim_dir.is_dir() or sim_dir.name.startswith("meta"):
                continue
            for bt_dir in sorted(sim_dir.iterdir()):
                if not bt_dir.is_dir():
                    continue
                meta_path = label_dir / f"meta_{sim_dir.name}_{bt_dir.name}.json"
                meta = {}
                if meta_path.is_file():
                    meta = json.loads(meta_path.read_text())

                counts = _count_slot_files(bt_dir)

                rows.append({
                    "label":               label,
                    "sim_type":            sim_dir.name,
                    "breeder_type":        bt_dir.name,
                    "saved_at":            meta.get("saved_at", "—"),
                    "armor":               meta.get("armor",    "—"),
                    "structural":          meta.get("structural", "—"),
                    "neutronics_profiles": counts["n_neutronics_profiles"],
                    "activity_csvs":       counts["n_activity_csvs"],
                    "decayheat_csvs":      counts["n_decayheat_csvs"],
                })
    return rows

def print_references() -> None:
    """Pretty-print all saved reference slots to stdout."""
    rows = list_references()
    if not rows:
        print("[reference] No reference runs found under", REFERENCE_ROOT)
        return

    keys = (
        "label", "sim_type", "breeder_type", "saved_at",
        "armor", "structural",
        "neutronics_profiles", "activity_csvs", "decayheat_csvs",
    )
    hdrs = (
        "label", "sim_type", "breeder", "saved_at",
        "armor", "structural",
        "n_profiles", "n_act", "n_dh",
    )
    col_w = [
        max(len(str(r[k])) for r in rows + [dict(zip(keys, hdrs))])
        for k in keys
    ]
    fmt = "  ".join(f"{{:<{w}}}" for w in col_w)
    print(fmt.format(*hdrs))
    print("─" * (sum(col_w) + 2 * (len(col_w) - 1)))
    for r in rows:
        print(fmt.format(*[str(r[k]) for k in keys]))

# ──────────────────────────────────────────────────────────────────────────────
# python reference_manager.py
# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print_references()
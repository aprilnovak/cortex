#!/usr/bin/env python3
"""
Post-process OpenMC statepoint and save plots/CSVs per requested chunk key.

Outputs saved to:
  neutronics_results/<chunk_key>/
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import openmc

# ----------------------------
# Load neutronics model
# ----------------------------
import neutronics_model as bm

# ----------------------------
# User settings / paths
# ----------------------------
BASE_DIR = Path.cwd()

RUN_DIR = (BASE_DIR / "neutronics_run")

# find last touched statepoint.h5
def find_latest_statepoint(run_dir: Path) -> Path:
    candidates = list(run_dir.glob("statepoint.*.h5"))
    if not candidates:
        raise FileNotFoundError(f"No statepoint files found in: {run_dir}")
    candidates = sorted(candidates, key=lambda p: p.stat().st_mtime)
    latest = candidates[-1]
    print(f"  Last modified: {__import__('datetime').datetime.fromtimestamp(latest.stat().st_mtime)}")
    return latest

STATEPOINT_FILE = find_latest_statepoint(RUN_DIR)
print(f"Using statepoint: {STATEPOINT_FILE.name}")

INPUT_JSON = (BASE_DIR / "EUDEMO_HCPB_inputs.json")

RESULTS_DIR = BASE_DIR / "neutronics_results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

if not STATEPOINT_FILE.is_file():
    raise FileNotFoundError(
        f"Statepoint file not found: {STATEPOINT_FILE}\n"
        f"Expected neutronics output inside: {RUN_DIR}"
    )

# =============================================================================
# Choose which chunks to process
# =============================================================================
KEYS_TO_PROCESS = [
    "IB_1_b4",
    "OB_1_b6",
]

# layer poloidal analysis (region 1 only)
vv_index = -1
layers = [("Armor", 0), ("First_Wall", 1), ("VV", vv_index)]

# =============================================================================
# Optional: exclude surfaces from albedo summary
# =============================================================================
EXCLUDED_SURFACES = set()

# =============================================================================
# Source intensity -> scaling
# =============================================================================
total_power = 2e9
number_sectors = 16
section_power = total_power / number_sectors
ev_to_joule = 1.60218e-19
ev_fusion = 17.6e6
convert_e = ev_to_joule * ev_fusion
neutron_source_rate = section_power / convert_e
s_in_y = 365 * 24 * 60 * 60

energies = openmc.mgxs.GROUP_STRUCTURES["CCFE-709"]
unit_lethargy = np.array(
    [np.log(energies[i + 1] / energies[i]) for i in range(len(energies) - 1)],
    dtype=float,
)
E_mid = 0.5 * (
    np.asarray(energies[:-1], dtype=float) + np.asarray(energies[1:], dtype=float)
)

# =============================================================================
# Helpers
# =============================================================================
def subset_by_cells(data: np.ndarray, cell_bins: list[int], desired_cells: list[int]) -> np.ndarray:
    idx = [cell_bins.index(int(cid)) for cid in desired_cells]
    return data[idx]

def get_cell_bins_from_tally(tally: openmc.Tally) -> list[int]:
    cf = next(f for f in tally.filters if isinstance(f, openmc.CellFilter))
    return [int(x) for x in cf.bins]

def scaling_for_cells(cell_ids: list[int]) -> dict[int, float]:
    scaling: dict[int, float] = {}
    all_cells = getattr(bm, "all_cells", None)
    if all_cells is None:
        all_cells = bm.model.geometry.get_all_cells()

    for cid in cell_ids:
        cid = int(cid)
        vol = all_cells[cid].volume
        if vol is None or vol <= 0.0:
            raise ValueError(f"Cell {cid} has no valid volume (got {vol}).")
        scaling[cid] = neutron_source_rate / float(vol)
    return scaling

def get_flux_spectrum_arrays(t_flux: openmc.Tally) -> tuple[np.ndarray, np.ndarray]:
    """
    Return flux_mean, flux_std with shape:
        (n_cells, n_particles, n_energy)

    OpenMC may return extra trailing singleton axes, e.g.
        (n_cells, n_particles, n_energy, 1, 1)
    so remove them here.
    """
    flux_mean = np.asarray(t_flux.get_reshaped_data(value="mean"), dtype=float)
    flux_std = np.asarray(t_flux.get_reshaped_data(value="std_dev"), dtype=float)

    while flux_mean.ndim > 3 and flux_mean.shape[-1] == 1:
        flux_mean = flux_mean[..., 0]
    while flux_std.ndim > 3 and flux_std.shape[-1] == 1:
        flux_std = flux_std[..., 0]

    if flux_mean.ndim != 3 or flux_std.ndim != 3:
        raise ValueError(
            f"Unexpected flux tally shape after cleanup: "
            f"mean={flux_mean.shape}, std={flux_std.shape}"
        )

    return flux_mean, flux_std

def set_ylim_and_ticks(
    ax,
    ydata,
    ratio=0.3,
    scale="auto",
    eps=1e-99,
    n_linear_ticks=6,
    log_threshold_decades=1.0,
    min_log_pad_decades=0.10,
    ):
    ratio = float(ratio)
    ratio = max(0.0, min(ratio, 0.95))

    y = np.asarray(ydata, dtype=float)
    y = y[np.isfinite(y)]
    if y.size == 0:
        return

    if scale == "auto":
        ypos = y[y > 0]
        if ypos.size < 2:
            scale_eff = "linear"
        else:
            lo = np.log10(np.min(ypos))
            hi = np.log10(np.max(ypos))
            scale_eff = "log" if (hi - lo) >= float(log_threshold_decades) else "linear"
    else:
        scale_eff = scale

    if scale_eff == "log":
        ypos = y[y > 0]
        if ypos.size == 0:
            scale_eff = "linear"
        else:
            logy = np.log10(ypos)
            lo = float(logy.min())
            hi = float(logy.max())
            span = max(hi - lo, 1e-12)

            pad = max(ratio * span, float(min_log_pad_decades))
            ymin = max(10 ** (lo - pad), eps)
            ymax = 10 ** (hi + pad)
            ax.set_ylim(ymin, ymax)

            ax.yaxis.set_major_locator(mticker.LogLocator(base=10.0, subs=(1.0, 2.0, 5.0)))
            ax.yaxis.set_major_formatter(mticker.LogFormatterMathtext(base=10.0, labelOnlyBase=False))
            ax.yaxis.set_minor_locator(mticker.LogLocator(base=10.0, subs=np.arange(1, 10) * 0.1))
            ax.yaxis.set_minor_formatter(mticker.NullFormatter())
            return

    ymin = float(np.min(y))
    ymax = float(np.max(y))
    if ymin == ymax:
        if ymin == 0.0:
            ymin, ymax = -1.0, 1.0
        else:
            ymin = ymin * (1.0 - ratio)
            ymax = ymax * (1.0 + ratio)
    else:
        span = ymax - ymin
        ymin -= ratio * span
        ymax += ratio * span

    ax.set_ylim(ymin, ymax)
    ax.yaxis.set_major_locator(mticker.MaxNLocator(n_linear_ticks))
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("{x:g}"))

# =============================================================================
# structural-origin fractions per cell/nuclide
# =============================================================================
def require_struct_maps(bm):
    need = [
        "cell_struct_origin_frac",
        "cell_struct_nuclide_atoms",
        "cell_total_atoms_struct",
    ]
    missing = [k for k in need if not hasattr(bm, k)]
    if missing:
        raise RuntimeError(
            f"Missing {missing} in neutronics_model module. "
            "Expected you to compute these in neutronics_model.py "
            "via build_structural_maps_vo/build_structural_maps."
        )

def write_struct_origin_csv(
    outdir: Path,
    chunk_key: str,
    cell_ids: list[int],
    labels: list[str],
    cell_struct_origin_frac: Dict[int, Dict[str, float]],
    ):
    rows = []
    for cid, lab in zip(cell_ids, labels):
        cid = int(cid)
        frac_map = cell_struct_origin_frac.get(cid, {}) or {}
        for nuc, f in sorted(frac_map.items()):
            rows.append(
                dict(
                    chunk_key=chunk_key,
                    cell_id=cid,
                    layer_label=lab,
                    nuclide=str(nuc),
                    f_struct_origin=float(f),
                )
            )
    df = pd.DataFrame(rows)
    df.to_csv(outdir / f"struct_origin_fractions_{chunk_key}.csv", index=False)

def tally_mean_std_for_nuclide(t, *, score: str, nuclide: str) -> tuple[float, float]:
    mean = float(t.get_values(scores=[score], nuclides=[nuclide], value="mean").ravel()[0])
    std = float(t.get_values(scores=[score], nuclides=[nuclide], value="std_dev").ravel()[0])
    return mean, std

def corrected_sum_mean_std_getvalues(
    t: openmc.Tally,
    *,
    score: str,
    nuclides: list[str],
    f_struct: Dict[str, float],
    ) -> tuple[float, float]:
    mean_tot = 0.0
    var_tot = 0.0

    for nuc in nuclides:
        nuc = str(nuc)
        w = float(f_struct.get(nuc, 0.0))
        if w == 0.0:
            continue

        try:
            m, sd = tally_mean_std_for_nuclide(t, score=score, nuclide=nuc)
        except Exception:
            continue

        mean_tot += w * m
        var_tot += (w * sd) ** 2

    return mean_tot, math.sqrt(max(var_tot, 0.0))

# =============================================================================
# Load chunking + bin helpers from neutronics_model.py
# =============================================================================
_cells_chunk = bm.build_breeder_chunks(
    INPUT_JSON,
    default_chunk_key="OB_1_b6",
    gap_cm=2.0,
    start_cm=0.0,
)

data = _cells_chunk["data"]
geom = _cells_chunk["geom"]

cell_ids_all = _cells_chunk["cell_ids_all"]
ob_by_key = _cells_chunk["ob_by_key"]
ib_by_key = _cells_chunk["ib_by_key"]
n_breeder = int(geom["n_breeder"])

cell_ids_for_key = _cells_chunk["cell_ids_for_key"]
radial_bins_for_key = _cells_chunk["radial_bins_for_key"]

def generate_colors(n):
    """Generates a smooth rainbow gradient of n RGB colors."""
    nc = n + 1
    cmap = plt.get_cmap('rainbow')
    color_range = cmap(np.linspace(1, 0, nc))
    return color_range

ob_n_layers = int(geom["ob_n_layers"])
colors = generate_colors(ob_n_layers + 3)

def make_chunk_base_df(
    *,
    chunk_key: str,
    cell_ids: list[int],
    labels: list[str],
    xcent: np.ndarray,
    xedges: np.ndarray,
    extra_cols: dict[str, object] | None = None,
    ) -> pd.DataFrame:
    cell_ids = [int(c) for c in cell_ids]
    xcent = np.asarray(xcent, float).ravel()
    xedges = np.asarray(xedges, float).ravel()

    n = len(cell_ids)
    if len(labels) != n:
        raise ValueError(f"{chunk_key}: labels length {len(labels)} != cell_ids length {n}")
    if len(xcent) != n:
        raise ValueError(f"{chunk_key}: len(xcent) {len(xcent)} != N {n}")
    if len(xedges) != n + 1:
        raise ValueError(f"{chunk_key}: len(xedges) {len(xedges)} != N+1 {n+1}")

    df = pd.DataFrame(
        {
            "chunk_key": [chunk_key] * n,
            "cell_id": cell_ids,
            "layer_label": list(labels),
            "x_center_cm": xcent,
            "x_left_cm": xedges[:-1],
            "x_right_cm": xedges[1:],
        }
    )

    if extra_cols:
        for k, v in extra_cols.items():
            if np.isscalar(v) or isinstance(v, str):
                df[k] = v
            else:
                vv = np.asarray(v)
                if vv.shape[0] != n:
                    raise ValueError(f"{chunk_key}: extra_cols[{k}] has length {vv.shape[0]} != {n}")
                df[k] = vv

    return df

def save_profile_from_base(
    base_df: pd.DataFrame,
    outdir: Path,
    *,
    quantity: str,
    mean: np.ndarray,
    std: np.ndarray,
    units: str,
    nonnegative_lower: bool = True,
    ) -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    n = len(base_df)
    mean = np.asarray(mean, float).ravel()
    std = np.asarray(std, float).ravel()
    if len(mean) != n or len(std) != n:
        raise ValueError(f"{quantity}: mean/std length mismatch with base_df ({len(mean)},{len(std)}) vs {n}")

    df = base_df.copy()
    df["quantity"] = quantity
    df["units"] = units
    df["mean"] = mean
    df["std"] = std

    lower = mean - std
    if nonnegative_lower:
        lower = np.maximum(lower, 0.0)

    df["lower_1sigma"] = lower
    df["upper_1sigma"] = mean + std

    if "x_left_cm" in df.columns:
        df = df.sort_values(["x_left_cm", "cell_id"]).reset_index(drop=True)

    chunk_key = str(df["chunk_key"].iloc[0])
    csv_path = outdir / f"profile_{quantity}_{chunk_key}.csv"
    df.to_csv(csv_path, index=False)
    return csv_path

# =============================================================================
# build tally maps from the statepoint
# =============================================================================
DPA_SCORES = {"damage-energy"}
GAS_SCORES_EXPLICIT = {"H1-production", "He3-production", "He4-production"}
GAS_SCORES_REACTION = {"(n,Xp)", "(n,X3He)", "(n,Xa)"}

def build_dpa_gas_map(sp: openmc.StatePoint) -> Dict[int, openmc.Tally]:
    out: Dict[int, openmc.Tally] = {}
    for t in sp.tallies.values():
        scores = {str(s) for s in (t.scores or [])}

        if "damage-energy" not in scores:
            continue

        has_explicit = GAS_SCORES_EXPLICIT.issubset(scores)
        has_reaction = GAS_SCORES_REACTION.issubset(scores)
        if not (has_explicit or has_reaction):
            continue

        cid = _get_cell_id_from_tally(t)
        if cid is None:
            continue
        out[int(cid)] = t
    return out

def _get_cell_id_from_tally(t: openmc.Tally) -> Optional[int]:
    cf = next((f for f in t.filters if isinstance(f, openmc.CellFilter)), None)
    if cf is None or not cf.bins:
        return None
    return int(cf.bins[0])

def require_tally_id(sp: openmc.StatePoint, bm_tally_obj: openmc.Tally, desc: str) -> int:
    tid = int(bm_tally_obj.id)
    if tid not in sp.tallies:
        raise RuntimeError(
            f"[fatal] {desc}: bm id={tid} not found in statepoint. "
            f"Your postprocessing module and the statepoint were produced by different tally sets."
        )
    return tid

# =============================================================================
# Per-chunk labels
# =============================================================================
def layer_names_for_chunk(chunk_cells: list[int], region_tag: str) -> list[str]:
    n = len(chunk_cells)
    labels = ["Armor", "First_Wall"]
    n_layers = n - 3
    labels += [f"Breeder Layer {i+1}" for i in range(n_layers)]
    labels += ["Vacuum Vessel"]
    return labels

# =============================================================================
# Albedo calculation
# =============================================================================
def _pydagmc_surface(pydagmc_model, sid: int):
    return pydagmc_model.surfaces_by_id[int(sid)]

def _surface_area(pydagmc_model, sid: int) -> float:
    return float(_pydagmc_surface(pydagmc_model, sid).area)

def _is_surface_unshared(pydagmc_model, sid: int) -> bool:
    s = _pydagmc_surface(pydagmc_model, sid)
    try:
        return len(s.volumes) == 1
    except Exception:
        senses = getattr(s, "senses", [None, None])
        vols = [v for v in senses if v is not None]
        return len(vols) == 1

def _cell_surface_ids_from_info_flat(bm_info_flat: dict, cell_id: int) -> List[int]:
    out = []
    for rec in bm_info_flat.get(int(cell_id), []):
        try:
            out.append(int(rec["surface_id"]))
        except Exception:
            pass
    return out

def _rank_sids_by_area(pydagmc_model, sids: List[int]) -> List[int]:
    pairs = []
    for sid in sids:
        try:
            pairs.append((_surface_area(pydagmc_model, sid), int(sid)))
        except Exception:
            continue
    pairs.sort(key=lambda x: x[0], reverse=True)
    return [sid for _, sid in pairs]

def build_special_surfaces(
    pydagmc_model,
    bm_info_flat: dict,
    *,
    armor_cell_id: int,
    breeder_back_cell_id: int,
    vv_cell_id: int,
    allowed_surface_ids: Set[int],
    external_surface_ids: Set[int],
    internal_surface_ids: Set[int],
) -> Tuple[Dict[str, int], Set[int], Dict[str, int], Dict[str, object]]:
    allowed = set(int(x) for x in allowed_surface_ids)
    ext_set = set(int(x) for x in external_surface_ids)
    int_set = set(int(x) for x in internal_surface_ids)

    special: Dict[str, int] = {}
    owner: Dict[str, int] = {}
    special_ids: Set[int] = set()
    debug: Dict[str, object] = {}

    armor_sids = [sid for sid in _cell_surface_ids_from_info_flat(bm_info_flat, armor_cell_id) if sid in allowed]
    armor_ranked = _rank_sids_by_area(pydagmc_model, armor_sids)

    armor_front = next((sid for sid in armor_ranked if sid in ext_set), None)
    armor_internal_ranked = [sid for sid in armor_ranked if sid in int_set]
    armor_back = armor_internal_ranked[0] if armor_internal_ranked else None

    debug["armor"] = {
        "cell_id": int(armor_cell_id),
        "ranked_top": armor_ranked[:8],
        "picked_front_ext": armor_front,
        "picked_back_int": armor_back,
        "internal_candidates": armor_internal_ranked[:8],
    }

    if armor_front is not None:
        special["Armor_front_ext"] = int(armor_front)
        owner["Armor_front_ext"] = int(armor_cell_id)
        special_ids.add(int(armor_front))

    if armor_back is not None:
        special["Armor_back_int"] = int(armor_back)
        owner["Armor_back_int"] = int(armor_cell_id)
        special_ids.add(int(armor_back))

    breeder_sids = [
        sid
        for sid in _cell_surface_ids_from_info_flat(bm_info_flat, breeder_back_cell_id)
        if sid in allowed
    ]
    breeder_ranked = _rank_sids_by_area(pydagmc_model, breeder_sids)
    breeder_ext = next((sid for sid in breeder_ranked if sid in ext_set), None)

    debug["breeder_back"] = {
        "cell_id": int(breeder_back_cell_id),
        "ranked_top": breeder_ranked[:8],
        "picked_largest_ext": breeder_ext,
    }

    if breeder_ext is not None:
        special["Breeder_back_ext"] = int(breeder_ext)
        owner["Breeder_back_ext"] = int(breeder_back_cell_id)
        special_ids.add(int(breeder_ext))

    vv_sids = [sid for sid in _cell_surface_ids_from_info_flat(bm_info_flat, vv_cell_id) if sid in allowed]
    vv_ranked = _rank_sids_by_area(pydagmc_model, vv_sids)

    vv_picks = []
    vv_rejects = []
    for sid in vv_ranked:
        if _is_surface_unshared(pydagmc_model, sid):
            vv_picks.append(sid)
            if len(vv_picks) == 2:
                break
        else:
            try:
                vols = [v.id for v in _pydagmc_surface(pydagmc_model, sid).volumes]
            except Exception:
                vols = []
            vv_rejects.append((sid, vols))

    debug["vv"] = {
        "cell_id": int(vv_cell_id),
        "ranked_top": vv_ranked[:12],
        "picked_unshared": vv_picks,
        "first_rejects": vv_rejects[:8],
    }

    if len(vv_picks) >= 2:
        special["VV_face_1"] = int(vv_picks[0])
        special["VV_face_2"] = int(vv_picks[1])
        owner["VV_face_1"] = int(vv_cell_id)
        owner["VV_face_2"] = int(vv_cell_id)
        special_ids.add(int(vv_picks[0]))
        special_ids.add(int(vv_picks[1]))

    return special, special_ids, owner, debug

def compute_albedo_for_chunk(
    sp: openmc.StatePoint,
    bm,
    *,
    cell_ids: List[int],
    labels: List[str],
    outdir: Path,
    chunk_key: str,
    pydagmc_model,
    EXCLUDED_SURFACES: Optional[set[int]] = None,
) -> Dict[str, object]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    EXCLUDED_SURFACES = set(EXCLUDED_SURFACES or set())

    required = ("info", "t_current_tally", "p_current_tallies", "external_surface_ids", "internal_surface_ids")
    missing = [x for x in required if not hasattr(bm, x)]
    if missing:
        print(f"[warn] {chunk_key}: skipping albedo (bm missing: {missing})")
        return {}

    if len(labels) != len(cell_ids):
        raise ValueError(f"{chunk_key}: labels length {len(labels)} != cell_ids length {len(cell_ids)}")

    PARTICLES = ("neutron", "photon")
    eps = 1e-15

    info = bm.info
    ext_ids = [int(x) for x in bm.external_surface_ids]
    int_ids = [int(x) for x in bm.internal_surface_ids]
    all_ids = [int(x) for x in getattr(bm, "all_surface_ids", [])] or sorted(set(ext_ids) | set(int_ids))

    ext_set = set(ext_ids)
    int_set = set(int_ids)
    all_set = set(all_ids)
    chunk_cell_set = set(int(x) for x in cell_ids)
    idx_map = {int(cid): i for i, cid in enumerate(cell_ids)}

    def _cell_surfaces(cid: int, which: str) -> List[dict]:
        x = info.get(int(cid), {})
        if isinstance(x, dict):
            return list(x.get(which, []) or [])
        return list(x or [])

    def _read_current_df_map(tally) -> Dict[tuple[int, str], tuple[float, float]]:
        df = tally.get_pandas_dataframe()
        out = {}
        for _, r in df.iterrows():
            sid = int(r["surface"])
            part = str(r["particle"]).lower()
            out[(sid, part)] = (float(r["mean"]), float(r["std. dev."]))
        return out

    surf_to_cells: Dict[int, set[int]] = {}
    for cid in cell_ids:
        cid = int(cid)
        for s in _cell_surfaces(cid, "all_surfaces"):
            try:
                sid = int(s["surface_id"])
            except Exception:
                continue
            surf_to_cells.setdefault(sid, set()).add(cid)

    if not surf_to_cells:
        print(f"[warn] {chunk_key}: no surfaces mapped from bm.info; skipping albedo.")
        return {}

    Jnet_by_sid_particle: Dict[tuple[int, str], tuple[float, float]] = {}
    try:
        t_tot = sp.get_tally(id=bm.t_current_tally.id)
        Jnet_by_sid_particle = _read_current_df_map(t_tot)
    except Exception as e:
        print(f"[warn] {chunk_key}: could not read total current tally; external albedo limited. ({e})")

    allowed_surface_ids = set(all_set)

    special: Dict[str, int] = {}
    special_surface_ids: set[int] = set()
    special_owner_cell: Dict[str, int] = {}
    special_debug: Dict[str, object] = {}

    if len(cell_ids) >= 3:
        armor_cell_id = int(cell_ids[0])
        breeder_back_cell_id = int(cell_ids[-2])
        vv_cell_id = int(cell_ids[-1])

        bm_info_flat = {int(cid): _cell_surfaces(int(cid), "all_surfaces") for cid in cell_ids}

        special, special_surface_ids, special_owner_cell, special_debug = build_special_surfaces(
            pydagmc_model,
            bm_info_flat,
            armor_cell_id=armor_cell_id,
            breeder_back_cell_id=breeder_back_cell_id,
            vv_cell_id=vv_cell_id,
            allowed_surface_ids=allowed_surface_ids,
            external_surface_ids=ext_set,
            internal_surface_ids=int_set,
        )
        print(f"[debug] {chunk_key} specials picked: {special}")
    else:
        print(f"[warn] {chunk_key}: chunk has only {len(cell_ids)} cells; skipping special selection.")

    partial_by_cell_sid_particle: Dict[tuple[int, int, str], tuple[float, float]] = {}

    for cid, tt in bm.p_current_tallies.items():
        cid = int(cid)
        if cid not in chunk_cell_set:
            continue
        try:
            tally = sp.get_tally(id=tt.id)
            d = _read_current_df_map(tally)
            for (sid, part), (m, s) in d.items():
                partial_by_cell_sid_particle[(cid, sid, part)] = (m, s)
        except Exception as e:
            print(f"[warn] {chunk_key}: failed reading partial tally for cell {cid}: {e}")

    special_labels = {
        "Armor_front_ext": "Armor Front",
        "Armor_back_int": "Armor Back",
        "Breeder_back_ext": "Breeder Back",
        "VV_face_1": "VV Back",
        "VV_face_2": "VV Front",
    }
    special_colors = {
        "Armor_front_ext": "tab:red",
        "Armor_back_int": "tab:orange",
        "Breeder_back_ext": "tab:purple",
        "VV_face_1": "tab:green",
        "VV_face_2": "tab:olive",
    }
    special_markers = {k: "x" for k in special_labels.keys()}
    MIN_VISIBLE_YERR = 1e-4

    results_by_particle: Dict[str, Dict[str, object]] = {}

    for particle in PARTICLES:
        rows: List[dict] = []

        for sid in ext_ids:
            cset = surf_to_cells.get(int(sid), set()) & chunk_cell_set
            if not cset:
                continue

            cid = sorted(cset)[0]

            if (sid, particle) not in Jnet_by_sid_particle:
                continue

            Jnet_m, Jnet_s = Jnet_by_sid_particle[(sid, particle)]
            Jout_m, Jout_s = partial_by_cell_sid_particle.get((cid, sid, particle), (0.0, 0.0))

            Jnet = abs(Jnet_m)
            Jout = abs(Jout_m)
            Jin = abs(Jout - Jnet)

            if Jout < eps:
                A_mean, A_std = 0.0, 0.0
            else:
                A_mean = Jin / Jout
                dA_dJout = Jnet / (Jout ** 2)
                dA_dJnet = -1.0 / Jout
                A_var = (dA_dJout ** 2) * (Jout_s ** 2) + (dA_dJnet ** 2) * (Jnet_s ** 2)
                A_std = math.sqrt(max(A_var, 0.0))

            rows.append({
                "particle": particle,
                "surface_id": sid,
                "cell_id": cid,
                "kind": "external",
                "J_total_mean": float(Jnet_m),
                "J_total_std": float(Jnet_s),
                "J_out_mean": float(Jout_m),
                "J_out_std": float(Jout_s),
                "albedo_mean": float(A_mean),
                "albedo_std": float(A_std),
            })

        for sid in int_ids:
            cset = surf_to_cells.get(int(sid), set()) & chunk_cell_set
            if len(cset) != 2:
                continue

            c1, c2 = sorted(cset)

            Jout1_m, Jout1_s = partial_by_cell_sid_particle.get((c1, sid, particle), (0.0, 0.0))
            Jout2_m, Jout2_s = partial_by_cell_sid_particle.get((c2, sid, particle), (0.0, 0.0))

            J1 = abs(Jout1_m)
            J2 = abs(Jout2_m)
            s1 = abs(Jout1_s)
            s2 = abs(Jout2_s)

            if J1 <= eps or J2 <= eps:
                A1 = np.nan
                A2 = np.nan
                A1_std = np.nan
                A2_std = np.nan
            else:
                A1 = J2 / J1
                dA1_dJ2 = 1.0 / J1
                dA1_dJ1 = -J2 / (J1 ** 2)
                A1_var = (dA1_dJ2 ** 2) * (s2 ** 2) + (dA1_dJ1 ** 2) * (s1 ** 2)
                A1_std = math.sqrt(max(A1_var, 0.0))

                A2 = J1 / J2
                dA2_dJ1 = 1.0 / J2
                dA2_dJ2 = -J1 / (J2 ** 2)
                A2_var = (dA2_dJ1 ** 2) * (s1 ** 2) + (dA2_dJ2 ** 2) * (s2 ** 2)
                A2_std = math.sqrt(max(A2_var, 0.0))

            rows.append({
                "particle": particle,
                "surface_id": sid,
                "cell_id": c1,
                "kind": "internal",
                "J_total_mean": np.nan,
                "J_total_std": np.nan,
                "J_out_mean": float(Jout1_m),
                "J_out_std": float(Jout1_s),
                "albedo_mean": float(A1) if np.isfinite(A1) else np.nan,
                "albedo_std": float(A1_std) if np.isfinite(A1_std) else np.nan,
            })
            rows.append({
                "particle": particle,
                "surface_id": sid,
                "cell_id": c2,
                "kind": "internal",
                "J_total_mean": np.nan,
                "J_total_std": np.nan,
                "J_out_mean": float(Jout2_m),
                "J_out_std": float(Jout2_s),
                "albedo_mean": float(A2) if np.isfinite(A2) else np.nan,
                "albedo_std": float(A2_std) if np.isfinite(A2_std) else np.nan,
            })

        if not rows:
            print(f"[warn] {chunk_key}: no albedo rows created for {particle}; skipping outputs.")
            continue

        df_alb = pd.DataFrame(rows).sort_values(["kind", "surface_id", "cell_id"])
        df_alb.to_csv(outdir / f"surface_albedo_all_{particle}_{chunk_key}.csv", index=False)

        special_rows = []
        for key, sid in special.items():
            want_cid = int(special_owner_cell.get(key, -1))
            cand = df_alb[df_alb["surface_id"] == int(sid)]
            if want_cid != -1:
                cand2 = cand[cand["cell_id"] == want_cid]
                if not cand2.empty:
                    cand = cand2
            if cand.empty:
                continue
            r = cand.iloc[0].to_dict()
            r["special_key"] = key
            special_rows.append(r)

        df_special = (
            pd.DataFrame(special_rows)
            if special_rows
            else pd.DataFrame(columns=list(df_alb.columns) + ["special_key"])
        )
        df_special.to_csv(outdir / f"surface_albedo_special_{particle}_{chunk_key}.csv", index=False)

        incident_records = []
        for _, r in df_special.iterrows():
            key = str(r.get("special_key", ""))
            if key != "Armor_front_ext":
                continue
            J_out = abs(float(r["J_out_mean"]))
            J_net = abs(float(r["J_total_mean"])) if pd.notna(r["J_total_mean"]) else 0.0
            J_in  = abs(J_out - J_net)
            J_out_std = abs(float(r["J_out_std"]))
            J_net_std = abs(float(r["J_total_std"])) if pd.notna(r["J_total_std"]) else 0.0
            J_in_std  = math.sqrt(J_out_std**2 + J_net_std**2)

            incident_records.append({
                "particle":      particle,
                "special_key":   key,
                "surface_id":    int(r["surface_id"]),
                "cell_id":       int(r["cell_id"]),
                "J_out_mean":    J_out,
                "J_out_std":     J_out_std,
                "J_net_mean":    float(r["J_total_mean"]) if pd.notna(r["J_total_mean"]) else None,
                "J_net_std":     J_net_std,
                "J_in_mean":     J_in,
                "J_in_std":      J_in_std,
                "chunk_key":     chunk_key,
            })

        if incident_records:
            out_path = outdir / f"armor_current_{particle}.json"
            with open(out_path, "w") as f:
                json.dump(incident_records, f, indent=2)
            print(f"[saved] {out_path}")

        df_external_rest = df_alb[
            (df_alb["kind"] == "external")
            & (~df_alb["surface_id"].isin(EXCLUDED_SURFACES))
            & (~df_alb["surface_id"].isin(special_surface_ids))
        ].copy()
        df_external_rest.to_csv(outdir / f"surface_albedo_external_filtered_{particle}_{chunk_key}.csv", index=False)

        cell_rows = []
        for cid, grp in df_external_rest.groupby("cell_id"):
            cid = int(cid)
            A_mean = float(grp["albedo_mean"].mean()) if len(grp) else float("nan")
            st = grp["albedo_std"].to_numpy(dtype=float)
            A_std = float(np.sqrt(np.nansum(st * st)) / max(len(st), 1)) if len(st) else float("nan")
            layer_label = labels[idx_map[cid]] if cid in idx_map else f"cell_{cid}"

            cell_rows.append({
                "particle": particle,
                "cell_id": cid,
                "layer_label": layer_label,
                "albedo_mean": A_mean,
                "albedo_std": A_std,
                "n_surfaces": int(len(grp)),
            })

        df_cell = (
            pd.DataFrame(cell_rows).sort_values("cell_id")
            if cell_rows
            else pd.DataFrame(columns=["particle", "cell_id", "layer_label", "albedo_mean", "albedo_std", "n_surfaces"])
        )
        df_cell.to_csv(outdir / f"cell_albedo_summary_external_{particle}_{chunk_key}.csv", index=False)

        xpos = np.arange(len(labels))
        y = np.full(len(labels), np.nan, dtype=float)
        e = np.full(len(labels), np.nan, dtype=float)

        for _, row in df_cell.iterrows():
            cid = int(row["cell_id"])
            if cid not in idx_map:
                continue
            i = idx_map[cid]
            y[i] = float(row["albedo_mean"]) if pd.notna(row["albedo_mean"]) else np.nan
            e[i] = float(row["albedo_std"]) if pd.notna(row["albedo_std"]) else np.nan

        plt.figure(figsize=(11, 6))
        plt.errorbar(xpos, y, yerr=e, fmt="none", elinewidth=1, capsize=4)
        plt.scatter(
            xpos, y, marker="x", s=90, linewidths=2,
            label=f"{particle.capitalize()} cell avg (external, excluding specials)"
        )

        used_labels = set()
        if not df_special.empty:
            for _, r in df_special.iterrows():
                key = str(r.get("special_key", "special"))
                cid = int(r["cell_id"])
                if cid not in idx_map:
                    continue

                x0 = idx_map[cid]
                y0 = float(r["albedo_mean"]) if pd.notna(r["albedo_mean"]) else np.nan
                e0 = float(r["albedo_std"]) if pd.notna(r["albedo_std"]) else np.nan
                if not np.isfinite(y0):
                    continue

                color = special_colors.get(key, "black")
                label = special_labels.get(key, key)
                marker = special_markers.get(key, "x")
                plot_label = None if label in used_labels else label
                used_labels.add(label)

                if np.isfinite(e0) and e0 > 0:
                    e_vis = max(e0, MIN_VISIBLE_YERR)
                    plt.errorbar(
                        [x0], [y0],
                        yerr=[[e_vis], [e_vis]],
                        fmt="none",
                        ecolor=color,
                        elinewidth=2,
                        capsize=6,
                        capthick=2,
                        zorder=5,
                    )

                plt.scatter([x0], [y0], marker=marker, s=120, color=color, zorder=6, label=plot_label)

        plt.xticks(xpos, labels, rotation=45, fontsize=11)
        plt.xlabel("Layer", fontsize=12)
        plt.ylabel("Albedo", fontsize=12)
        plt.title(f"{particle.capitalize()} Albedo per Layer: {chunk_key}", fontsize=15)
        plt.grid(axis="y", linestyle="--", alpha=0.6)
        plt.tight_layout()
        plt.legend()
        plt.savefig(outdir / f"albedo_layers_{particle}_{chunk_key}.png", dpi=300)
        plt.close()

        results_by_particle[particle] = {
            "df_albedo": df_alb,
            "df_albedo_external_filtered": df_external_rest,
            "df_special": df_special,
            "df_cell_summary": df_cell,
        }

    if not results_by_particle:
        print(f"[warn] {chunk_key}: no particle-specific albedo outputs were created.")
        return {}

    print(f"[ok] {chunk_key}: albedo done for particles = {list(results_by_particle.keys())}")

    return {
        "by_particle": results_by_particle,
        "special": special,
        "special_surface_ids": special_surface_ids,
        "special_owner_cell": special_owner_cell,
        "special_debug": special_debug,
        "excluded_surfaces": EXCLUDED_SURFACES,
    }

# =============================================================================
# DPA helper
# =============================================================================
def element_from_nuclide(nuc: str) -> str:
    m = re.match(r"[A-Za-z]+", nuc)
    return m.group(0) if m else nuc

# =============================================================================
# Core per-chunk
# =============================================================================
def process_chunk(
    sp: openmc.StatePoint,
    chunk_key: str,
    cell_ids: list[int],
    pydagmc_model,
    *,
    dpa_gas_map: Dict[int, openmc.Tally],
    xcentroids,
    xedges,
):
    outdir = RESULTS_DIR / chunk_key
    outdir.mkdir(parents=True, exist_ok=True)

    is_ob = chunk_key.startswith("OB_")
    region_tag = "OB" if is_ob else "IB"

    xcent = np.asarray(xcentroids, float)
    xedges = np.asarray(xedges, float)

    labels = layer_names_for_chunk(cell_ids, region_tag)
    scaling = scaling_for_cells(cell_ids)

    base_df = make_chunk_base_df(
        chunk_key=chunk_key,
        cell_ids=cell_ids,
        labels=labels,
        xcent=xcent,
        xedges=xedges,
    )

    cell_struct_origin_frac = bm.cell_struct_origin_frac
    write_struct_origin_csv(outdir, chunk_key, cell_ids, labels, cell_struct_origin_frac)

    # ==========================================
    # Neutron & photon spectra (flux_tally)
    # ==========================================
    t_flux_id = require_tally_id(sp, bm.flux_tally, "flux spectrum tally")
    t_flux = sp.get_tally(id=t_flux_id)
    cell_bins_flux = get_cell_bins_from_tally(t_flux)

    flux_mean, flux_std = get_flux_spectrum_arrays(t_flux)

    neutron_flux = flux_mean[:, 0, :]
    neutron_flux_std = flux_std[:, 0, :]
    photon_flux = flux_mean[:, 1, :]
    photon_flux_std = flux_std[:, 1, :]

    neutron_flux_chunk = subset_by_cells(neutron_flux, cell_bins_flux, cell_ids)
    neutron_flux_std_chunk = subset_by_cells(neutron_flux_std, cell_bins_flux, cell_ids)
    photon_flux_chunk = subset_by_cells(photon_flux, cell_bins_flux, cell_ids)
    photon_flux_std_chunk = subset_by_cells(photon_flux_std, cell_bins_flux, cell_ids)

    plt.figure()
    for i, cid in enumerate(cell_ids):
        flux_scaled = neutron_flux_chunk[i].flatten() * scaling[int(cid)] / unit_lethargy
        plt.loglog(energies[:-1], flux_scaled, label=f"{labels[i]} ({xcent[i]:.1f} cm)", color=colors[i])
    plt.legend(fontsize=8, ncol=2)
    plt.grid(True, which="both")
    plt.ylabel("Neutron flux per unit lethargy [1/cm$^2$/s]")
    plt.xlabel("Energy [eV]")
    plt.xlim([1e-3, 100e6])
    plt.ylim([1e6, 1e17])
    plt.savefig(outdir / f"n_flux_spectrum_{chunk_key}.png", dpi=300, bbox_inches="tight")
    plt.close()

    plt.figure()
    for i, cid in enumerate(cell_ids):
        flux_scaled = photon_flux_chunk[i].flatten() * scaling[int(cid)] / unit_lethargy
        plt.loglog(energies[:-1], flux_scaled, label=f"{labels[i]} ({xcent[i]:.1f} cm)", color=colors[i])
    plt.legend(fontsize=8, loc="lower left")
    plt.grid(True, which="both")
    plt.ylabel("Photon flux per unit lethargy [1/cm$^2$/s]")
    plt.xlabel("Energy [eV]")
    plt.xlim([1e3, 100e6])
    plt.ylim([1e6, 1e17])
    plt.savefig(outdir / f"p_flux_spectrum_{chunk_key}.png", dpi=300, bbox_inches="tight")
    plt.close()

    neutron_df = pd.DataFrame({"E_mid_eV": E_mid})
    for i, lab in enumerate(labels):
        s = scaling[int(cell_ids[i])]
        neutron_df[lab] = neutron_flux_chunk[i].flatten() * s / unit_lethargy
        neutron_df[f"{lab}_std"] = neutron_flux_std_chunk[i].flatten() * s / unit_lethargy
    neutron_df.to_csv(outdir / f"neutron_spectrum_{chunk_key}.csv", index=False)

    photon_df = pd.DataFrame({"E_mid_eV": E_mid})
    for i, lab in enumerate(labels):
        s = scaling[int(cell_ids[i])]
        photon_df[lab] = photon_flux_chunk[i].flatten() * s / unit_lethargy
        photon_df[f"{lab}_std"] = photon_flux_std_chunk[i].flatten() * s / unit_lethargy
    photon_df.to_csv(outdir / f"photon_spectrum_{chunk_key}.csv", index=False)

    # ==========================================
    # Total flux reconstructed from spectrum
    # ==========================================
    direct_total_neut = np.array(
        [np.sum(neutron_flux_chunk[i]) * scaling[int(cid)] for i, cid in enumerate(cell_ids)],
        dtype=float,
    )
    direct_total_neut_std = np.array(
        [np.sqrt(np.sum(neutron_flux_std_chunk[i] ** 2)) * scaling[int(cid)] for i, cid in enumerate(cell_ids)],
        dtype=float,
    )

    direct_total_phot = np.array(
        [np.sum(photon_flux_chunk[i]) * scaling[int(cid)] for i, cid in enumerate(cell_ids)],
        dtype=float,
    )
    direct_total_phot_std = np.array(
        [np.sqrt(np.sum(photon_flux_std_chunk[i] ** 2)) * scaling[int(cid)] for i, cid in enumerate(cell_ids)],
        dtype=float,
    )

    fig, ax = plt.subplots()
    neut_lower = np.maximum(direct_total_neut - direct_total_neut_std, 1e-30)
    neut_upper = direct_total_neut + direct_total_neut_std
    phot_lower = np.maximum(direct_total_phot - direct_total_phot_std, 1e-30)
    phot_upper = direct_total_phot + direct_total_phot_std

    ax.step(xedges, np.r_[direct_total_neut, direct_total_neut[-1]], where="post", label="neutron")
    ax.fill_between(
        xedges,
        np.r_[neut_lower, neut_lower[-1]],
        np.r_[neut_upper, neut_upper[-1]],
        step="post",
        alpha=0.15,
    )

    ax.step(xedges, np.r_[direct_total_phot, direct_total_phot[-1]], where="post", label="photon")
    ax.fill_between(
        xedges,
        np.r_[phot_lower, phot_lower[-1]],
        np.r_[phot_upper, phot_upper[-1]],
        step="post",
        alpha=0.15,
    )

    ax.set_yscale("log")
    ax.set_ylabel("Total Flux [1/cm$^2$/s]")
    ax.set_xlabel("Radial Position [cm]")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.legend()
    fig.savefig(outdir / f"flux_total_{chunk_key}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    save_profile_from_base(
        base_df, outdir, quantity="flux_total_neutron",
        mean=direct_total_neut, std=direct_total_neut_std, units="1/cm2/s",
    )
    save_profile_from_base(
        base_df, outdir, quantity="flux_total_photon",
        mean=direct_total_phot, std=direct_total_phot_std, units="1/cm2/s",
    )

    # ==========================================
    # Heating
    # ==========================================
    t_heat_id = require_tally_id(sp, bm.heating_tally, "heating tally")
    t_heat = sp.get_tally(id=t_heat_id)
    cell_bins_heat = get_cell_bins_from_tally(t_heat)

    heat_mean = t_heat.get_values(value="mean").flatten()
    heat_std = t_heat.get_values(value="std_dev").flatten()

    heat_mean_chunk = subset_by_cells(heat_mean, cell_bins_heat, cell_ids)
    heat_std_chunk = subset_by_cells(heat_std, cell_bins_heat, cell_ids)

    heat_w = np.array(
        [heat_mean_chunk[i] * scaling[int(cid)] * ev_to_joule for i, cid in enumerate(cell_ids)],
        dtype=float,
    )
    heat_w_std = np.array(
        [heat_std_chunk[i] * scaling[int(cid)] * ev_to_joule for i, cid in enumerate(cell_ids)],
        dtype=float,
    )

    lower_h = np.maximum(heat_w - heat_w_std, 1e-30)
    upper_h = heat_w + heat_w_std

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.step(xedges, np.r_[heat_w, heat_w[-1]], where="post", label="Heating")
    ax.fill_between(xedges, np.r_[lower_h, lower_h[-1]], np.r_[upper_h, upper_h[-1]], step="post", alpha=0.3)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.set_ylabel("Heating [W/cm³]")
    ax.set_xlabel("Radial Position [cm]")
    fig.savefig(outdir / f"heating_{chunk_key}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    save_profile_from_base(base_df, outdir, quantity="heating", mean=heat_w, std=heat_w_std, units="W/cm3")

    # ==========================================
    # (n,gamma) reaction rate
    # ==========================================
    if hasattr(bm, "ngamma_tally"): # Not mandatory tally
        t_ngamma_id = require_tally_id(sp, bm.ngamma_tally, "(n,gamma) tally")
        t_ngamma = sp.get_tally(id=t_ngamma_id)
        cell_bins_ngamma = get_cell_bins_from_tally(t_ngamma)

        ngamma_mean = t_ngamma.get_values(value="mean").flatten()
        ngamma_std  = t_ngamma.get_values(value="std_dev").flatten()

        ngamma_mean_chunk = subset_by_cells(ngamma_mean, cell_bins_ngamma, cell_ids)
        ngamma_std_chunk  = subset_by_cells(ngamma_std,  cell_bins_ngamma, cell_ids)

        # Scale to reaction rate
        ngamma_rate = np.array(
            [ngamma_mean_chunk[i] * scaling[int(cid)] for i, cid in enumerate(cell_ids)],
            dtype=float,
        )
        ngamma_rate_std = np.array(
            [ngamma_std_chunk[i] * scaling[int(cid)] for i, cid in enumerate(cell_ids)],
            dtype=float,
        )

        lower_ng = np.maximum(ngamma_rate - ngamma_rate_std, 1e-30)
        upper_ng = ngamma_rate + ngamma_rate_std

        fig, ax = plt.subplots()
        ax.set_yscale("log")
        ax.step(xedges, np.r_[ngamma_rate, ngamma_rate[-1]], where="post", label="(n,γ)")
        ax.fill_between(
            xedges,
            np.r_[lower_ng, lower_ng[-1]],
            np.r_[upper_ng, upper_ng[-1]],
            step="post",
            alpha=0.3,
        )
        ax.grid(True, which="both", linestyle="--", linewidth=0.5)
        ax.set_ylabel("(n,γ) reaction rate [reactions/cm³/s]")
        ax.set_xlabel("Radial Position [cm]")
        ax.set_title(f"(n,γ) reaction rate — {chunk_key}")
        ax.legend()
        fig.savefig(outdir / f"ngamma_{chunk_key}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

        save_profile_from_base(
            base_df, outdir,
            quantity="ngamma_reaction_rate",
            mean=ngamma_rate, std=ngamma_rate_std, units="reactions/cm3/s",
        )
    else:
        print(f"[warn] {chunk_key}: bm.ngamma_tally not found, skipping (n,gamma) plot.")

    # ==========================================
    # H production
    # ==========================================
    h_appm_y = np.zeros(len(cell_ids), dtype=float)
    h_appm_y_std = np.zeros(len(cell_ids), dtype=float)

    for i, cid in enumerate(cell_ids):
        cid = int(cid)
        t = dpa_gas_map.get(cid, None)
        if t is None:
            continue

        scores = {str(s) for s in (t.scores or [])}
        if GAS_SCORES_EXPLICIT.issubset(scores):
            h1_score = "H1-production"
            h2_score = "H2-production"
            h3_score = "H3-production"
        elif GAS_SCORES_REACTION.issubset(scores):
            h1_score = "(n,Xp)"
            h2_score = "(n,Xd)"
            h3_score = "(n,Xt)"
        else:
            continue

        f_struct = cell_struct_origin_frac.get(cid, {}) or {}
        nuclides = list(t.nuclides or [])

        h1_mean_corr, h1_std_corr = corrected_sum_mean_std_getvalues(
            t, score=h1_score, nuclides=nuclides, f_struct=f_struct
        )
        h2_mean_corr, h2_std_corr = corrected_sum_mean_std_getvalues(
            t, score=h2_score, nuclides=nuclides, f_struct=f_struct
        )
        h3_mean_corr, h3_std_corr = corrected_sum_mean_std_getvalues(
            t, score=h3_score, nuclides=nuclides, f_struct=f_struct
        )

        h_mean_corr = h1_mean_corr + h2_mean_corr + h3_mean_corr
        h_std_corr = math.sqrt(h1_std_corr ** 2 + h2_std_corr ** 2 + h3_std_corr ** 2)

        denom = float(bm.cell_total_atoms_struct.get(cid, 0.0))
        if denom > 0.0:
            factor = neutron_source_rate * s_in_y / denom * 1e6
            h_appm_y[i] = h_mean_corr * factor
            h_appm_y_std[i] = h_std_corr * factor

    lower = np.maximum(h_appm_y - h_appm_y_std, 1e-30)
    upper = h_appm_y + h_appm_y_std

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.step(xedges, np.r_[h_appm_y, h_appm_y[-1]], where="post")
    ax.fill_between(xedges, np.r_[lower, lower[-1]], np.r_[upper, upper[-1]], step="post", alpha=0.3)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.set_ylabel("H [appm/fpy]")
    ax.set_xlabel("Radial Position [cm]")
    fig.savefig(outdir / f"h1_{chunk_key}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    save_profile_from_base(
        base_df, outdir,
        quantity="H_appm_fpy_struct_origin",
        mean=h_appm_y, std=h_appm_y_std, units="appm/fpy",
    )

    # ==========================================
    # He production
    # ==========================================
    he_appm_y = np.zeros(len(cell_ids), dtype=float)
    he_appm_y_std = np.zeros(len(cell_ids), dtype=float)

    for i, cid in enumerate(cell_ids):
        cid = int(cid)
        t = dpa_gas_map.get(cid, None)
        if t is None:
            continue

        scores = {str(s) for s in (t.scores or [])}
        if GAS_SCORES_EXPLICIT.issubset(scores):
            he3_score = "He3-production"
            he4_score = "He4-production"
        elif GAS_SCORES_REACTION.issubset(scores):
            he3_score = "(n,X3He)"
            he4_score = "(n,Xa)"
        else:
            continue

        f_struct = cell_struct_origin_frac.get(cid, {}) or {}
        nuclides = list(t.nuclides or [])

        he3_mean_corr, he3_std_corr = corrected_sum_mean_std_getvalues(
            t, score=he3_score, nuclides=nuclides, f_struct=f_struct
        )
        he4_mean_corr, he4_std_corr = corrected_sum_mean_std_getvalues(
            t, score=he4_score, nuclides=nuclides, f_struct=f_struct
        )

        he_mean_corr = he3_mean_corr + he4_mean_corr
        he_std_corr = math.sqrt(he3_std_corr ** 2 + he4_std_corr ** 2)

        denom = float(bm.cell_total_atoms_struct.get(cid, 0.0))
        if denom > 0.0:
            factor = neutron_source_rate * s_in_y / denom * 1e6
            he_appm_y[i] = he_mean_corr * factor
            he_appm_y_std[i] = he_std_corr * factor

    lower = np.maximum(he_appm_y - he_appm_y_std, 1e-30)
    upper = he_appm_y + he_appm_y_std

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.step(xedges, np.r_[he_appm_y, he_appm_y[-1]], where="post")
    ax.fill_between(xedges, np.r_[lower, lower[-1]], np.r_[upper, upper[-1]], step="post", alpha=0.3)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.set_ylabel("He [appm/fpy]")
    ax.set_xlabel("Radial Position [cm]")
    fig.savefig(outdir / f"he_{chunk_key}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    save_profile_from_base(
        base_df, outdir,
        quantity="He_appm_fpy_struct_origin",
        mean=he_appm_y, std=he_appm_y_std, units="appm/fpy",
    )

    # ==========================================
    # DPA
    # ==========================================
    dpa_y = np.zeros(len(cell_ids), dtype=float)
    dpa_y_std = np.zeros(len(cell_ids), dtype=float)

    for i, cid in enumerate(cell_ids):
        cid = int(cid)
        t = dpa_gas_map.get(cid, None)
        if t is None:
            continue

        denom_atoms = float(bm.cell_total_atoms_struct.get(cid, 0.0))
        if denom_atoms <= 0.0:
            continue

        f_struct = cell_struct_origin_frac.get(cid, {}) or {}
        nuclides = list(t.nuclides or [])

        sum_mean = 0.0
        sum_var = 0.0

        for nuc in nuclides:
            nuc = str(nuc)
            dmg_mean, dmg_std = tally_mean_std_for_nuclide(t, score="damage-energy", nuclide=nuc)

            fs = float(f_struct.get(nuc, 0.0))
            dmg_mean *= fs
            dmg_std *= fs

            el = element_from_nuclide(nuc)
            Ed = float(bm.materials.Ed(el))

            disp_scale = 0.8 / (2.0 * Ed)
            source_scale = neutron_source_rate * s_in_y

            disp_mean = disp_scale * source_scale * dmg_mean
            disp_std = disp_scale * source_scale * dmg_std

            dpa_mean_per_nuclide = disp_mean / denom_atoms
            dpa_std_per_nuclide = disp_std / denom_atoms

            sum_mean += dpa_mean_per_nuclide
            sum_var += dpa_std_per_nuclide ** 2

        dpa_y[i] = sum_mean
        dpa_y_std[i] = math.sqrt(sum_var)

    lower = np.maximum(dpa_y - dpa_y_std, 1e-30)
    upper = dpa_y + dpa_y_std

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.step(xedges, np.r_[dpa_y, dpa_y[-1]], where="post")
    ax.fill_between(xedges, np.r_[lower, lower[-1]], np.r_[upper, upper[-1]], step="post", alpha=0.3)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.set_ylabel("NRT-dpa/fpy")
    ax.set_xlabel("Radial Position [cm]")
    fig.savefig(outdir / f"dpa_{chunk_key}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    save_profile_from_base(
        base_df, outdir,
        quantity="dpa_fpy_struct_origin",
        mean=dpa_y, std=dpa_y_std, units="DPA/fpy",
    )

    if DO_ALBEDO:
        compute_albedo_for_chunk(
            sp,
            bm,
            cell_ids=cell_ids,
            labels=labels,
            outdir=outdir,
            chunk_key=chunk_key,
            pydagmc_model=pydagmc_model,
            EXCLUDED_SURFACES=EXCLUDED_SURFACES,
        )

# =============================================================================
# Region-1-only per-layer analysis helpers
# =============================================================================
def layer_cells_first_region(ob_by_key, ib_by_key, layer_index, n_breeder, order="OB_THEN_IB"):
    keys_ob = [f"OB_1_b{b}" for b in range(1, n_breeder + 1)]
    keys_ib = [f"IB_1_b{b}" for b in range(1, n_breeder + 1)]

    if order == "OB_THEN_IB":
        keys = keys_ob + keys_ib
    elif order == "IB_THEN_OB":
        keys = keys_ib + keys_ob
    else:
        raise ValueError("order must be 'OB_THEN_IB' or 'IB_THEN_OB'")

    cell_ids_layer = []
    chunk_labels = []

    for k in keys:
        chunk = ob_by_key.get(k) if k.startswith("OB_") else ib_by_key.get(k)
        if not chunk:
            continue
        if layer_index >= len(chunk):
            continue
        cell_ids_layer.append(int(chunk[layer_index]))
        chunk_labels.append(k)

    return cell_ids_layer, chunk_labels

def compute_dpa_layer_region1_only(
    sp: openmc.StatePoint,
    cell_ids_layer: list[int],
    neutron_source_rate: float,
    s_in_y: float,
    *,
    dpa_gas_map: Dict[int, openmc.Tally],
):
    idx_map = {int(cid): i for i, cid in enumerate(cell_ids_layer)}
    n = len(cell_ids_layer)

    dpa_y = np.zeros(n, dtype=float)
    dpa_y_std = np.zeros(n, dtype=float)

    cell_struct_origin_frac = bm.cell_struct_origin_frac
    source_scale = float(neutron_source_rate) * float(s_in_y)

    for cid in cell_ids_layer:
        cid = int(cid)
        t = dpa_gas_map.get(cid, None)
        if t is None:
            continue

        denom_atoms = float(bm.cell_total_atoms_struct.get(cid, 0.0))
        if denom_atoms <= 0.0:
            continue

        f_struct = cell_struct_origin_frac.get(cid, {}) or {}
        nuclides = list(t.nuclides or [])

        sum_mean = 0.0
        sum_var = 0.0

        for nuc in nuclides:
            nuc = str(nuc)

            try:
                dmg_mean, dmg_std = tally_mean_std_for_nuclide(t, score="damage-energy", nuclide=nuc)
            except Exception:
                continue

            fs = float(f_struct.get(nuc, 0.0))
            if fs == 0.0:
                continue

            dmg_mean *= fs
            dmg_std *= fs

            el = element_from_nuclide(nuc)
            Ed = float(bm.materials.Ed(el))

            disp_scale = 0.8 / (2.0 * Ed)

            disp_mean = disp_scale * source_scale * dmg_mean
            disp_std = disp_scale * source_scale * dmg_std

            dpa_mean = disp_mean / denom_atoms
            dpa_std = disp_std / denom_atoms

            sum_mean += dpa_mean
            sum_var += dpa_std ** 2

        i = idx_map[cid]
        dpa_y[i] = sum_mean
        dpa_y_std[i] = math.sqrt(max(sum_var, 0.0))

    return dpa_y, dpa_y_std

def empty_poloidal_plot(x_data, num_pts):
    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    split = n_breeder - 0.5
    xticks = np.arange(num_pts)

    ax.set_yscale("linear")
    ax.set_xlim([-0.5, len(x_data) - 2.5])
    ax.axvline(split, linestyle="--", linewidth=1, color="k")
    ax.text(
        split - 0.15,
        0.95,
        "$\\leftarrow$ outboard   inboard $\\rightarrow$",
        transform=ax.get_xaxis_transform(),
        ha="center",
        va="top",
        fontsize=10,
    )

    ax.set_xticks(xticks)
    xlabels = [f"OB{i+1}" if i < n_breeder else f"IB{i+1-n_breeder}" for i in range(num_pts)]
    ax.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=9)
    ax.set_xlabel("Poloidal regions", fontsize=11)

    ax.grid(True, which="major", linestyle="--", linewidth=0.6, alpha=0.8)
    ax.yaxis.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.5)
    ax.minorticks_on()
    plt.gca().tick_params(axis="x", which="minor", bottom=False)
    fig.tight_layout()

    return fig, ax

def pad_poloidal_data(data):
    return np.concatenate([[data[0]], data, [data[-1]]])

def process_layer_region1_only(
    sp: openmc.StatePoint,
    layer_index: int,
    layer_tag: str,
    ob_by_key: dict,
    ib_by_key: dict,
    n_breeder: int,
    outdir: Path,
    flux_tally_id: int,
    heating_tally_id: int,
    *,
    dpa_gas_map: Dict[int, openmc.Tally],
    dpa: bool = True,
    hprod: bool = True,
    heprod: bool = True,
    ylim_ratio: float = 0.3,
    yscale_mode: str = "auto",
    log_threshold_decades: float = 1.0,
):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cell_ids_layer, _chunk_labels = layer_cells_first_region(
        ob_by_key, ib_by_key, layer_index=layer_index, n_breeder=n_breeder, order="OB_THEN_IB"
    )
    n = len(cell_ids_layer)
    if n == 0:
        print(f"[warn] No cells found for layer_index={layer_index} (region 1 only). Skipping.")
        return

    scaling = scaling_for_cells(cell_ids_layer)

    x = np.arange(n)
    x = np.concatenate([[x[0] - 1], x, [x[-1] + 1]])

    cell_struct_origin_frac = bm.cell_struct_origin_frac

    def _apply_scale(ax, data, eps_for_log):
        if yscale_mode == "log":
            scale_eff = "log"
        elif yscale_mode == "linear":
            scale_eff = "linear"
        else:
            d = np.asarray(data, dtype=float)
            d = d[np.isfinite(d)]
            dpos = d[d > 0]
            if dpos.size < 2:
                scale_eff = "linear"
            else:
                lo = np.log10(np.min(dpos))
                hi = np.log10(np.max(dpos))
                scale_eff = "log" if (hi - lo) >= float(log_threshold_decades) else "linear"

        ax.set_yscale(scale_eff)
        set_ylim_and_ticks(
            ax,
            data,
            ratio=ylim_ratio,
            scale=scale_eff,
            eps=eps_for_log,
            log_threshold_decades=log_threshold_decades,
        )
        return scale_eff

    # ==========================================
    # TOTAL FLUX reconstructed from spectrum tally
    # ==========================================
    t_flux = sp.get_tally(id=int(flux_tally_id))
    cell_bins_flux = get_cell_bins_from_tally(t_flux)

    flux_mean, flux_std = get_flux_spectrum_arrays(t_flux)

    neutron_flux = flux_mean[:, 0, :]
    neutron_flux_std = flux_std[:, 0, :]
    photon_flux = flux_mean[:, 1, :]
    photon_flux_std = flux_std[:, 1, :]

    neutron_flux_layer = subset_by_cells(neutron_flux, cell_bins_flux, cell_ids_layer)
    neutron_flux_std_layer = subset_by_cells(neutron_flux_std, cell_bins_flux, cell_ids_layer)
    photon_flux_layer = subset_by_cells(photon_flux, cell_bins_flux, cell_ids_layer)
    photon_flux_std_layer = subset_by_cells(photon_flux_std, cell_bins_flux, cell_ids_layer)

    neut = np.zeros(n, dtype=float)
    neut_std = np.zeros(n, dtype=float)
    phot = np.zeros(n, dtype=float)
    phot_std = np.zeros(n, dtype=float)

    for i, cid in enumerate(cell_ids_layer):
        s = scaling[int(cid)]
        neut[i] = np.sum(neutron_flux_layer[i]) * s
        neut_std[i] = np.sqrt(np.sum(neutron_flux_std_layer[i] ** 2)) * s
        phot[i] = np.sum(photon_flux_layer[i]) * s
        phot_std[i] = np.sqrt(np.sum(photon_flux_std_layer[i] ** 2)) * s

    fig, ax = empty_poloidal_plot(x, n)

    eps_flux = 1e-99
    neut = pad_poloidal_data(neut)
    phot = pad_poloidal_data(phot)
    neut_std = pad_poloidal_data(neut_std)
    phot_std = pad_poloidal_data(phot_std)

    neut_plot = np.clip(neut, eps_flux, None)
    phot_plot = np.clip(phot, eps_flux, None)
    neut_low = np.clip(neut - neut_std, eps_flux, None)
    phot_low = np.clip(phot - phot_std, eps_flux, None)
    neut_high = np.clip(neut + neut_std, eps_flux, None)
    phot_high = np.clip(phot + phot_std, eps_flux, None)

    _apply_scale(ax, np.r_[neut_plot, phot_plot], eps_for_log=eps_flux)
    ax.step(x, neut_plot, where="mid", linewidth=2, label="neutron")
    ax.step(x, phot_plot, where="mid", linewidth=2, label="photon")
    ax.fill_between(x, neut_low, neut_high, step="mid", alpha=0.15)
    ax.fill_between(x, phot_low, phot_high, step="mid", alpha=0.15)

    ax.set_ylabel("Total Flux [1/cm$^2$/s]")
    ax.set_title(f"Total Flux — {layer_tag}")
    plt.legend()
    fig.savefig(outdir / f"flux_total_{layer_tag}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ==========================================
    # HEATING
    # ==========================================
    t_heat = sp.get_tally(id=int(heating_tally_id))
    cell_bins_heat = get_cell_bins_from_tally(t_heat)

    heat_mean = t_heat.get_values(value="mean").flatten()
    heat_std = t_heat.get_values(value="std_dev").flatten()
    heat_mean_layer = subset_by_cells(heat_mean, cell_bins_heat, cell_ids_layer)
    heat_std_layer = subset_by_cells(heat_std, cell_bins_heat, cell_ids_layer)

    heat_w = np.zeros(n, dtype=float)
    heat_w_std = np.zeros(n, dtype=float)
    for i, cid in enumerate(cell_ids_layer):
        s = scaling[int(cid)]
        heat_w[i] = heat_mean_layer[i] * s * ev_to_joule
        heat_w_std[i] = heat_std_layer[i] * s * ev_to_joule

    heat_w = pad_poloidal_data(heat_w)
    heat_w_std = pad_poloidal_data(heat_w_std)

    fig, ax = empty_poloidal_plot(x, n)

    eps_heat = 1e-99
    heat_plot = np.clip(heat_w, eps_heat, None)
    heat_low = np.clip(heat_w - heat_w_std, eps_heat, None)
    heat_high = np.clip(heat_w + heat_w_std, eps_heat, None)

    _apply_scale(ax, heat_plot, eps_for_log=eps_heat)
    ax.step(x, heat_plot, where="mid", linewidth=2)
    ax.fill_between(x, heat_low, heat_high, step="mid", alpha=0.15)

    ax.set_ylabel("Heating [W/cm³]")
    ax.set_title(f"Heating — {layer_tag}")
    fig.savefig(outdir / f"heating_{layer_tag}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ==========================================
    # DPA/y
    # ==========================================
    if dpa:
        dpa_y, dpa_y_std = compute_dpa_layer_region1_only(
            sp=sp,
            cell_ids_layer=cell_ids_layer,
            neutron_source_rate=neutron_source_rate,
            s_in_y=s_in_y,
            dpa_gas_map=dpa_gas_map,
        )

        fig, ax = empty_poloidal_plot(x, n)
        dpa_plot = pad_poloidal_data(np.asarray(dpa_y, dtype=float))
        dpa_plot_std = pad_poloidal_data(np.asarray(dpa_y_std, dtype=float))

        dpa_low = np.clip(dpa_plot - dpa_plot_std, eps_heat, None)
        dpa_high = np.clip(dpa_plot + dpa_plot_std, eps_heat, None)

        set_ylim_and_ticks(ax, dpa_plot, ratio=ylim_ratio, scale="linear")
        ax.step(x, dpa_plot, where="mid", linewidth=2)
        ax.fill_between(x, dpa_low, dpa_high, step="mid", alpha=0.15)

        ax.set_ylabel("NRT-dpa/fpy", fontsize=11)
        ax.set_title(f"NRT-dpa/fpy — {layer_tag}", fontsize=12)

        fig.savefig(outdir / f"dpa_{layer_tag}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # ==========================================
    # H production
    # ==========================================
    if hprod:
        h_appm_y = np.zeros(n, dtype=float)
        h_appm_y_std = np.zeros(n, dtype=float)

        for i, cid in enumerate(cell_ids_layer):
            cid = int(cid)
            t = dpa_gas_map.get(cid, None)
            if t is None:
                continue

            scores = {str(s) for s in (t.scores or [])}
            if GAS_SCORES_EXPLICIT.issubset(scores):
                h1_score = "H1-production"
                h2_score = "H2-production"
                h3_score = "H3-production"
            elif GAS_SCORES_REACTION.issubset(scores):
                h1_score = "(n,Xp)"
                h2_score = "(n,Xd)"
                h3_score = "(n,Xt)"
            else:
                continue

            f_struct = cell_struct_origin_frac.get(cid, {}) or {}
            nuclides = list(t.nuclides or [])

            h1_mean_corr, h1_std_corr = corrected_sum_mean_std_getvalues(
                t, score=h1_score, nuclides=nuclides, f_struct=f_struct
            )
            h2_mean_corr, h2_std_corr = corrected_sum_mean_std_getvalues(
                t, score=h2_score, nuclides=nuclides, f_struct=f_struct
            )
            h3_mean_corr, h3_std_corr = corrected_sum_mean_std_getvalues(
                t, score=h3_score, nuclides=nuclides, f_struct=f_struct
            )

            h_mean_corr = h1_mean_corr + h2_mean_corr + h3_mean_corr
            h_std_corr = math.sqrt(h1_std_corr ** 2 + h2_std_corr ** 2 + h3_std_corr ** 2)

            denom = float(bm.cell_total_atoms_struct.get(cid, 0.0))
            if denom > 0.0:
                factor = neutron_source_rate * s_in_y / denom * 1e6
                h_appm_y[i] = h_mean_corr * factor
                h_appm_y_std[i] = h_std_corr * factor

        fig, ax = empty_poloidal_plot(x, n)
        h_appm_y = pad_poloidal_data(h_appm_y)
        h_appm_y_std = pad_poloidal_data(h_appm_y_std)

        eps_h = 1e-30
        h_plot = np.clip(h_appm_y, eps_h, None)
        h_low = np.clip(h_appm_y - h_appm_y_std, eps_h, None)
        h_high = np.clip(h_appm_y + h_appm_y_std, eps_h, None)

        _apply_scale(ax, h_plot, eps_for_log=eps_h)
        ax.step(x, h_plot, where="mid", linewidth=2)
        ax.fill_between(x, h_low, h_high, step="mid", alpha=0.15)

        ax.set_ylabel("H production [appm/fpy]")
        ax.set_title(f"H production — {layer_tag}")
        fig.savefig(outdir / f"h1_{layer_tag}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # ==========================================
    # He production
    # ==========================================
    if heprod:
        he_appm_y = np.zeros(n, dtype=float)
        he_appm_y_std = np.zeros(n, dtype=float)

        for i, cid in enumerate(cell_ids_layer):
            cid = int(cid)
            t = dpa_gas_map.get(cid, None)
            if t is None:
                continue

            scores = {str(s) for s in (t.scores or [])}
            if GAS_SCORES_EXPLICIT.issubset(scores):
                he3_score = "He3-production"
                he4_score = "He4-production"
            elif GAS_SCORES_REACTION.issubset(scores):
                he3_score = "(n,X3He)"
                he4_score = "(n,Xa)"
            else:
                continue

            f_struct = cell_struct_origin_frac.get(cid, {}) or {}
            nuclides = list(t.nuclides or [])
            he3_mean_corr, he3_std_corr = corrected_sum_mean_std_getvalues(
                t, score=he3_score, nuclides=nuclides, f_struct=f_struct
            )
            he4_mean_corr, he4_std_corr = corrected_sum_mean_std_getvalues(
                t, score=he4_score, nuclides=nuclides, f_struct=f_struct
            )

            he_mean_corr = he3_mean_corr + he4_mean_corr
            he_std_corr = math.sqrt(he3_std_corr ** 2 + he4_std_corr ** 2)

            denom = float(bm.cell_total_atoms_struct.get(cid, 0.0))
            if denom > 0.0:
                factor = neutron_source_rate * s_in_y / denom * 1e6
                he_appm_y[i] = he_mean_corr * factor
                he_appm_y_std[i] = he_std_corr * factor

        fig, ax = empty_poloidal_plot(x, n)
        he_appm_y = pad_poloidal_data(he_appm_y)
        he_appm_y_std = pad_poloidal_data(he_appm_y_std)

        eps_he = 1e-30
        he_plot = np.clip(he_appm_y, eps_he, None)
        he_low = np.clip(he_appm_y - he_appm_y_std, eps_he, None)
        he_high = np.clip(he_appm_y + he_appm_y_std, eps_he, None)

        ax.step(x, he_plot, where="mid", linewidth=2)
        ax.fill_between(x, he_low, he_high, step="mid", alpha=0.15)

        _apply_scale(ax, he_plot, eps_for_log=eps_he)

        ax.set_ylabel("He production [appm/fpy]")
        ax.set_title(f"He production — {layer_tag}")
        fig.savefig(outdir / f"he_{layer_tag}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

# =============================================================================
# Run
# =============================================================================
with openmc.StatePoint(str(STATEPOINT_FILE)) as sp:
    dpa_gas_map = build_dpa_gas_map(sp)
    require_struct_maps(bm)
    DO_ALBEDO = True

    for key in KEYS_TO_PROCESS:
        if key in ob_by_key:
            chunk_cells = ob_by_key[key]
        elif key in ib_by_key:
            chunk_cells = ib_by_key[key]
        else:
            print(f"[warn] Unknown chunk key: {key} (skipping)")
            continue

        try:
            xcentroids, widths, x_edges = radial_bins_for_key(key)
        except KeyError:
            print(f"[warn] No radial bins found for key={key} (skipping)")
            continue

        process_chunk(
            sp=sp,
            chunk_key=key,
            cell_ids=chunk_cells,
            pydagmc_model=bm.pydagmc_model,
            dpa_gas_map=dpa_gas_map,
            xcentroids=xcentroids,
            xedges=x_edges,
        )

    for layer_tag, layer_index in layers:
        process_layer_region1_only(
            sp=sp,
            layer_index=layer_index,
            layer_tag=layer_tag,
            ob_by_key=ob_by_key,
            ib_by_key=ib_by_key,
            n_breeder=n_breeder,
            outdir=RESULTS_DIR / f"layer_{layer_tag}_region1_only",
            flux_tally_id=bm.flux_tally.id,
            heating_tally_id=bm.heating_tally.id,
            dpa_gas_map=dpa_gas_map,
            dpa=True,
            hprod=True,
            heprod=True,
            ylim_ratio=0.3,
            yscale_mode="auto",
            log_threshold_decades=1.0,
        )
    
    # ------------------------------------------------------------------
    # Print VV port fill flux tally results
    # ------------------------------------------------------------------
    if bm.test_VV_port_fill:
        t_vv = sp.get_tally(name="flux_tally_VV_port_fill")
        
        flux_mean = t_vv.get_values(scores=["flux"], value="mean").flatten()
        flux_std  = t_vv.get_values(scores=["flux"], value="std_dev").flatten()
        
        # particle_filter bins are ["neutron", "photon"]
        neutron_flux_mean = flux_mean[0]
        neutron_flux_std  = flux_std[0]
        photon_flux_mean  = flux_mean[1]
        photon_flux_std   = flux_std[1]

        # Get cell volume for flux rate conversion
        vv_cell = bm.all_cells[bm.cell_id_VV_port_fill]
        vol = float(vv_cell.volume)

        n_flux_rate = neutron_flux_mean * neutron_source_rate / vol
        p_flux_rate = photon_flux_mean  * neutron_source_rate / vol

        print("\n--- VV Port Fill Flux ---")
        print(f"  Cell ID : {bm.cell_id_VV_port_fill}")
        print(f"  Volume  : {vol:.4e} cm³")
        print(f"  Neutron flux : {n_flux_rate:.4e} +/- {neutron_flux_std * neutron_source_rate / vol:.4e} n/cm²/s")
        print(f"  Photon  flux : {p_flux_rate:.4e} +/- {photon_flux_std  * neutron_source_rate / vol:.4e} γ/cm²/s")

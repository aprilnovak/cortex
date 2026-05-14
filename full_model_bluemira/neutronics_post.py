#!/usr/bin/env python3
"""
neutronics_post.py
==================
Post-process an OpenMC statepoint for both tokamak and slab simulations.

Imports
-------
  inputs.py   : user configuration flags and paths
  geometry.py : pydagmc model, materials module, chunk helpers

Does NOT import neutronics_model.py.

Load
-------
  statepoint.h5         : Tally scoring
  run_meta.json         : Tally IDs and neutron_source_rate
  structural_maps.json  : all structural data

Behaviour differences between sim types
----------------------------------------
tokamak:
  - All keys in cfg.NEUTRONICS_KEYS_TO_PROCESS are processed (OB + IB)
  - Albedo calculation is enabled (cfg.DO_ALBEDO)
  - Poloidal layer sweeps (process_layer_region1_only) are run
  - Neutron source rate default for (1/16) tokamak

slab:
  - Only the single slab chunk key is processed (cfg.ALBEDO_CHUNK_KEY)
  - Albedo is disabled
  - Poloidal layer sweeps are disabled
  - Neutron source rate is scaled by neutron_ratio_source (from run_meta.json)
    x J_in_mean from armor_current_neutron.json
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import openmc

import inputs as cfg
import geometry as geo

# ──────────────────────────────────────────────────────────────────────────────
# Directories
# ──────────────────────────────────────────────────────────────────────────────
RESULTS_DIR = cfg.NEUTRONICS_RESULTS_DIR
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ──────────────────────────────────────────────────────────────────────────────
# Statepoint detection
# ──────────────────────────────────────────────────────────────────────────────
def _find_latest_statepoint(run_dir: Path) -> Path:
    candidates = sorted(run_dir.glob("statepoint.*.h5"), key=lambda p: p.stat().st_mtime)
    if not candidates:
        raise FileNotFoundError(f"No statepoint files found in: {run_dir}")
    latest = candidates[-1]
    print(f"  Last modified: {__import__('datetime').datetime.fromtimestamp(latest.stat().st_mtime)}")
    return latest


STATEPOINT_FILE = _find_latest_statepoint(cfg.NEUTRONICS_RUN_DIR)
print(f"Using statepoint: {STATEPOINT_FILE.name}")

# ──────────────────────────────────────────────────────────────────────────────
# Load Reference Run (gold)
# ──────────────────────────────────────────────────────────────────────────────
from reference_manager import load_reference

_ref_data = None
if cfg.REFERENCE_LABEL:
    try:
        _ref_data = load_reference(cfg.REFERENCE_LABEL)
        print(f"[reference] Loaded gold run: {cfg.REFERENCE_LABEL}")
    except FileNotFoundError:
        print(f"[reference] WARNING: reference '{cfg.REFERENCE_LABEL}' not found — skipping overlay")

# ──────────────────────────────────────────────────────────────────────────────
# Structural maps — loaded from JSON written by neutronics_model.py
# ──────────────────────────────────────────────────────────────────────────────
def _load_structural_maps_json(json_path: Path) -> dict:
    '''
    Load structural map information from .json file produced in neutronics_model.py
    (Can be removed once model/post run together)
    '''
    if not json_path.is_file():
        raise FileNotFoundError(
            f"Structural maps JSON not found: {json_path}\n"
            f"Run neutronics_model.py first to generate it."
        )
    with open(json_path, encoding="utf-8") as f:
        raw = json.load(f)

    def _int_keys(d: dict) -> dict:
        return {int(k): v for k, v in d.items()}

    return {
        "breeder_type":              raw["breeder_type"],
        "structural_nuclides":       raw["structural_nuclides"],
        "cell_volumes":              _int_keys(raw["cell_volumes"]),
        "cell_nuclide_atoms":        _int_keys(raw["cell_nuclide_atoms"]),
        "cell_struct_nuclide_atoms": _int_keys(raw["cell_struct_nuclide_atoms"]),
        "cell_total_atoms_struct":   _int_keys(raw["cell_total_atoms_struct"]),
        "cell_struct_origin_frac":   _int_keys(raw["cell_struct_origin_frac"]),
    }


STRUCTURAL_MAPS_JSON = cfg.NEUTRONICS_RESULTS_DIR / "structural_maps.json"
_struct_maps         = _load_structural_maps_json(STRUCTURAL_MAPS_JSON)

structural_nuclides       = _struct_maps["structural_nuclides"]
cell_volumes              = _struct_maps["cell_volumes"]
cell_nuclide_atoms        = _struct_maps["cell_nuclide_atoms"]
cell_struct_nuclide_atoms = _struct_maps["cell_struct_nuclide_atoms"]
cell_total_atoms_struct   = _struct_maps["cell_total_atoms_struct"]
cell_struct_origin_frac   = _struct_maps["cell_struct_origin_frac"]

# ──────────────────────────────────────────────────────────────────────────────
# Run metadata — loaded from run_meta.json written by neutronics_model.py
# ──────────────────────────────────────────────────────────────────────────────
_run_meta_path = cfg.NEUTRONICS_RESULTS_DIR / "run_meta.json"
if not _run_meta_path.is_file():
    raise FileNotFoundError(
        f"run_meta.json not found: {_run_meta_path}\n"
        f"Run neutronics_model.py first to generate it."
    )
with open(_run_meta_path, encoding="utf-8") as f:
    _run_meta = json.load(f)

neutron_ratio_source: float = float(_run_meta["neutron_ratio_source"])
_tally_ids: dict             = _run_meta.get("tally_ids", {})

if not _tally_ids:
    raise RuntimeError(
        "run_meta.json has no 'tally_ids' section.\n"
        "Re-run neutronics_model.py to regenerate run_meta.json with tally IDs."
    )

# ──────────────────────────────────────────────────────────────────────────────
# Tallies — loaded by ID from run_meta.json
# ──────────────────────────────────────────────────────────────────────────────
_sp_for_tallies = openmc.StatePoint(str(STATEPOINT_FILE))
_all_tallies    = _sp_for_tallies.tallies   # dict: id -> Tally

# Use this function for single tallies (flux, heating, ngamma, etc..)
def _get_tally_by_id(tid, desc: str) -> openmc.Tally:
    if tid is None:
        raise RuntimeError(
            f"Tally ID for '{desc}' is None in run_meta.json. "
            f"Re-run neutronics_model.py."
        )
    tid = int(tid)
    if tid not in _all_tallies:
        raise RuntimeError(
            f"Tally '{desc}' (id={tid}) not found in statepoint. "
            f"Available IDs: {sorted(_all_tallies.keys())}"
        )
    return _all_tallies[tid]

flux_tally    = _get_tally_by_id(_tally_ids["flux_spectrum"], "flux_spectrum")
heating_tally = _get_tally_by_id(_tally_ids["heating"],       "heating")

# total_current is tokamak-only — None for slab
t_current_tally: Optional[openmc.Tally] = None
p_current_tallies: Dict[int, openmc.Tally] = {}

if cfg.SIM_TYPE == "tokamak":
    t_current_tally = (
        _all_tallies.get(int(_tally_ids["total_current"]))
        if _tally_ids.get("total_current") is not None
        else None
    )
    p_current_tallies = {
        int(cid): _all_tallies[int(tid)]
        for cid, tid in _tally_ids.get("partial_current", {}).items()
        if int(tid) in _all_tallies
    }

# DPA/gas tally IDs 
_dpa_gas_ids: Dict[int, int] = {
    int(cid): int(tid)
    for cid, tid in _tally_ids.get("dpa_gas", {}).items()
}

_sp_for_tallies.close()

# ──────────────────────────────────────────────────────────────────────────────
# Source scaling constants
# ──────────────────────────────────────────────────────────────────────────────
ev_to_joule: float = geo.ev_to_joule
s_in_y:      float = geo.s_in_y

if cfg.SIM_TYPE in ("slab", "fast_slab"):
    if not cfg.ARMOR_CURRENT_JSON.is_file():
        raise FileNotFoundError(
            f"Armor current JSON not found: {cfg.ARMOR_CURRENT_JSON}\n"
            f"Run tokamak neutronics_post.py first to generate it."
        )
    with open(cfg.ARMOR_CURRENT_JSON) as f:
        _records = json.load(f)
    _J_in = float(_records[0]["J_in_mean"])
    neutron_source_rate: float = (
        neutron_ratio_source * _J_in
        * (cfg.TOTAL_FUSION_POWER_W / cfg.NUMBER_OF_SECTORS)
        / (ev_to_joule * cfg.EV_PER_FUSION)
    )
    print(f"[slab] neutron_source_rate (surface-scaled) = {neutron_source_rate:.4e}")
else:
    neutron_source_rate: float = geo.neutron_source_rate

# ──────────────────────────────────────────────────────────────────────────────
# Energy / lethargy bins
# ──────────────────────────────────────────────────────────────────────────────
energies      = openmc.mgxs.GROUP_STRUCTURES["CCFE-709"]
unit_lethargy = np.array(
    [np.log(energies[i + 1] / energies[i]) for i in range(len(energies) - 1)],
    dtype=float,
)
E_mid = 0.5 * (np.asarray(energies[:-1], float) + np.asarray(energies[1:], float))

# ──────────────────────────────────────────────────────────────────────────────
# Chunk helpers — from geometry
# ──────────────────────────────────────────────────────────────────────────────
ob_by_key           = geo.ob_by_key
ib_by_key           = geo.ib_by_key
n_breeder           = geo.n_breeder
OB_CHUNK_SIZE       = geo.OB_CHUNK_SIZE
radial_bins_for_key = geo._chunk_result["radial_bins_for_key"]

# ──────────────────────────────────────────────────────────────────────────────
# Colour helpers
# ──────────────────────────────────────────────────────────────────────────────
def _generate_colors(n: int):
    cmap = plt.get_cmap("rainbow")
    return cmap(np.linspace(1, 0, max(n, 1)))

colors = _generate_colors(OB_CHUNK_SIZE)

# ──────────────────────────────────────────────────────────────────────────────
# Tally / array helpers
# ──────────────────────────────────────────────────────────────────────────────
GAS_SCORES_EXPLICIT = {"H1-production", "He3-production", "He4-production"}
GAS_SCORES_REACTION = {"(n,Xp)", "(n,X3He)", "(n,Xa)"}

def build_dpa_gas_map(sp: openmc.StatePoint) -> Dict[int, openmc.Tally]:
    """
    Load dpa/gas tallies directly by ID from run_meta.json.
    Falls back to score-based scanning if run_meta has no dpa_gas entries
    (e.g. running against an old statepoint).
    """
    if _dpa_gas_ids:
        out: Dict[int, openmc.Tally] = {}
        for cid, tid in _dpa_gas_ids.items():
            if tid in sp.tallies:
                out[cid] = sp.tallies[tid]
            else:
                print(f"[warn] dpa_gas tally id={tid} for cell {cid} "
                      f"not found in statepoint, skipping.")
        return out

def subset_by_cells(
    data: np.ndarray, cell_bins: List[int], desired_cells: List[int]
    ) -> np.ndarray:
    idx = [cell_bins.index(int(cid)) for cid in desired_cells]
    return data[idx]

def get_cell_bins(tally: openmc.Tally) -> List[int]:
    cf = next(f for f in tally.filters if isinstance(f, openmc.CellFilter))
    return [int(x) for x in cf.bins]

def get_flux_spectrum_arrays(
    t: openmc.Tally,
    ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Return flux_mean, flux_std with shape:
        (n_cells, n_particles, n_energy)
    """
    flux_mean = np.asarray(t.get_reshaped_data(value="mean"),    dtype=float)
    flux_std  = np.asarray(t.get_reshaped_data(value="std_dev"), dtype=float)

    while flux_mean.ndim > 3 and flux_mean.shape[-1] == 1:
        flux_mean = flux_mean[..., 0]
    while flux_std.ndim > 3 and flux_std.shape[-1] == 1:
        flux_std = flux_std[..., 0]
    if flux_mean.ndim != 3:
        raise ValueError(f"Unexpected flux shape: {flux_mean.shape}")
    return flux_mean, flux_std

def scaling_for_cells(cell_ids: List[int]) -> Dict[int, float]:
    """
    Returns neutron_source_rate / cell_volume for each cell.
    Volumes are read from structural_maps.json — no live OpenMC objects needed.
    """
    out: Dict[int, float] = {}
    for cid in cell_ids:
        cid = int(cid)
        vol = cell_volumes.get(cid, 0.0)
        if not vol or vol <= 0.0:
            raise ValueError(
                f"Cell {cid} has no valid volume in structural_maps.json ({vol}). "
                f"Re-run neutronics_model.py to regenerate."
            )
        out[cid] = neutron_source_rate / float(vol)
    return out

def require_tally_id(sp: openmc.StatePoint, tally: openmc.Tally, desc: str) -> int:
    tid = int(tally.id)
    if tid not in sp.tallies:
        raise RuntimeError(f"{desc}: tally id={tid} not in statepoint")
    return tid

def tally_mean_std(
    t: openmc.Tally, *, score: str, nuclide: str
    ) -> Tuple[float, float]:
    m = float(t.get_values(scores=[score], nuclides=[nuclide],
                            value="mean").ravel()[0])
    s = float(t.get_values(scores=[score], nuclides=[nuclide],
                            value="std_dev").ravel()[0])
    return m, s

def corrected_sum_detailed(
    t: openmc.Tally,
    *,
    score: str,
    nuclides: List[str],
    f_struct: Dict[str, float],
    ) -> Tuple[float, float, Dict[str, Tuple[float, float]]]:
    """
    Returns (mean_total, std_total, per_nuc) where
    per_nuc[nuc] = (weighted_mean, weighted_std)
    """
    mean_tot = 0.0
    var_tot  = 0.0
    per_nuc: Dict[str, Tuple[float, float]] = {}

    for nuc in nuclides:
        nuc = str(nuc)
        w   = float(f_struct.get(nuc, 0.0))
        if w == 0.0:
            continue
        try:
            m, sd = tally_mean_std(t, score=score, nuclide=nuc)
        except Exception:
            continue
        wm = w * m
        ws = w * sd
        mean_tot        += wm
        var_tot         += ws ** 2
        per_nuc[nuc]     = (wm, ws)

    return mean_tot, math.sqrt(max(var_tot, 0.0)), per_nuc

def element_from_nuclide(nuc: str) -> str:
    m = re.match(r"[A-Za-z]+", nuc)
    return m.group(0) if m else nuc

# ──────────────────────────────────────────────────────────────────────────────
# CSV helpers
# ──────────────────────────────────────────────────────────────────────────────
def make_chunk_base_df(
    *,
    chunk_key: str,
    cell_ids: List[int],
    labels: List[str],
    xcent: np.ndarray,
    xedges: np.ndarray,
    ) -> pd.DataFrame:
    n = len(cell_ids)
    return pd.DataFrame({
        "chunk_key":   [chunk_key] * n,
        "cell_id":     [int(c) for c in cell_ids],
        "layer_label": list(labels),
        "x_center_cm": np.asarray(xcent,  float).ravel(),
        "x_left_cm":   np.asarray(xedges, float).ravel()[:-1],
        "x_right_cm":  np.asarray(xedges, float).ravel()[1:],
    })

def make_layer_base_df(
    *,
    layer_tag: str,
    cell_ids: List[int],
    chunk_labels: List[str],
    poloidal_index: List[int],
    ) -> pd.DataFrame:
    """
    Base DataFrame for poloidal-layer sweeps (process_layer_region1_only).

    No radial x_cm columns — one row per poloidal region.
    ``chunk_key`` holds the poloidal chunk label (e.g. 'OB_1_b2'),
    ``layer_tag`` identifies which radial layer this file represents
    (e.g. 'armor', 'fw', 'vv').
    ``poloidal_index`` is the integer position along the poloidal sweep
    (0-based), used for stable ordering.
    """
    n = len(cell_ids)
    return pd.DataFrame({
        "layer_tag":      [layer_tag] * n,
        "chunk_key":      list(chunk_labels),
        "cell_id":        [int(c) for c in cell_ids],
        "poloidal_index": list(poloidal_index),
    })

def save_profile(
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

    n    = len(base_df)
    mean = np.asarray(mean, float).ravel()
    std  = np.asarray(std,  float).ravel()
    if len(mean) != n or len(std) != n:
        raise ValueError(
            f"{quantity}: mean/std length mismatch with base_df "
            f"({len(mean)}, {len(std)}) vs {n}"
        )

    df             = base_df.copy()
    df["quantity"] = quantity
    df["units"]    = units
    df["mean"]     = mean
    df["std"]      = std

    lower = mean - std
    if nonnegative_lower:
        lower = np.maximum(lower, 0.0)

    df["lower_1sigma"] = lower
    df["upper_1sigma"] = mean + std

    if "x_left_cm" in df.columns:
        df = df.sort_values(["x_left_cm", "cell_id"]).reset_index(drop=True)
    elif "poloidal_index" in df.columns:
        df = df.sort_values("poloidal_index").reset_index(drop=True)

    # Poloidal-layer base_df uses layer_tag as the file discriminator;
    # radial-chunk base_df uses chunk_key.
    if "layer_tag" in df.columns:
        file_key = str(df["layer_tag"].iloc[0])
    else:
        file_key = str(df["chunk_key"].iloc[0])

    path = outdir / f"profile_{quantity}_{file_key}.csv"
    df.to_csv(path, index=False)
    return path

def write_struct_origin_csv(
    outdir: Path,
    chunk_key: str,
    cell_ids: List[int],
    labels: List[str],
    cell_struct_origin_frac: Dict[int, Dict[str, float]],
    ) -> None:
    rows = []
    for cid, lab in zip(cell_ids, labels):
        cid      = int(cid)
        frac_map = cell_struct_origin_frac.get(cid, {}) or {}
        for nuc, f in sorted(frac_map.items()):
            rows.append(dict(
                chunk_key=chunk_key,
                cell_id=cid,
                layer_label=lab,
                nuclide=str(nuc),
                f_struct_origin=float(f),
            ))
    pd.DataFrame(rows).to_csv(
        outdir / f"struct_origin_fractions_{chunk_key}.csv", index=False
    )

# ──────────────────────────────────────────────────────────────────────────────
# y-axis helper
# ──────────────────────────────────────────────────────────────────────────────
def set_ylim_and_ticks(
    ax,
    ydata,
    *,
    ratio: float = 0.3,
    scale: str = "auto",
    eps: float = 1e-99,
    n_linear_ticks: int = 6,
    log_threshold_decades: float = 1.0,
    min_log_pad_decades: float = 0.10,
    ) -> None:
    ratio = max(0.0, min(float(ratio), 0.95))

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
            lo   = float(logy.min())
            hi   = float(logy.max())
            span = max(hi - lo, 1e-12)
            pad  = max(ratio * span, float(min_log_pad_decades))
            ax.set_ylim(max(10 ** (lo - pad), eps), 10 ** (hi + pad))
            ax.yaxis.set_major_locator(
                mticker.LogLocator(base=10.0, subs=(1.0, 2.0, 5.0))
            )
            ax.yaxis.set_major_formatter(
                mticker.LogFormatterMathtext(base=10.0, labelOnlyBase=False)
            )
            ax.yaxis.set_minor_locator(
                mticker.LogLocator(base=10.0, subs=np.arange(1, 10) * 0.1)
            )
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
        span  = ymax - ymin
        ymin -= ratio * span
        ymax += ratio * span

    ax.set_ylim(ymin, ymax)
    ax.yaxis.set_major_locator(mticker.MaxNLocator(n_linear_ticks))
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("{x:g}"))

# ──────────────────────────────────────────────────────────────────────────────
# Overlay Gold Run
# ──────────────────────────────────────────────────────────────────────────────

def _overlay_ref(ax: plt.Axes, chunk_key: str, quantity: str,
                 label: str = "W Armor + Eurofer Structural",
                 color: str = "black",
                 ls: str = "--") -> None:
    if _ref_data is None or _ref_data.get("neutronics_dir") is None:
        return
    csv = (
        _ref_data["neutronics_dir"]
        / chunk_key / "profiles" / f"profile_{quantity}_{chunk_key}.csv"
    )
    if not csv.is_file():
        return
    df = pd.read_csv(csv)
    if "x_left_cm" not in df.columns or "mean" not in df.columns:
        return
    df = df.sort_values("x_left_cm").reset_index(drop=True)
    x  = df["x_left_cm"].to_numpy(float)
    y  = df["mean"].to_numpy(float)
    lo = df["lower_1sigma"].to_numpy(float)
    hi = df["upper_1sigma"].to_numpy(float)

    if "x_right_cm" in df.columns:
        x_final = float(df["x_right_cm"].iloc[-1])
    else:
        x_final = x[-1] + (x[-1] - x[-2])

    xp = np.append(x, x_final)
    ax.step(xp, np.append(y, y[-1]),
            where="post", color=color, ls=ls, lw=1.2,
            alpha=0.65, label=label, zorder=0)
    ax.fill_between(xp, np.append(lo, lo[-1]), np.append(hi, hi[-1]),
                    step="post", color=color, alpha=0.08, zorder=0)

def _overlay_ref_poloidal(ax: plt.Axes, layer_tag: str, quantity: str,
                           label: str = "W Armor + Eurofer Structural",
                           color: str = "black",
                           ls: str = "--") -> None:
    if _ref_data is None or _ref_data.get("neutronics_dir") is None:
        return
    csv = (
        _ref_data["neutronics_dir"]
        / f"layer_{layer_tag}_region1_only"
        / "profiles"
        / f"profile_{quantity}_{layer_tag}.csv"
    )
    if not csv.is_file():
        return
    try:
        df = pd.read_csv(csv)
        if "poloidal_index" not in df.columns or "mean" not in df.columns:
            return
        df = df.sort_values("poloidal_index").reset_index(drop=True)
        y  = df["mean"].to_numpy(float)
        lo = df["lower_1sigma"].to_numpy(float)
        hi = df["upper_1sigma"].to_numpy(float)
        yp  = _pad(y)
        lop = _pad(lo)
        hip = _pad(hi)
        xp = _pad(np.arange(len(y), dtype=float))
        ax.step(xp, yp, where="mid", color=color, ls=ls, lw=1.2,
                alpha=0.65, label=label, zorder=0)
        ax.fill_between(xp, lop, hip, step="mid",
                        color=color, alpha=0.08, zorder=0)
    except Exception as e:
        print(f"[reference] poloidal overlay failed for {quantity}/{layer_tag}: {e}")

# ──────────────────────────────────────────────────────────────────────────────
# Layer naming
# ──────────────────────────────────────────────────────────────────────────────
def layer_names_for_chunk(chunk_cells: List[int], region_tag: str) -> List[str]:
    n        = len(chunk_cells)
    n_layers = n - 3
    labels   = ["Armor", "First_Wall"]
    labels  += [f"Breeder layer {i+1}" for i in range(n_layers)]
    labels  += ["Vacuum Vessel"]
    return labels

# ──────────────────────────────────────────────────────────────────────────────
# ALBEDO  (tokamak only)
# ──────────────────────────────────────────────────────────────────────────────
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
        vols   = [v for v in senses if v is not None]
        return len(vols) == 1

def _cell_surface_ids_from_info_flat(info_flat: dict, cell_id: int) -> List[int]:
    out = []
    for rec in info_flat.get(int(cell_id), []):
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

def _build_special_surfaces(
    pydagmc_model,
    info_flat: dict,
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

    special:     Dict[str, int]    = {}
    owner:       Dict[str, int]    = {}
    special_ids: Set[int]          = set()
    debug:       Dict[str, object] = {}

    armor_sids   = [sid for sid in _cell_surface_ids_from_info_flat(info_flat, armor_cell_id) if sid in allowed]
    armor_ranked = _rank_sids_by_area(pydagmc_model, armor_sids)

    armor_front             = next((sid for sid in armor_ranked if sid in ext_set), None)
    armor_internal_ranked   = [sid for sid in armor_ranked if sid in int_set]
    armor_back              = armor_internal_ranked[0] if armor_internal_ranked else None

    debug["armor"] = {
        "cell_id":           int(armor_cell_id),
        "ranked_top":        armor_ranked[:8],
        "picked_front_ext":  armor_front,
        "picked_back_int":   armor_back,
        "internal_candidates": armor_internal_ranked[:8],
    }

    if armor_front is not None:
        special["Armor_front_ext"] = int(armor_front)
        owner["Armor_front_ext"]   = int(armor_cell_id)
        special_ids.add(int(armor_front))

    if armor_back is not None:
        special["Armor_back_int"] = int(armor_back)
        owner["Armor_back_int"]   = int(armor_cell_id)
        special_ids.add(int(armor_back))

    breeder_sids   = [sid for sid in _cell_surface_ids_from_info_flat(info_flat, breeder_back_cell_id) if sid in allowed]
    breeder_ranked = _rank_sids_by_area(pydagmc_model, breeder_sids)
    breeder_ext    = next((sid for sid in breeder_ranked if sid in ext_set), None)

    debug["breeder_back"] = {
        "cell_id":           int(breeder_back_cell_id),
        "ranked_top":        breeder_ranked[:8],
        "picked_largest_ext": breeder_ext,
    }

    if breeder_ext is not None:
        special["Breeder_back_ext"] = int(breeder_ext)
        owner["Breeder_back_ext"]   = int(breeder_back_cell_id)
        special_ids.add(int(breeder_ext))

    vv_sids   = [sid for sid in _cell_surface_ids_from_info_flat(info_flat, vv_cell_id) if sid in allowed]
    vv_ranked = _rank_sids_by_area(pydagmc_model, vv_sids)

    vv_picks   = []
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
        "cell_id":          int(vv_cell_id),
        "ranked_top":       vv_ranked[:12],
        "picked_unshared":  vv_picks,
        "first_rejects":    vv_rejects[:8],
    }

    if len(vv_picks) >= 2:
        special["VV_face_1"] = int(vv_picks[0])
        special["VV_face_2"] = int(vv_picks[1])
        owner["VV_face_1"]   = int(vv_cell_id)
        owner["VV_face_2"]   = int(vv_cell_id)
        special_ids.add(int(vv_picks[0]))
        special_ids.add(int(vv_picks[1]))

    return special, special_ids, owner, debug

def compute_albedo_for_chunk(
    sp: openmc.StatePoint,
    *,
    cell_ids: List[int],
    labels: List[str],
    outdir: Path,
    chunk_key: str,
    excluded_surfaces: Optional[Set[int]] = None,
) -> Dict[str, object]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    EXCLUDED_SURFACES = set(excluded_surfaces or set())

    if len(labels) != len(cell_ids):
        raise ValueError(
            f"{chunk_key}: labels length {len(labels)} != cell_ids length {len(cell_ids)}"
        )

    if t_current_tally is None or not p_current_tallies:
        print(f"[warn] {chunk_key}: skipping albedo – current tallies not available")
        return {}

    info        = geo.info
    ext_ids     = [int(x) for x in geo.external_surface_ids]
    int_ids     = [int(x) for x in geo.internal_surface_ids]
    all_ids     = [int(x) for x in getattr(geo, "all_surface_ids", [])] or sorted(set(ext_ids) | set(int_ids))
    ext_set     = set(ext_ids)
    int_set     = set(int_ids)
    all_set     = set(all_ids)

    PARTICLES      = ("neutron", "photon")
    eps            = 1e-15
    chunk_cell_set = set(int(c) for c in cell_ids)
    idx_map        = {int(cid): i for i, cid in enumerate(cell_ids)}

    def _cell_surfaces(cid: int, which: str) -> List[dict]:
        x = info.get(int(cid), {})
        if isinstance(x, dict):
            return list(x.get(which, []) or [])
        return list(x or [])

    def _read_current_df_map(tally) -> Dict[Tuple[int, str], Tuple[float, float]]:
        df = tally.get_pandas_dataframe()
        out = {}
        for _, r in df.iterrows():
            sid  = int(r["surface"])
            part = str(r["particle"]).lower()
            out[(sid, part)] = (float(r["mean"]), float(r["std. dev."]))
        return out

    surf_to_cells: Dict[int, Set[int]] = {}
    for cid in cell_ids:
        cid = int(cid)
        for s in _cell_surfaces(cid, "all_surfaces"):
            try:
                sid = int(s["surface_id"])
            except Exception:
                continue
            surf_to_cells.setdefault(sid, set()).add(cid)

    if not surf_to_cells:
        print(f"[warn] {chunk_key}: no surfaces mapped from geo.info; skipping albedo.")
        return {}

    Jnet_by_sid_particle: Dict[Tuple[int, str], Tuple[float, float]] = {}
    try:
        t_tot                = sp.get_tally(id=t_current_tally.id)
        Jnet_by_sid_particle = _read_current_df_map(t_tot)
    except Exception as e:
        print(f"[warn] {chunk_key}: could not read total current tally; "
              f"external albedo limited. ({e})")

    allowed_surface_ids = set(all_set)

    special:             Dict[str, int]    = {}
    special_surface_ids: Set[int]          = set()
    special_owner_cell:  Dict[str, int]    = {}
    special_debug:       Dict[str, object] = {}

    if len(cell_ids) >= 3:
        bm_info_flat = {int(cid): _cell_surfaces(int(cid), "all_surfaces") for cid in cell_ids}
        special, special_surface_ids, special_owner_cell, special_debug = _build_special_surfaces(
            geo.pydagmc_model,
            bm_info_flat,
            armor_cell_id=int(cell_ids[0]),
            breeder_back_cell_id=int(cell_ids[-2]),
            vv_cell_id=int(cell_ids[-1]),
            allowed_surface_ids=allowed_surface_ids,
            external_surface_ids=ext_set,
            internal_surface_ids=int_set,
        )
        print(f"[debug] {chunk_key} specials picked: {special}")
    else:
        print(f"[warn] {chunk_key}: chunk has only {len(cell_ids)} cells; "
              f"skipping special selection.")

    partial_by_cell_sid_particle: Dict[Tuple[int, int, str], Tuple[float, float]] = {}
    for cid, tt in p_current_tallies.items():
        cid = int(cid)
        if cid not in chunk_cell_set:
            continue
        try:
            tally = sp.get_tally(id=tt.id)
            d     = _read_current_df_map(tally)
            for (sid, part), (m, s) in d.items():
                partial_by_cell_sid_particle[(cid, sid, part)] = (m, s)
        except Exception as e:
            print(f"[warn] {chunk_key}: failed reading partial tally for cell {cid}: {e}")

    SPECIAL_LABELS = {
        "Armor_front_ext":  "Armor Front",
        "Armor_back_int":   "Armor Back",
        "Breeder_back_ext": "Breeder Back",
        "VV_face_1":        "VV Back",
        "VV_face_2":        "VV Front",
    }
    SPECIAL_COLORS = {
        "Armor_front_ext":  "tab:red",
        "Armor_back_int":   "tab:orange",
        "Breeder_back_ext": "tab:purple",
        "VV_face_1":        "tab:green",
        "VV_face_2":        "tab:olive",
    }
    SPECIAL_MARKERS  = {k: "x" for k in SPECIAL_LABELS}
    MIN_VISIBLE_YERR = 1e-4

    results_by_particle: Dict[str, Dict[str, object]] = {}

    for particle in PARTICLES:
        rows: List[dict] = []

        # External surfaces
        for sid in ext_ids:
            cset = surf_to_cells.get(int(sid), set()) & chunk_cell_set
            if not cset:
                continue
            cid = sorted(cset)[0]
            if (sid, particle) not in Jnet_by_sid_particle:
                continue

            Jnet_m, Jnet_s = Jnet_by_sid_particle[(sid, particle)]
            Jout_m, Jout_s = partial_by_cell_sid_particle.get((cid, sid, particle), (0.0, 0.0))
            Jnet_a = abs(Jnet_m)
            Jout_a = abs(Jout_m)
            Jin    = abs(Jout_a - Jnet_a)

            if Jout_a < eps:
                A_mean = A_std = 0.0
            else:
                A_mean   = Jin / Jout_a
                dA_dJout = Jnet_a / (Jout_a ** 2)
                dA_dJnet = -1.0   / Jout_a
                A_var    = (dA_dJout ** 2) * (Jout_s ** 2) + (dA_dJnet ** 2) * (Jnet_s ** 2)
                A_std    = math.sqrt(max(A_var, 0.0))

            rows.append({
                "particle":     particle,
                "surface_id":   sid,
                "cell_id":      cid,
                "kind":         "external",
                "J_total_mean": float(Jnet_m),
                "J_total_std":  float(Jnet_s),
                "J_out_mean":   float(Jout_m),
                "J_out_std":    float(Jout_s),
                "albedo_mean":  float(A_mean),
                "albedo_std":   float(A_std),
            })

        # Internal surfaces
        for sid in int_ids:
            cset = surf_to_cells.get(int(sid), set()) & chunk_cell_set
            if len(cset) != 2:
                continue
            c1, c2   = sorted(cset)
            Jout1_m, Jout1_s = partial_by_cell_sid_particle.get((c1, sid, particle), (0.0, 0.0))
            Jout2_m, Jout2_s = partial_by_cell_sid_particle.get((c2, sid, particle), (0.0, 0.0))
            J1, J2   = abs(Jout1_m), abs(Jout2_m)
            s1, s2   = abs(Jout1_s), abs(Jout2_s)

            if J1 <= eps or J2 <= eps:
                A1 = A2 = A1_std = A2_std = float("nan")
            else:
                A1      = J2 / J1
                dA1_dJ2 = 1.0 / J1
                dA1_dJ1 = -J2 / (J1 ** 2)
                A1_var  = (dA1_dJ2 ** 2) * (s2 ** 2) + (dA1_dJ1 ** 2) * (s1 ** 2)
                A1_std  = math.sqrt(max(A1_var, 0.0))

                A2      = J1 / J2
                dA2_dJ1 = 1.0 / J2
                dA2_dJ2 = -J1 / (J2 ** 2)
                A2_var  = (dA2_dJ1 ** 2) * (s1 ** 2) + (dA2_dJ2 ** 2) * (s2 ** 2)
                A2_std  = math.sqrt(max(A2_var, 0.0))

            rows.append({
                "particle":     particle,
                "surface_id":   sid,
                "cell_id":      c1,
                "kind":         "internal",
                "J_total_mean": np.nan,
                "J_total_std":  np.nan,
                "J_out_mean":   float(Jout1_m),
                "J_out_std":    float(Jout1_s),
                "albedo_mean":  float(A1)    if np.isfinite(A1)    else np.nan,
                "albedo_std":   float(A1_std) if np.isfinite(A1_std) else np.nan,
            })
            rows.append({
                "particle":     particle,
                "surface_id":   sid,
                "cell_id":      c2,
                "kind":         "internal",
                "J_total_mean": np.nan,
                "J_total_std":  np.nan,
                "J_out_mean":   float(Jout2_m),
                "J_out_std":    float(Jout2_s),
                "albedo_mean":  float(A2)    if np.isfinite(A2)    else np.nan,
                "albedo_std":   float(A2_std) if np.isfinite(A2_std) else np.nan,
            })

        if not rows:
            print(f"[warn] {chunk_key}: no albedo rows created for {particle}; skipping outputs.")
            continue

        df_alb = pd.DataFrame(rows).sort_values(["kind", "surface_id", "cell_id"])
        df_alb.to_csv(
            outdir / f"surface_albedo_all_{particle}_{chunk_key}.csv", index=False
        )

        # Save special surface rows
        # (Armor_front/back, VV_front/back, Last_breeder_back)
        special_rows: List[dict] = []
        for key, sid in special.items():
            want_cid = int(special_owner_cell.get(key, -1))
            cand     = df_alb[df_alb["surface_id"] == int(sid)]
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
        df_special.to_csv(
            outdir / f"surface_albedo_special_{particle}_{chunk_key}.csv", index=False
        )

        # Armor current JSON (save armor_current.json file for slab source)
        incident_records: List[dict] = []
        for _, r in df_special.iterrows():
            key = str(r.get("special_key", ""))
            if key != "Armor_front_ext":
                continue
            J_out     = abs(float(r["J_out_mean"]))
            J_net     = abs(float(r["J_total_mean"])) if pd.notna(r["J_total_mean"]) else 0.0
            J_in      = abs(J_out - J_net)
            J_out_std = abs(float(r["J_out_std"]))
            J_net_std = abs(float(r["J_total_std"])) if pd.notna(r["J_total_std"]) else 0.0
            J_in_std  = math.sqrt(J_out_std ** 2 + J_net_std ** 2)
            incident_records.append({
                "particle":    particle,
                "special_key": key,
                "surface_id":  int(r["surface_id"]),
                "cell_id":     int(r["cell_id"]),
                "J_out_mean":  J_out,
                "J_out_std":   J_out_std,
                "J_net_mean":  float(r["J_total_mean"]) if pd.notna(r["J_total_mean"]) else None,
                "J_net_std":   J_net_std,
                "J_in_mean":   J_in,
                "J_in_std":    J_in_std,
                "chunk_key":   chunk_key,
            })
        if incident_records:
            out_path = outdir / f"armor_current_{particle}.json"
            with open(out_path, "w") as f:
                json.dump(incident_records, f, indent=2)
            print(f"[completed] {out_path}")

        # Save external filtered
        df_external_rest = df_alb[
            (df_alb["kind"] == "external") &
            (~df_alb["surface_id"].isin(EXCLUDED_SURFACES)) &
            (~df_alb["surface_id"].isin(special_surface_ids))
        ].copy()
        df_external_rest.to_csv(
            outdir / f"surface_albedo_external_filtered_{particle}_{chunk_key}.csv",
            index=False,
        )

        # Save cell summary 
        cell_rows: List[dict] = []
        for cid, grp in df_external_rest.groupby("cell_id"):
            cid    = int(cid)
            A_mean = float(grp["albedo_mean"].mean()) if len(grp) else float("nan")
            st     = grp["albedo_std"].to_numpy(dtype=float)
            A_std  = float(np.sqrt(np.nansum(st * st)) / max(len(st), 1)) if len(st) else float("nan")
            cell_rows.append({
                "particle":    particle,
                "cell_id":     cid,
                "layer_label": labels[idx_map[cid]] if cid in idx_map else f"cell_{cid}",
                "albedo_mean": A_mean,
                "albedo_std":  A_std,
                "n_surfaces":  int(len(grp)),
            })

        df_cell = (
            pd.DataFrame(cell_rows).sort_values("cell_id")
            if cell_rows
            else pd.DataFrame(columns=[
                "particle", "cell_id", "layer_label",
                "albedo_mean", "albedo_std", "n_surfaces",
            ])
        )
        df_cell.to_csv(
            outdir / f"cell_albedo_summary_external_{particle}_{chunk_key}.csv",
            index=False,
        )

        # Plot
        xpos = np.arange(len(labels))
        y    = np.full(len(labels), np.nan, dtype=float)
        e    = np.full(len(labels), np.nan, dtype=float)
        for _, row in df_cell.iterrows():
            cid = int(row["cell_id"])
            if cid not in idx_map:
                continue
            i    = idx_map[cid]
            y[i] = float(row["albedo_mean"]) if pd.notna(row["albedo_mean"]) else np.nan
            e[i] = float(row["albedo_std"])  if pd.notna(row["albedo_std"])  else np.nan

        plt.figure(figsize=(11, 6))
        plt.errorbar(xpos, y, yerr=e, fmt="none", elinewidth=1, capsize=4)
        plt.scatter(
            xpos, y, marker="x", s=90, linewidths=2,
            label=f"{particle.capitalize()} cell avg (external, excluding specials)",
        )

        used_labels: Set[str] = set()
        if not df_special.empty:
            for _, r in df_special.iterrows():
                key = str(r.get("special_key", "special"))
                cid = int(r["cell_id"])
                if cid not in idx_map:
                    continue
                x0  = idx_map[cid]
                y0  = float(r["albedo_mean"]) if pd.notna(r["albedo_mean"]) else np.nan
                e0  = float(r["albedo_std"])  if pd.notna(r["albedo_std"])  else np.nan
                if not np.isfinite(y0):
                    continue
                color      = SPECIAL_COLORS.get(key, "black")
                label      = SPECIAL_LABELS.get(key, key)
                marker     = SPECIAL_MARKERS.get(key, "x")
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
                plt.scatter([x0], [y0], marker=marker, s=120,
                            color=color, zorder=6, label=plot_label)

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
            "df_albedo":                   df_alb,
            "df_albedo_external_filtered": df_external_rest,
            "df_special":                  df_special,
            "df_cell_summary":             df_cell,
        }

    if not results_by_particle:
        print(f"[warn] {chunk_key}: no particle-specific albedo outputs were created.")
        return {}

    print(f"[completed] {chunk_key}: albedo done for particles = {list(results_by_particle.keys())}")

    return {
        "by_particle":        results_by_particle,
        "special":            special,
        "special_surface_ids": special_surface_ids,
        "special_owner_cell": special_owner_cell,
        "special_debug":      special_debug,
        "excluded_surfaces":  EXCLUDED_SURFACES,
    }

def save_summary_neutronics_results(
    *,
    chunk_key: str,
    out_csv: "Path",
    results_dir: "Path" = None,
) -> "pd.DataFrame":
    """
    Read the profile CSVs already written by process_chunk() and extract the
    summary (mean ± std) DPA, H, and He appm/fpy for the Armor, First Wall,
    and Vacuum Vessel layers. Save to ``out_csv``.

    Reads (from results_dir / chunk_key /)
    ----------------------------------------
    profile_dpa_fpy_struct_origin_<chunk_key>.csv
    profile_He_appm_fpy_struct_origin_<chunk_key>.csv
    profile_H_appm_fpy_struct_origin_<chunk_key>.csv

    Parameters
    ----------
    chunk_key : str
        The reference chunk, e.g. cfg.TRIGGER_CHUNK_KEY = "OB_1_b6".
    out_csv : Path
        Destination file (cfg.PEAK_RESULTS_CSV).
    results_dir : Path, optional
        Top-level neutronics results directory. Defaults to
        cfg.NEUTRONICS_RESULTS_DIR.
    """
    results_dir = Path(results_dir) if results_dir else cfg.NEUTRONICS_RESULTS_DIR
    chunk_dir = results_dir / chunk_key / "profiles"

    # Label normalisation (just in case)
    LAYER_ALIASES = {
        "Armor":         "Armor",
        "First_Wall":    "First_Wall",
        "First wall":    "First_Wall",
        "Vacuum Vessel": "Vacuum Vessel",
        "VV":            "Vacuum Vessel",
    }

    TARGET_LAYERS = list(cfg.SUMMARY_TARGET_LAYERS)

    # Profile CSV paths
    _profiles = {
        "dpa": f"profile_dpa_fpy_struct_origin_{chunk_key}.csv",
        "He":  f"profile_He_appm_fpy_struct_origin_{chunk_key}.csv",
        "H":   f"profile_H_appm_fpy_struct_origin_{chunk_key}.csv",
    }

    dfs: "Dict[str, pd.DataFrame]" = {}
    for name, fname in _profiles.items():
        fpath = chunk_dir / fname
        if not fpath.is_file():
            raise FileNotFoundError(
                f"Profile CSV not found: {fpath}\n"
                f"Make sure process_chunk() has already run for '{chunk_key}'."
            )
        df = pd.read_csv(fpath)
        df["layer_label"] = df["layer_label"].map(lambda x: LAYER_ALIASES.get(x, x))
        dfs[name] = df

    # Cell lookup helper
    def _get(df: "pd.DataFrame", cid: int, col_mean: str = "mean", col_std: str = "std"):
        """Extract mean and std for cid from a profile dataframe."""
        r = df[df["cell_id"] == cid]
        if r.empty:
            return float("nan"), float("nan")
        return float(r.iloc[0][col_mean]), float(r.iloc[0][col_std])

    # Per-layer extraction 
    rows: "List[dict]" = []

    for layer_label in TARGET_LAYERS:
        ref   = dfs["dpa"]
        match = ref[ref["layer_label"] == layer_label]

        if match.empty:
            print(
                f"[warn] save_summary_neutronics_results: layer '{layer_label}' "
                f"not found in {chunk_key} DPA profile — skipping."
            )
            continue

        if len(match) > 1:
            print(
                f"[warn] save_summary_neutronics_results: {len(match)} rows found "
                f"for layer '{layer_label}' in {chunk_key} — using the first."
            )

        cell_id = int(match.iloc[0]["cell_id"])

        dpa_mean, dpa_std = _get(dfs["dpa"], cell_id)
        he_mean,  he_std  = _get(dfs["He"],  cell_id)
        h_mean,   h_std   = _get(dfs["H"],   cell_id)

        rows.append(dict(
            chunk_key       = chunk_key,
            layer           = layer_label,
            cell_id         = cell_id,
            dpa_fpy         = dpa_mean,
            dpa_fpy_std     = dpa_std,
            H_appm_fpy      = h_mean,
            H_appm_fpy_std  = h_std,
            He_appm_fpy     = he_mean,
            He_appm_fpy_std = he_std,
        ))

        # inline print
        #print(
        #    f"  [{chunk_key}] {layer_label:<14}  cell={cell_id}  "
        #    f"DPA={dpa_mean:.3f}±{dpa_std:.4f}  "
        #    f"H={h_mean:.3f}±{h_std:.4f} appm/fpy  "
        #    f"He={he_mean:.3f}±{he_std:.4f} appm/fpy"
        #)

    if not rows:
        raise RuntimeError(
            f"No rows found for any target layer {TARGET_LAYERS} "
            f"in {chunk_key} profile CSVs. "
            f"Check that process_chunk() completed successfully."
        )

    # Build and save DataFrame
    col_order = [
        "chunk_key", "layer", "cell_id",
        "dpa_fpy",      "dpa_fpy_std",
        "H_appm_fpy",   "H_appm_fpy_std",
        "He_appm_fpy",  "He_appm_fpy_std",
    ]
    df = pd.DataFrame(rows)
    df = df[[c for c in col_order if c in df.columns]]

    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"[completed] summary neutronics results → {out_csv}")
    return df

# ──────────────────────────────────────────────────────────────────────────────
# Gas-production breakdown
# ──────────────────────────────────────────────────────────────────────────────
def _save_nuclide_breakdown(
    df: pd.DataFrame,
    col_total: str,
    gas_or_dpa: str,
    cid: int,
    outdir: Path,
) -> None:
    """Save nuclide CSV, element CSV, and print console summary."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Determine percentage column 
    pct_col = next((c for c in df.columns if c.endswith("_pct") and "cumulative" not in c), None)

    elem_summary = (
        df.groupby("element", sort=False)
          .agg(val_total=(col_total, "sum"),
               val_pct=(pct_col, "sum") if pct_col else (col_total, "sum"))
          .sort_values("val_total", ascending=False)
          .reset_index()
    )

    N_total    = float(df["N_total_struct"].iloc[0]) if "N_total_struct" in df.columns else float("nan")
    factor     = float(df["factor"].iloc[0])         if "factor"         in df.columns else float("nan")
    total_val  = float(df[col_total].sum())
    units      = "appm/fpy" if "appm" in col_total else "DPA/fpy"

    #inline print
    #print(f"\n{'='*70}")
    #print(f"  {gas_or_dpa} diagnostic  |  cell {cid}")
    #print(f"{'='*70}")
    #print(f"  N_total_struct         = {N_total:.6e} atoms")
    #print(f"  factor / source_scale  = {factor:.6e}")
    #print(f"  Total {gas_or_dpa:<6}           = {total_val:.6f} {units}")
    #print(f"\n  {'Element':<8} {gas_or_dpa:>14}  {'%':>7}")
    #print(f"  {'-'*8}  {'-'*14}  {'-'*7}")
    #for _, r in elem_summary.iterrows():
    #    print(f"  {r['element']:<8} {r['val_total']:>14.4f}  {r['val_pct']:>6.2f}%")

    nuclide_path = outdir / f"{gas_or_dpa.lower()}_nuclide.csv"
    df.to_csv(nuclide_path, index=False)
    #print(f"[saved] {nuclide_path}")

    element_path = outdir / f"{gas_or_dpa.lower()}_element.csv"
    elem_summary.to_csv(element_path, index=False)
    #print(f"[saved] {element_path}")

def compute_gas_cell_detailed(
    cid: int,
    t: openmc.Tally,
    *,
    he3_score: str,
    he4_score: str,
    h1_score: str,
    h2_score: str,
    h3_score: str,
) -> Dict[str, object]:
    """
    Compute He and H appm/fpy for a single cell, returning both
    scalar totals and per-nuclide breakdown DataFrames.

    Returns
    -------
    {
        "he_appm":     float,
        "he_appm_std": float,
        "h_appm":      float,
        "h_appm_std":  float,
        "df_He":       pd.DataFrame,   # per-nuclide He breakdown
        "df_H":        pd.DataFrame,   # per-nuclide H breakdown
    }
    """
    cid          = int(cid)
    f_struct     = cell_struct_origin_frac.get(cid, {}) or {}
    N_total      = float(cell_total_atoms_struct.get(cid, 0.0))
    nuclides     = list(t.nuclides or [])
    factor       = neutron_source_rate * s_in_y / N_total * 1e6 if N_total > 0 else 0.0
    struct_atoms = cell_struct_nuclide_atoms.get(cid, {}) or {}

    gas_configs = {
        "He": {
            "scores":    [he3_score,          he4_score         ],
            "col_raws":  [("he3_raw","he3_raw_std"), ("he4_raw","he4_raw_std")],
            "col_total": "he_appm_contrib",
        },
        "H": {
            "scores":    [h1_score,   h2_score,   h3_score  ],
            "col_raws":  [("h1_raw","h1_raw_std"), ("h2_raw","h2_raw_std"),
                          ("h3_raw","h3_raw_std")],
            "col_total": "h_appm_contrib",
        },
    }

    out = {}
    for gas, gcfg in gas_configs.items():
        rows = []
        total_mean = 0.0
        total_var  = 0.0

        for nuc in nuclides:
            nuc = str(nuc)
            f   = float(f_struct.get(nuc, 0.0))
            N   = float(struct_atoms.get(nuc, 0.0))
            raw_vals      = {}
            total_weighted = 0.0

            for (col_m, col_s), score in zip(gcfg["col_raws"], gcfg["scores"]):
                try:
                    m, s = tally_mean_std(t, score=score, nuclide=nuc)
                except Exception:
                    m, s = 0.0, 0.0
                raw_vals[col_m]  = m
                raw_vals[col_s]  = s
                total_weighted  += f * m

            appm_contrib  = total_weighted * factor
            total_mean   += appm_contrib
            # propagate std: recompute weighted std per score
            for (col_m, col_s), score in zip(gcfg["col_raws"], gcfg["scores"]):
                ws = f * raw_vals.get(col_s, 0.0) * factor
                total_var += ws ** 2

            row = dict(
                nuclide        = nuc,
                element        = element_from_nuclide(nuc),
                N_atoms_struct = N,
                f_struct       = f,
                **raw_vals,
                weighted_total  = total_weighted,
                N_total_struct  = N_total,
                factor          = factor,
            )
            row[gcfg["col_total"]] = appm_contrib
            rows.append(row)

        col_total  = gcfg["col_total"]
        df         = pd.DataFrame(rows)
        df         = df.sort_values(col_total, ascending=False).reset_index(drop=True)
        total_appm = df[col_total].sum()
        pct_col    = f"{gas.lower()}_appm_pct"
        df[pct_col]         = df[col_total] / total_appm * 100.0 if total_appm > 0 else 0.0
        df["cumulative_pct"] = df[pct_col].cumsum()

        out[f"{gas.lower()}_appm"]     = total_mean
        out[f"{gas.lower()}_appm_std"] = math.sqrt(max(total_var, 0.0))
        out[f"df_{gas}"]               = df

    return out

def compute_dpa_cell_detailed(
    cid: int,
    t: openmc.Tally,
) -> Dict[str, object]:
    """
    Compute NRT-dpa/fpy for a single cell, returning both
    scalar total and per-nuclide breakdown DataFrame.

    Returns
    -------
    {
        "dpa":     float,
        "dpa_std": float,
        "df_dpa":  pd.DataFrame,
    }
    """
    cid          = int(cid)
    f_struct     = cell_struct_origin_frac.get(cid, {}) or {}
    N_total      = float(cell_total_atoms_struct.get(cid, 0.0))
    nuclides     = list(t.nuclides or [])
    source_scale = neutron_source_rate * s_in_y
    struct_atoms = cell_struct_nuclide_atoms.get(cid, {}) or {}

    rows     = []
    sum_mean = 0.0
    sum_var  = 0.0

    for nuc in nuclides:
        nuc = str(nuc)
        f   = float(f_struct.get(nuc, 0.0))
        N   = float(struct_atoms.get(nuc, 0.0))

        try:
            dmg_m, dmg_s = tally_mean_std(t, score="damage-energy", nuclide=nuc)
        except Exception:
            dmg_m, dmg_s = 0.0, 0.0

        if f == 0.0 or N_total <= 0.0:
            dpa_contrib = 0.0
            Ed         = float("nan")
            disp_scale = float("nan")
        else:
            el          = element_from_nuclide(nuc)
            Ed          = float(geo.materials.Ed(el))
            disp_scale  = 0.8 / (2.0 * Ed)
            dpa_contrib = disp_scale * source_scale * f * dmg_m / N_total
            dpa_std_c   = disp_scale * source_scale * f * dmg_s  / N_total
            sum_mean   += dpa_contrib
            sum_var    += dpa_std_c ** 2

        rows.append(dict(
            nuclide        = nuc,
            element        = element_from_nuclide(nuc),
            N_atoms_struct = N,
            f_struct       = f,
            dmg_raw        = dmg_m,
            dmg_raw_std    = dmg_s,
            Ed_eV          = Ed,
            disp_scale     = disp_scale,
            N_total_struct = N_total,
            dpa_contrib    = dpa_contrib,
        ))

    df        = pd.DataFrame(rows)
    df        = df.sort_values("dpa_contrib", ascending=False).reset_index(drop=True)
    total_dpa = df["dpa_contrib"].sum()
    df["dpa_pct"]        = df["dpa_contrib"] / total_dpa * 100.0 if total_dpa > 0 else 0.0
    df["cumulative_pct"] = df["dpa_pct"].cumsum()

    return {
        "dpa":     sum_mean,
        "dpa_std": math.sqrt(max(sum_var, 0.0)),
        "df_dpa":  df,
    }

# ──────────────────────────────────────────────────────────────────────────────
# Core per-chunk processor
# ──────────────────────────────────────────────────────────────────────────────
def process_chunk(
    sp: openmc.StatePoint,
    chunk_key: str,
    cell_ids: List[int],
    *,
    dpa_gas_map: Dict[int, openmc.Tally],
    xcentroids: np.ndarray,
    xedges: np.ndarray,
    do_albedo: bool = False,
) -> None:
    outdir = RESULTS_DIR / chunk_key
    outdir.mkdir(parents=True, exist_ok=True)

    profiles_dir  = outdir / "profiles"
    spectra_dir   = outdir / "spectra"
    plots_dir     = outdir / "plots"
    albedo_dir    = outdir / "albedo"
    breakdown_dir = outdir / "nuclide_breakdown"

    for d in (profiles_dir, spectra_dir, plots_dir, albedo_dir, breakdown_dir):
        d.mkdir(parents=True, exist_ok=True)

    is_ob      = chunk_key.startswith("OB_")
    region_tag = "OB" if is_ob else "IB"

    xcent  = np.asarray(xcentroids, float)
    xedges = np.asarray(xedges, float)

    labels  = layer_names_for_chunk(cell_ids, region_tag)
    scaling = scaling_for_cells(cell_ids)
    base_df = make_chunk_base_df(
        chunk_key=chunk_key, cell_ids=cell_ids,
        labels=labels, xcent=xcent, xedges=xedges,
    )
    write_struct_origin_csv(outdir, chunk_key, cell_ids, labels,
                            cell_struct_origin_frac)

    # ==========================================
    # Neutron & photon spectra (flux_tally)
    # ==========================================
    t_flux_id         = require_tally_id(sp, flux_tally, "flux spectrum tally")
    t_flux            = sp.get_tally(id=t_flux_id)
    cell_bins_flux    = get_cell_bins(t_flux)
    flux_mean, flux_std = get_flux_spectrum_arrays(t_flux)

    neutron_flux          = flux_mean[:, 0, :]
    neutron_flux_std      = flux_std[:, 0, :]
    photon_flux           = flux_mean[:, 1, :]
    photon_flux_std       = flux_std[:, 1, :]

    neutron_flux_chunk     = subset_by_cells(neutron_flux,     cell_bins_flux, cell_ids)
    neutron_flux_std_chunk = subset_by_cells(neutron_flux_std, cell_bins_flux, cell_ids)
    photon_flux_chunk      = subset_by_cells(photon_flux,      cell_bins_flux, cell_ids)
    photon_flux_std_chunk  = subset_by_cells(photon_flux_std,  cell_bins_flux, cell_ids)

    # Neutron spectrum plot
    fig, ax = plt.subplots()
    for i, cid in enumerate(cell_ids):
        fl = neutron_flux_chunk[i].flatten() * scaling[int(cid)] / unit_lethargy
        sl = neutron_flux_std_chunk[i].flatten() * scaling[int(cid)] / unit_lethargy
        ax.loglog(energies[:-1], fl,
                  label=f"{labels[i]} ({xcent[i]:.1f} cm)", color=colors[i])
        ax.fill_between(energies[:-1], np.maximum(fl - sl, 1e-99), fl + sl,
                        alpha=0.25, color=colors[i], lw=0)
    ax.set(xlim=[1e-3, 1e8], ylim=[1e6, 1e17],
           xlabel="Energy [eV]",
           ylabel="Neutron flux per unit lethargy [1/cm²/s]")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.legend(fontsize=8, ncol=2)
    fig.savefig(spectra_dir / f"n_flux_spectrum_{chunk_key}.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)


    # Photon spectrum plot
    fig, ax = plt.subplots()
    for i, cid in enumerate(cell_ids):
        fl = photon_flux_chunk[i].flatten() * scaling[int(cid)] / unit_lethargy
        sl = photon_flux_std_chunk[i].flatten() * scaling[int(cid)] / unit_lethargy
        ax.loglog(energies[:-1], fl,
                  label=f"{labels[i]} ({xcent[i]:.1f} cm)", color=colors[i])
        ax.fill_between(energies[:-1], np.maximum(fl - sl, 1e-99), fl + sl,
                        alpha=0.25, color=colors[i], lw=0)
    ax.set(xlim=[1e3, 1e8], ylim=[1e6, 1e17],
           xlabel="Energy [eV]",
           ylabel="Photon flux per unit lethargy [1/cm²/s]")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.legend(fontsize=8, ncol=2)
    fig.savefig(spectra_dir / f"p_flux_spectrum_{chunk_key}.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Spectrum CSVs
    for name, fl_arr, fl_std in [
        ("neutron", neutron_flux_chunk, neutron_flux_std_chunk),
        ("photon",  photon_flux_chunk,  photon_flux_std_chunk),
    ]:
        df = pd.DataFrame({"E_mid_eV": E_mid})
        for i, lab in enumerate(labels):
            s             = scaling[int(cell_ids[i])]
            df[lab]       = fl_arr[i].flatten() * s / unit_lethargy
            df[f"{lab}_std"] = fl_std[i].flatten() * s / unit_lethargy
        df.to_csv(spectra_dir / f"{name}_spectrum_{chunk_key}.csv", index=False)

    # ==========================================
    # Total flux reconstructed from spectrum
    # ==========================================
    def _total(arr, sarr, cids):
        mn  = np.array([arr[i].sum() * scaling[int(c)]
                        for i, c in enumerate(cids)], dtype=float)
        std = np.array([np.sqrt((sarr[i] ** 2).sum()) * scaling[int(c)]
                        for i, c in enumerate(cids)], dtype=float)
        return mn, std

    direct_total_neut, direct_total_neut_std = _total(
        neutron_flux_chunk, neutron_flux_std_chunk, cell_ids
    )
    direct_total_phot, direct_total_phot_std = _total(
        photon_flux_chunk, photon_flux_std_chunk, cell_ids
    )

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    for y, s, lab in [
        (direct_total_neut, direct_total_neut_std, "neutron"),
        (direct_total_phot, direct_total_phot_std, "photon"),
    ]:
        lo = np.maximum(y - s, 1e-30)
        hi = y + s
        ax.step(xedges, np.r_[y, y[-1]], where="post", label=lab)
        ax.fill_between(xedges, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]],
                        step="post", alpha=0.15)
    ax.set(xlabel="Radial Position [cm]", ylabel="Total Flux [1/cm²/s]")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.legend()
    fig.savefig(plots_dir / f"flux_total_{chunk_key}.png",
                dpi=300, bbox_inches="tight")

    if _ref_data:
        _overlay_ref(ax, chunk_key, "flux_total_neutron",
                     label="neutron (W Armor + Eurofer Structural)", color="tab:blue", ls=(0, (5, 2)))
        _overlay_ref(ax, chunk_key, "flux_total_photon",
                     label="photon (W Armor + Eurofer Structural)",  color="tab:orange", ls=(0, (5, 2)))
        handles, _ = ax.get_legend_handles_labels()
        if handles:
            ax.legend()
        fig.savefig(plots_dir / f"flux_total_{chunk_key}_vs_ref.png",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)

    save_profile(base_df, profiles_dir, quantity="flux_total_neutron",
                 mean=direct_total_neut, std=direct_total_neut_std, units="1/cm2/s")
    save_profile(base_df, profiles_dir, quantity="flux_total_photon",
                 mean=direct_total_phot, std=direct_total_phot_std, units="1/cm2/s")

    # ==========================================
    # Heating
    # ==========================================
    t_heat_id      = require_tally_id(sp, heating_tally, "heating tally")
    t_heat         = sp.get_tally(id=t_heat_id)
    cell_bins_heat = get_cell_bins(t_heat)
    heat_mean      = t_heat.get_values(value="mean").flatten()
    heat_std       = t_heat.get_values(value="std_dev").flatten()

    heat_mean_chunk = subset_by_cells(heat_mean, cell_bins_heat, cell_ids)
    heat_std_chunk  = subset_by_cells(heat_std,  cell_bins_heat, cell_ids)

    heat_w = np.array(
        [heat_mean_chunk[i] * scaling[int(cid)] * ev_to_joule
         for i, cid in enumerate(cell_ids)], dtype=float,
    )
    heat_w_std = np.array(
        [heat_std_chunk[i] * scaling[int(cid)] * ev_to_joule
         for i, cid in enumerate(cell_ids)], dtype=float,
    )

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    lo = np.maximum(heat_w - heat_w_std, 1e-30)
    hi = heat_w + heat_w_std
    ax.step(xedges, np.r_[heat_w, heat_w[-1]], where="post", label="Heating")
    ax.fill_between(xedges, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]],
                    step="post", alpha=0.3)
    ax.set(xlabel="Radial Position [cm]", ylabel="Heating [W/cm³]")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    fig.savefig(plots_dir / f"heating_{chunk_key}.png",
                dpi=300, bbox_inches="tight")

    if _ref_data:
        _overlay_ref(ax, chunk_key, "heating",
                     color="black", ls=(0, (5, 2)))
        handles, _ = ax.get_legend_handles_labels()
        if handles:
            ax.legend()
        fig.savefig(plots_dir / f"heating_{chunk_key}_vs_ref.png",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)

    save_profile(base_df, profiles_dir, quantity="heating",
                 mean=heat_w, std=heat_w_std, units="W/cm3")



    # ==========================================
    # H, He, DPA — computed once per cell,
    # nuclide breakdown saved automatically
    # ==========================================
    h_appm_y      = np.zeros(len(cell_ids), dtype=float)
    h_appm_y_std  = np.zeros(len(cell_ids), dtype=float)
    he_appm_y     = np.zeros(len(cell_ids), dtype=float)
    he_appm_y_std = np.zeros(len(cell_ids), dtype=float)
    dpa_y         = np.zeros(len(cell_ids), dtype=float)
    dpa_y_std     = np.zeros(len(cell_ids), dtype=float)

    for i, cid in enumerate(cell_ids):
        cid = int(cid)
        t   = dpa_gas_map.get(cid)
        if t is None:
            continue

        scores = {str(s) for s in (t.scores or [])}
        if GAS_SCORES_EXPLICIT.issubset(scores):
            he3_score, he4_score         = "He3-production", "He4-production"
            h1_score, h2_score, h3_score = "H1-production", "H2-production", "H3-production"
        elif GAS_SCORES_REACTION.issubset(scores):
            he3_score, he4_score         = "(n,X3He)", "(n,Xa)"
            h1_score, h2_score, h3_score = "(n,Xp)", "(n,Xd)", "(n,Xt)"
        else:
            continue

        cell_breakdown_dir = breakdown_dir / f"cell_{cid}"

        # Gas (He + H) 
        gas_result = compute_gas_cell_detailed(
            cid, t,
            he3_score=he3_score, he4_score=he4_score,
            h1_score=h1_score,   h2_score=h2_score, h3_score=h3_score,
        )
        h_appm_y[i]      = gas_result["h_appm"]
        h_appm_y_std[i]  = gas_result["h_appm_std"]
        he_appm_y[i]     = gas_result["he_appm"]
        he_appm_y_std[i] = gas_result["he_appm_std"]

        for gas, col in [("He", "he_appm_contrib"), ("H", "h_appm_contrib")]:
            _save_nuclide_breakdown(gas_result[f"df_{gas}"], col, gas, cid, cell_breakdown_dir)

        # DPA 
        dpa_result   = compute_dpa_cell_detailed(cid, t)
        dpa_y[i]     = dpa_result["dpa"]
        dpa_y_std[i] = dpa_result["dpa_std"]

        _save_nuclide_breakdown(dpa_result["df_dpa"], "dpa_contrib", "DPA", cid, cell_breakdown_dir)

    # ==========================================
    # H (H1 + H2 + H3) production plot + CSV
    # ==========================================
    lo = np.maximum(h_appm_y - h_appm_y_std, 1e-30)
    hi = h_appm_y + h_appm_y_std

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.step(xedges, np.r_[h_appm_y, h_appm_y[-1]], where="post")
    ax.fill_between(xedges, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]],
                    step="post", alpha=0.3)
    ax.set(xlabel="Radial Position [cm]", ylabel="H [appm/fpy]")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    fig.savefig(plots_dir / f"h1_{chunk_key}.png", dpi=300, bbox_inches="tight")

    if _ref_data:
        _overlay_ref(ax, chunk_key, "H_appm_fpy_struct_origin",
                    color="black", ls=(0, (5, 2)))
        handles, _ = ax.get_legend_handles_labels()
        if handles:
            ax.legend()
        fig.savefig(plots_dir / f"h1_{chunk_key}_vs_ref.png",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)

    save_profile(base_df, profiles_dir, quantity="H_appm_fpy_struct_origin",
                 mean=h_appm_y, std=h_appm_y_std, units="appm/fpy")



    # ==========================================
    # He (He3 + He4) production plot + CSV
    # ==========================================
    lo = np.maximum(he_appm_y - he_appm_y_std, 1e-30)
    hi = he_appm_y + he_appm_y_std

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.step(xedges, np.r_[he_appm_y, he_appm_y[-1]], where="post")
    ax.fill_between(xedges, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]],
                    step="post", alpha=0.3)
    ax.set(xlabel="Radial Position [cm]", ylabel="He [appm/fpy]")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    fig.savefig(plots_dir / f"he_{chunk_key}.png", dpi=300, bbox_inches="tight")

    if _ref_data:
        _overlay_ref(ax, chunk_key, "He_appm_fpy_struct_origin",
                        color="black", ls=(0, (5, 2)))
        handles, _ = ax.get_legend_handles_labels()
        if handles:
            ax.legend()
        fig.savefig(plots_dir / f"he_{chunk_key}_vs_ref.png",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)

    save_profile(base_df, profiles_dir, quantity="He_appm_fpy_struct_origin",
                 mean=he_appm_y, std=he_appm_y_std, units="appm/fpy")



    # ==========================================
    # DPA plot + CSV
    # ==========================================
    lo = np.maximum(dpa_y - dpa_y_std, 1e-30)
    hi = dpa_y + dpa_y_std

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.step(xedges, np.r_[dpa_y, dpa_y[-1]], where="post")
    ax.fill_between(xedges, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]],
                    step="post", alpha=0.3)
    ax.set(xlabel="Radial Position [cm]", ylabel="NRT-dpa/fpy")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    fig.savefig(plots_dir / f"dpa_{chunk_key}.png", dpi=300, bbox_inches="tight")

    if _ref_data:
        _overlay_ref(ax, chunk_key, "dpa_fpy_struct_origin",
                            color="black", ls=(0, (5, 2)))
        handles, _ = ax.get_legend_handles_labels()
        if handles:
            ax.legend()
        fig.savefig(plots_dir / f"dpa_{chunk_key}_vs_ref.png",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
 
    save_profile(base_df, profiles_dir, quantity="dpa_fpy_struct_origin",
                 mean=dpa_y, std=dpa_y_std, units="DPA/fpy")

    if do_albedo:
        compute_albedo_for_chunk(
            sp,
            cell_ids=cell_ids,
            labels=labels,
            outdir=albedo_dir,
            chunk_key=chunk_key,
            excluded_surfaces=set(),
        )

    print(f"[completed] {chunk_key} → {outdir}")

# ──────────────────────────────────────────────────────────────────────────────
# Poloidal layer sweep  (tokamak only)
# ──────────────────────────────────────────────────────────────────────────────
def _layer_cells_region1(
    layer_index: int,
    order: str = "OB_THEN_IB",
    ) -> Tuple[List[int], List[str]]:
    keys_ob = [f"OB_1_b{b}" for b in range(1, n_breeder + 1)]
    keys_ib = [f"IB_1_b{b}" for b in range(1, n_breeder + 1)]

    if order == "OB_THEN_IB":
        keys = keys_ob + keys_ib
    elif order == "IB_THEN_OB":
        keys = keys_ib + keys_ob
    else:
        raise ValueError("order must be 'OB_THEN_IB' or 'IB_THEN_OB'")

    cell_ids_layer: List[int] = []
    chunk_labels:   List[str] = []

    for k in keys:
        chunk = ob_by_key.get(k) if k.startswith("OB_") else ib_by_key.get(k)
        if not chunk or layer_index >= len(chunk):
            continue
        cell_ids_layer.append(int(chunk[layer_index]))
        chunk_labels.append(k)

    return cell_ids_layer, chunk_labels

def compute_dpa_layer_region1_only(
    cell_ids_layer: List[int],
    *,
    dpa_gas_map: Dict[int, openmc.Tally],
    ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute NRT-dpa/fpy for each cell in cell_ids_layer.
    Returns (dpa_y, dpa_y_std) arrays of length n = len(cell_ids_layer).
    """
    n       = len(cell_ids_layer)
    idx_map = {int(cid): i for i, cid in enumerate(cell_ids_layer)}

    dpa_y     = np.zeros(n, dtype=float)
    dpa_y_std = np.zeros(n, dtype=float)

    source_scale = float(neutron_source_rate) * float(s_in_y)

    for cid in cell_ids_layer:
        cid = int(cid)
        t   = dpa_gas_map.get(cid)
        if t is None:
            continue

        denom_atoms = float(cell_total_atoms_struct.get(cid, 0.0))
        if denom_atoms <= 0.0:
            continue

        f_struct = cell_struct_origin_frac.get(cid, {}) or {}
        nuclides = list(t.nuclides or [])

        sum_mean = 0.0
        sum_var  = 0.0

        for nuc in nuclides:
            nuc = str(nuc)
            try:
                dmg_mean, dmg_std = tally_mean_std(t, score="damage-energy", nuclide=nuc)
            except Exception:
                continue

            fs = float(f_struct.get(nuc, 0.0))
            if fs == 0.0:
                continue
            dmg_mean *= fs
            dmg_std  *= fs

            el         = element_from_nuclide(nuc)
            Ed         = float(geo.materials.Ed(el))
            disp_scale = 0.8 / (2.0 * Ed)

            disp_mean = disp_scale * source_scale * dmg_mean
            disp_std  = disp_scale * source_scale * dmg_std

            dpa_mean = disp_mean / denom_atoms
            dpa_std  = disp_std  / denom_atoms

            sum_mean += dpa_mean
            sum_var  += dpa_std ** 2

        i             = idx_map[cid]
        dpa_y[i]     = sum_mean
        dpa_y_std[i] = math.sqrt(max(sum_var, 0.0))

    return dpa_y, dpa_y_std

def _empty_poloidal_plot(x_data, num_pts):
    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    split  = n_breeder - 0.5
    xticks = np.arange(num_pts)

    ax.set_yscale("linear")
    ax.set_xlim([-0.5, len(x_data) - 2.5])
    ax.axvline(split, linestyle="--", linewidth=1, color="k")
    ax.text(
        split - 0.15, 0.95,
        "$\\leftarrow$ outboard   inboard $\\rightarrow$",
        transform=ax.get_xaxis_transform(),
        ha="center", va="top", fontsize=10,
    )
    ax.set_xticks(xticks)
    ax.set_xticklabels(
        [f"OB{i+1}" if i < n_breeder else f"IB{i+1-n_breeder}"
         for i in range(num_pts)],
        rotation=45, ha="right", fontsize=9,
    )
    ax.set_xlabel("Poloidal regions", fontsize=11)
    ax.grid(True, which="major", linestyle="--", linewidth=0.6, alpha=0.8)
    ax.yaxis.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.5)
    ax.minorticks_on()
    ax.tick_params(axis="x", which="minor", bottom=False)
    fig.tight_layout()

    return fig, ax

def _pad(data: np.ndarray) -> np.ndarray:
    return np.concatenate([[data[0]], data, [data[-1]]])

def process_layer_region1_only(
    sp: openmc.StatePoint,
    layer_index: int,
    layer_tag: str,
    outdir: Path,
    *,
    dpa_gas_map: Dict[int, openmc.Tally],
    ) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Subdirectories
    profiles_dir = outdir / "profiles"
    plots_dir    = outdir / "plots"
    for d in (profiles_dir, plots_dir):
        d.mkdir(parents=True, exist_ok=True)

    cell_ids_layer, _chunk_labels = _layer_cells_region1(layer_index)
    n = len(cell_ids_layer)
    if n == 0:
        print(f"[warn] layer_index={layer_index} ({layer_tag}): "
              f"no cells found, skipping.")
        return

    scaling = scaling_for_cells(cell_ids_layer)
    x = _pad(np.arange(n, dtype=float))

    base_df = make_layer_base_df(
        layer_tag=layer_tag,
        cell_ids=cell_ids_layer,
        chunk_labels=_chunk_labels,
        poloidal_index=list(range(n)),
    )

    def _plot_step(y, ys, ylabel, fname, yscale="auto", save_dir=None):
        y  = np.asarray(y,  float)
        ys = np.asarray(ys, float)
        yp  = _pad(y)
        ysp = _pad(ys)
        lo  = np.clip(yp - ysp, 1e-99, None)
        hi  = yp + ysp
        fig, ax = _empty_poloidal_plot(x, n)
        set_ylim_and_ticks(ax, yp, scale=yscale)
        ax.set_yscale("log" if yscale == "log" else ax.get_yscale())
        ax.step(x, yp, where="mid", linewidth=2)
        ax.fill_between(x, lo, hi, step="mid", alpha=0.15)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(f"{ylabel} — {layer_tag}", fontsize=12)
        save_to = (save_dir or plots_dir) / fname
        fig.savefig(save_to, dpi=300, bbox_inches="tight")
        return fig, ax

    # ==========================================
    # Total flux reconstructed from spectrum
    # ==========================================
    t_flux         = sp.get_tally(id=require_tally_id(sp, flux_tally, "flux spectrum tally"))
    cell_bins_flux = get_cell_bins(t_flux)
    flux_mean, flux_std = get_flux_spectrum_arrays(t_flux)

    neutron_flux           = flux_mean[:, 0, :]
    neutron_flux_std       = flux_std[:, 0, :]
    photon_flux            = flux_mean[:, 1, :]
    photon_flux_std        = flux_std[:, 1, :]

    neutron_flux_layer     = subset_by_cells(neutron_flux,     cell_bins_flux, cell_ids_layer)
    neutron_flux_std_layer = subset_by_cells(neutron_flux_std, cell_bins_flux, cell_ids_layer)
    photon_flux_layer      = subset_by_cells(photon_flux,      cell_bins_flux, cell_ids_layer)
    photon_flux_std_layer  = subset_by_cells(photon_flux_std,  cell_bins_flux, cell_ids_layer)

    neut     = np.array([neutron_flux_layer[i].sum() * scaling[int(c)]
                         for i, c in enumerate(cell_ids_layer)], dtype=float)
    neut_std = np.array([np.sqrt((neutron_flux_std_layer[i] ** 2).sum()) * scaling[int(c)]
                         for i, c in enumerate(cell_ids_layer)], dtype=float)
    phot     = np.array([photon_flux_layer[i].sum() * scaling[int(c)]
                         for i, c in enumerate(cell_ids_layer)], dtype=float)
    phot_std = np.array([np.sqrt((photon_flux_std_layer[i] ** 2).sum()) * scaling[int(c)]
                         for i, c in enumerate(cell_ids_layer)], dtype=float)

    save_profile(base_df, profiles_dir, quantity="flux_total_neutron",
                 mean=neut, std=neut_std, units="1/cm2/s")
    save_profile(base_df, profiles_dir, quantity="flux_total_photon",
                 mean=phot, std=phot_std, units="1/cm2/s")

    fig, ax = _empty_poloidal_plot(x, n)
    set_ylim_and_ticks(ax, np.r_[_pad(neut), _pad(phot)], scale="auto")
    ax.step(x, _pad(np.clip(neut, 1e-99, None)),
            where="mid", linewidth=2, label="neutron")
    ax.step(x, _pad(np.clip(phot, 1e-99, None)),
            where="mid", linewidth=2, label="photon")
    ax.fill_between(x, _pad(np.clip(neut - neut_std, 1e-99, None)),
                    _pad(neut + neut_std), step="mid", alpha=0.15)
    ax.fill_between(x, _pad(np.clip(phot - phot_std, 1e-99, None)),
                    _pad(phot + phot_std), step="mid", alpha=0.15)
    ax.set_ylabel("Total Flux [1/cm²/s]", fontsize=11)
    ax.set_title(f"Total Flux — {layer_tag}", fontsize=12)
    ax.legend()
    fig.savefig(plots_dir / f"flux_total_{layer_tag}.png",
                dpi=300, bbox_inches="tight")
                
    if _ref_data:
        _overlay_ref_poloidal(ax, layer_tag, "flux_total_neutron",
                            label="neutron (W Armor + Eurofer Structural)", color="tab:blue", ls=(0, (3, 1)))
        _overlay_ref_poloidal(ax, layer_tag, "flux_total_photon",
                            label="photon (W Armor + Eurofer Structural)", color="tab:orange", ls=(0, (3, 1)))
        handles, _ = ax.get_legend_handles_labels()
        if handles:
            ax.legend(fontsize=8)
        fig.savefig(plots_dir / f"flux_total_{layer_tag}_vs_ref.png",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ==========================================
    # Heating
    # ==========================================
    t_heat         = sp.get_tally(id=require_tally_id(sp, heating_tally, "heating tally"))
    cell_bins_heat = get_cell_bins(t_heat)
    heat_mean      = t_heat.get_values(value="mean").flatten()
    heat_std       = t_heat.get_values(value="std_dev").flatten()

    heat_mean_layer = subset_by_cells(heat_mean, cell_bins_heat, cell_ids_layer)
    heat_std_layer  = subset_by_cells(heat_std,  cell_bins_heat, cell_ids_layer)

    heat_w = np.array(
        [heat_mean_layer[i] * scaling[int(c)] * ev_to_joule
         for i, c in enumerate(cell_ids_layer)], dtype=float,
    )
    heat_w_std = np.array(
        [heat_std_layer[i] * scaling[int(c)] * ev_to_joule
         for i, c in enumerate(cell_ids_layer)], dtype=float,
    )

    fig, ax = _plot_step(heat_w, heat_w_std, "Heating [W/cm³]",
                         f"heating_{layer_tag}.png", save_dir=plots_dir)
    if _ref_data:
        _overlay_ref_poloidal(ax, layer_tag, "heating")
        handles, _ = ax.get_legend_handles_labels()
        if handles:
            ax.legend(fontsize=8)
        fig.savefig(plots_dir / f"heating_{layer_tag}_vs_ref.png",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)

    save_profile(base_df, profiles_dir, quantity="heating",
                 mean=heat_w, std=heat_w_std, units="W/cm3")

    # ==========================================
    # H production (H1 + H2 + H3)
    # ==========================================
    h_appm_y     = np.zeros(n, dtype=float)
    h_appm_y_std = np.zeros(n, dtype=float)

    for i, cid in enumerate(cell_ids_layer):
        cid = int(cid)
        t   = dpa_gas_map.get(cid)
        if t is None:
            continue
        scores = {str(s) for s in (t.scores or [])}
        if GAS_SCORES_EXPLICIT.issubset(scores):
            h1_score, h2_score, h3_score = (
                "H1-production", "H2-production", "H3-production"
            )
        elif GAS_SCORES_REACTION.issubset(scores):
            h1_score, h2_score, h3_score = "(n,Xp)", "(n,Xd)", "(n,Xt)"
        else:
            continue

        f_struct = cell_struct_origin_frac.get(cid, {}) or {}
        nuclides = list(t.nuclides or [])

        m1, s1, _ = corrected_sum_detailed(t, score=h1_score, nuclides=nuclides, f_struct=f_struct)
        m2, s2, _ = corrected_sum_detailed(t, score=h2_score, nuclides=nuclides, f_struct=f_struct)
        m3, s3, _ = corrected_sum_detailed(t, score=h3_score, nuclides=nuclides, f_struct=f_struct)

        h_mean_corr = m1 + m2 + m3
        h_std_corr  = math.sqrt(s1 ** 2 + s2 ** 2 + s3 ** 2)

        denom = float(cell_total_atoms_struct.get(cid, 0.0))
        if denom > 0.0:
            factor           = neutron_source_rate * s_in_y / denom * 1e6
            h_appm_y[i]     = h_mean_corr * factor
            h_appm_y_std[i] = h_std_corr  * factor

    fig, ax = _plot_step(h_appm_y, h_appm_y_std,
                         "H production [appm/fpy]", f"h1_{layer_tag}.png",
                         save_dir=plots_dir)
    if _ref_data:
        _overlay_ref_poloidal(ax, layer_tag, "H_appm_fpy_struct_origin")
        handles, _ = ax.get_legend_handles_labels()
        if handles:
            ax.legend(fontsize=8)
        fig.savefig(plots_dir / f"h1_{layer_tag}_vs_ref.png",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)

    save_profile(base_df, profiles_dir, quantity="H_appm_fpy_struct_origin",
                 mean=h_appm_y, std=h_appm_y_std, units="appm/fpy")

    # ==========================================
    # He production (He3 + He4)
    # ==========================================
    he_appm_y     = np.zeros(n, dtype=float)
    he_appm_y_std = np.zeros(n, dtype=float)

    for i, cid in enumerate(cell_ids_layer):
        cid = int(cid)
        t   = dpa_gas_map.get(cid)
        if t is None:
            continue
        scores = {str(s) for s in (t.scores or [])}
        if GAS_SCORES_EXPLICIT.issubset(scores):
            he3_score, he4_score = "He3-production", "He4-production"
        elif GAS_SCORES_REACTION.issubset(scores):
            he3_score, he4_score = "(n,X3He)", "(n,Xa)"
        else:
            continue

        f_struct = cell_struct_origin_frac.get(cid, {}) or {}
        nuclides = list(t.nuclides or [])

        m3, s3, _ = corrected_sum_detailed(t, score=he3_score, nuclides=nuclides, f_struct=f_struct)
        m4, s4, _ = corrected_sum_detailed(t, score=he4_score, nuclides=nuclides, f_struct=f_struct)

        he_mean_corr = m3 + m4
        he_std_corr  = math.sqrt(s3 ** 2 + s4 ** 2)

        denom = float(cell_total_atoms_struct.get(cid, 0.0))
        if denom > 0.0:
            factor            = neutron_source_rate * s_in_y / denom * 1e6
            he_appm_y[i]     = he_mean_corr * factor
            he_appm_y_std[i] = he_std_corr  * factor

    fig, ax = _plot_step(he_appm_y, he_appm_y_std,
                         "He production [appm/fpy]", f"he_{layer_tag}.png",
                         save_dir=plots_dir)
    if _ref_data:
        _overlay_ref_poloidal(ax, layer_tag, "He_appm_fpy_struct_origin")
        handles, _ = ax.get_legend_handles_labels()
        if handles:
            ax.legend(fontsize=8)
        fig.savefig(plots_dir / f"he_{layer_tag}_vs_ref.png",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)

    save_profile(base_df, profiles_dir, quantity="He_appm_fpy_struct_origin",
                 mean=he_appm_y, std=he_appm_y_std, units="appm/fpy")

    # ==========================================
    # DPA/fpy
    # ==========================================
    dpa_y, dpa_y_std = compute_dpa_layer_region1_only(
        cell_ids_layer,
        dpa_gas_map=dpa_gas_map,
    )

    fig, ax = _plot_step(dpa_y, dpa_y_std, "NRT-dpa/fpy", f"dpa_{layer_tag}.png",
                         yscale="linear", save_dir=plots_dir)
    if _ref_data:
        _overlay_ref_poloidal(ax, layer_tag, "dpa_fpy_struct_origin")
        handles, _ = ax.get_legend_handles_labels()
        if handles:
            ax.legend(fontsize=8)
        fig.savefig(plots_dir / f"dpa_{layer_tag}_vs_ref.png",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)

    save_profile(base_df, profiles_dir, quantity="dpa_fpy_struct_origin",
                 mean=dpa_y, std=dpa_y_std, units="DPA/fpy")

    print(f"[completed] layer {layer_tag} → {outdir}")

# ──────────────────────────────────────────────────────────────────────────────
# MAIN RUN
# ──────────────────────────────────────────────────────────────────────────────
with openmc.StatePoint(str(STATEPOINT_FILE)) as sp:
    dpa_gas_map = build_dpa_gas_map(sp)

    if cfg.SIM_TYPE in ("slab", "fast_slab"):
        keys_to_process = [cfg.ALBEDO_CHUNK_KEY]
    else:
        keys_to_process = cfg.NEUTRONICS_KEYS_TO_PROCESS

    do_albedo = (cfg.DO_ALBEDO and cfg.SIM_TYPE == "tokamak")

    for key in keys_to_process:
        if key in ob_by_key:
            chunk_cells = ob_by_key[key]
        elif key in ib_by_key:
            chunk_cells = ib_by_key[key]
        else:
            print(f"[skipped] Unknown chunk key: {key} (skipping)")
            continue

        try:
            xcent, _, xedges = radial_bins_for_key(key)
        except Exception as e:
            print(f"[skipped] No radial bins for {key}: {e} (skipping)")
            continue

        process_chunk(
            sp, key, chunk_cells,
            dpa_gas_map=dpa_gas_map,
            xcentroids=xcent,
            xedges=xedges,
            do_albedo=do_albedo,
        )

    save_summary_neutronics_results(
        results_dir = RESULTS_DIR,
        chunk_key   = cfg.TRIGGER_CHUNK_KEY,
        out_csv     = cfg.SUMMARY_RESULTS_CSV,
    )
      
    if cfg.SIM_TYPE == "tokamak":
        for layer_tag, layer_index in cfg.POLOIDAL_LAYERS:
            layer_outdir = RESULTS_DIR / f"layer_{layer_tag}_region1_only"
            process_layer_region1_only(
                sp, layer_index, layer_tag, layer_outdir,
                dpa_gas_map=dpa_gas_map,
            )
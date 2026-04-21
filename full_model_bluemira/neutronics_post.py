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
# model.xml presence check  (DROP THIS)
# ──────────────────────────────────────────────────────────────────────────────
_model_xml = cfg.NEUTRONICS_RUN_DIR / "model.xml"
if not _model_xml.is_file():
    raise FileNotFoundError(
        f"model.xml not found: {_model_xml}\n"
        f"Run neutronics_model.py first to generate it."
    )

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

if cfg.SIM_TYPE == "slab":
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

def corrected_sum(
    t: openmc.Tally,
    *,
    score: str,
    nuclides: List[str],
    f_struct: Dict[str, float],
    ) -> Tuple[float, float]:
    mean_tot = 0.0
    var_tot = 0.0

    for nuc in nuclides:
        nuc = str(nuc)
        w   = float(f_struct.get(nuc, 0.0))
        if w == 0.0:
            continue

        try:
            m, sd = tally_mean_std(t, score=score, nuclide=nuc)
        except Exception:
            continue

        mean_tot += w * m
        var_tot  += (w * sd) ** 2

    return mean_tot, math.sqrt(max(var_tot, 0.0))

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

    chunk_key = str(df["chunk_key"].iloc[0])
    path      = outdir / f"profile_{quantity}_{chunk_key}.csv"
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
    log_threshold_decades: float = 1.0,
) -> None:
    y = np.asarray(ydata, float)
    y = y[np.isfinite(y)]
    if y.size == 0:
        return

    if scale == "auto":
        ypos  = y[y > 0]
        scale = ("log" if ypos.size >= 2 and
                 (np.log10(ypos.max()) - np.log10(ypos.min()))
                 >= log_threshold_decades
                 else "linear")

    if scale == "log":
        ypos = y[y > 0]
        if ypos.size == 0:
            return
        lo, hi = np.log10(ypos.min()), np.log10(ypos.max())
        span   = max(hi - lo, 1e-12)
        pad    = max(ratio * span, 0.10)
        ax.set_ylim(max(10 ** (lo - pad), eps), 10 ** (hi + pad))
        ax.yaxis.set_major_locator(mticker.LogLocator(base=10, subs=(1., 2., 5.)))
        ax.yaxis.set_major_formatter(mticker.LogFormatterMathtext())
        ax.yaxis.set_minor_locator(
            mticker.LogLocator(base=10, subs=np.arange(1, 10) * 0.1)
        )
        ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    else:
        ymin, ymax = y.min(), y.max()
        if ymin == ymax:
            ymin, ymax = (
                (ymin * (1 - ratio), ymax * (1 + ratio))
                if ymin != 0 else (-1.0, 1.0)
            )
        else:
            span   = ymax - ymin
            ymin  -= ratio * span
            ymax  += ratio * span
        ax.set_ylim(ymin, ymax)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(6))
        ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
        ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("{x:g}"))


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
def _pydagmc_surface(sid: int):
    return geo.pydagmc_model.surfaces_by_id[int(sid)]


def _surface_area(sid: int) -> float:
    return float(_pydagmc_surface(sid).area)


def _is_surface_unshared(sid: int) -> bool:
    s = _pydagmc_surface(sid)
    try:
        return len(s.volumes) == 1
    except Exception:
        return len([v for v in getattr(s, "senses", []) if v is not None]) == 1


def _rank_by_area(sids: List[int]) -> List[int]:
    pairs = []
    for sid in sids:
        try:
            pairs.append((_surface_area(sid), sid))
        except Exception:
            pass
    return [s for _, s in sorted(pairs, reverse=True)]


def _cell_surface_ids(info_flat: dict, cid: int) -> List[int]:
    return [int(r["surface_id"]) for r in info_flat.get(int(cid), [])
            if "surface_id" in r]


def _build_special_surfaces(
    info_flat: dict,
    *,
    armor_id: int,
    breeder_back_id: int,
    vv_id: int,
    allowed: Set[int],
    ext_set: Set[int],
    int_set: Set[int],
) -> Tuple[Dict[str, int], Set[int], Dict[str, int]]:
    special: Dict[str, int] = {}
    owner:   Dict[str, int] = {}
    s_ids:   Set[int]       = set()

    def _pick(cid: int, which_ext: bool) -> Optional[int]:
        sids   = [s for s in _cell_surface_ids(info_flat, cid) if s in allowed]
        ranked = _rank_by_area(sids)
        pool   = ext_set if which_ext else int_set
        return next((s for s in ranked if s in pool), None)

    af = _pick(armor_id, True)
    ab = next(
        (s for s in _rank_by_area(
            [s for s in _cell_surface_ids(info_flat, armor_id) if s in allowed])
         if s in int_set), None
    )

    if af:
        special["Armor_front_ext"] = af
        owner["Armor_front_ext"]   = armor_id
        s_ids.add(af)
    if ab:
        special["Armor_back_int"]  = ab
        owner["Armor_back_int"]    = armor_id
        s_ids.add(ab)

    bb = _pick(breeder_back_id, True)
    if bb:
        special["Breeder_back_ext"] = bb
        owner["Breeder_back_ext"]   = breeder_back_id
        s_ids.add(bb)

    vv_ranked = _rank_by_area(
        [s for s in _cell_surface_ids(info_flat, vv_id) if s in allowed]
    )
    vv_picks = [s for s in vv_ranked if _is_surface_unshared(s)][:2]
    if len(vv_picks) >= 2:
        special["VV_face_1"] = vv_picks[0]; owner["VV_face_1"] = vv_id
        special["VV_face_2"] = vv_picks[1]; owner["VV_face_2"] = vv_id
        s_ids.add(vv_picks[0]); s_ids.add(vv_picks[1])

    return special, s_ids, owner


def compute_albedo_for_chunk(
    sp: openmc.StatePoint,
    *,
    cell_ids: List[int],
    labels: List[str],
    outdir: Path,
    chunk_key: str,
    excluded_surfaces: Optional[Set[int]] = None,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    excluded = set(excluded_surfaces or set())

    if t_current_tally is None or not p_current_tallies:
        print(f"[warn] {chunk_key}: skipping albedo – current tallies not available")
        return

    info_   = geo.info
    ext_ids = [int(x) for x in geo.external_surface_ids]
    int_ids = [int(x) for x in geo.internal_surface_ids]
    ext_set = set(ext_ids)
    int_set = set(int_ids)
    all_set = ext_set | int_set

    PARTICLES = ("neutron", "photon")
    eps       = 1e-15
    chunk_set = set(int(c) for c in cell_ids)
    idx_map   = {int(c): i for i, c in enumerate(cell_ids)}

    def _cell_surfs(cid: int) -> List[dict]:
        x = info_.get(int(cid), {})
        return list(x.get("all_surfaces", [])) if isinstance(x, dict) else list(x or [])

    def _read_current_map(
        tally,
    ) -> Dict[Tuple[int, str], Tuple[float, float]]:
        df = tally.get_pandas_dataframe()
        return {
            (int(r["surface"]), str(r["particle"]).lower()):
            (float(r["mean"]), float(r["std. dev."]))
            for _, r in df.iterrows()
        }

    surf_to_cells: Dict[int, Set[int]] = {}
    for cid in cell_ids:
        for s in _cell_surfs(int(cid)):
            try:
                surf_to_cells.setdefault(
                    int(s["surface_id"]), set()
                ).add(int(cid))
            except Exception:
                pass

    Jnet_map: Dict[Tuple[int, str], Tuple[float, float]] = {}
    try:
        Jnet_map = _read_current_map(sp.get_tally(id=t_current_tally.id))
    except Exception as e:
        print(f"[warn] {chunk_key}: can't read total current tally ({e})")

    partial_map: Dict[Tuple[int, int, str], Tuple[float, float]] = {}
    for cid, tt in p_current_tallies.items():
        if int(cid) not in chunk_set:
            continue
        try:
            d = _read_current_map(sp.get_tally(id=tt.id))
            for (sid, part), v in d.items():
                partial_map[(int(cid), sid, part)] = v
        except Exception:
            pass

    info_flat = {int(c): _cell_surfs(int(c)) for c in cell_ids}
    special, special_ids, special_owner = (
        _build_special_surfaces(
            info_flat,
            armor_id=int(cell_ids[0]),
            breeder_back_id=int(cell_ids[-2]),
            vv_id=int(cell_ids[-1]),
            allowed=all_set,
            ext_set=ext_set,
            int_set=int_set,
        )
        if len(cell_ids) >= 3
        else ({}, set(), {})
    )
    print(f"[debug] {chunk_key} specials: {special}")

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

    for particle in PARTICLES:
        rows: List[dict] = []

        # ── External surfaces ─────────────────────────────────────────────────
        for sid in ext_ids:
            cset = surf_to_cells.get(sid, set()) & chunk_set
            if not cset:
                continue
            cid = sorted(cset)[0]
            if (sid, particle) not in Jnet_map:
                continue
            Jnet_m, Jnet_s = Jnet_map[(sid, particle)]
            Jout_m, Jout_s = partial_map.get((cid, sid, particle), (0.0, 0.0))
            Jnet_a = abs(Jnet_m); Jout_a = abs(Jout_m)

            if Jout_a < eps:
                A_mean = A_std = 0.0
            else:
                A_mean = abs(Jout_a - Jnet_a) / Jout_a
                A_var  = (
                    (Jnet_a / Jout_a**2)**2 * Jout_s**2
                    + (1.0 / Jout_a)**2    * Jnet_s**2
                )
                A_std = math.sqrt(max(A_var, 0.0))

            rows.append(dict(
                particle=particle, surface_id=sid, cell_id=cid, kind="external",
                J_total_mean=Jnet_m, J_total_std=Jnet_s,
                J_out_mean=Jout_m,   J_out_std=Jout_s,
                albedo_mean=A_mean,  albedo_std=A_std,
            ))

            if special.get("Armor_front_ext") == sid:
                J_in     = abs(Jout_a - Jnet_a)
                J_in_std = math.sqrt(Jout_s**2 + Jnet_s**2)
                rec = [{
                    "particle":    particle,
                    "special_key": "Armor_front_ext",
                    "surface_id":  sid,
                    "cell_id":     cid,
                    "J_out_mean":  Jout_a,  "J_out_std": abs(Jout_s),
                    "J_net_mean":  Jnet_m,  "J_net_std": abs(Jnet_s),
                    "J_in_mean":   J_in,    "J_in_std":  J_in_std,
                    "chunk_key":   chunk_key,
                }]
                with open(outdir / f"armor_current_{particle}.json", "w") as f:
                    json.dump(rec, f, indent=2)

        # ── Internal surfaces ─────────────────────────────────────────────────
        for sid in int_ids:
            cset = surf_to_cells.get(sid, set()) & chunk_set
            if len(cset) != 2:
                continue
            c1, c2   = sorted(cset)
            J1m, J1s = partial_map.get((c1, sid, particle), (0.0, 0.0))
            J2m, J2s = partial_map.get((c2, sid, particle), (0.0, 0.0))
            J1, J2   = abs(J1m), abs(J2m)

            if J1 <= eps or J2 <= eps:
                A1 = A2 = A1s = A2s = float("nan")
            else:
                A1  = J2 / J1
                A1s = math.sqrt((J2s/J1)**2 + (J2*J1s/J1**2)**2)
                A2  = J1 / J2
                A2s = math.sqrt((J1s/J2)**2 + (J1*J2s/J2**2)**2)

            for c, Am, As, Jom, Jos in [
                (c1, A1, A1s, J1m, J1s),
                (c2, A2, A2s, J2m, J2s),
            ]:
                rows.append(dict(
                    particle=particle, surface_id=sid, cell_id=c, kind="internal",
                    J_total_mean=float("nan"), J_total_std=float("nan"),
                    J_out_mean=float(Jom),     J_out_std=float(Jos),
                    albedo_mean=Am if np.isfinite(Am) else float("nan"),
                    albedo_std=As  if np.isfinite(As) else float("nan"),
                ))

        if not rows:
            continue

        df_alb = pd.DataFrame(rows)
        df_alb.to_csv(
            outdir / f"surface_albedo_all_{particle}_{chunk_key}.csv", index=False
        )

        special_rows = []
        for key, sid in special.items():
            want_cid = int(special_owner.get(key, -1))
            cand     = df_alb[df_alb["surface_id"] == sid]
            if want_cid != -1:
                c2 = cand[cand["cell_id"] == want_cid]
                if not c2.empty:
                    cand = c2
            if cand.empty:
                continue
            r = cand.iloc[0].to_dict(); r["special_key"] = key
            special_rows.append(r)
        pd.DataFrame(special_rows).to_csv(
            outdir / f"surface_albedo_special_{particle}_{chunk_key}.csv",
            index=False,
        )

        df_ext = df_alb[
            (df_alb["kind"] == "external") &
            (~df_alb["surface_id"].isin(excluded)) &
            (~df_alb["surface_id"].isin(special_ids))
        ].copy()
        df_ext.to_csv(
            outdir / f"surface_albedo_external_filtered_{particle}_{chunk_key}.csv",
            index=False,
        )

        cell_rows = []
        for cid, grp in df_ext.groupby("cell_id"):
            cid = int(cid)
            cell_rows.append(dict(
                particle=particle,
                cell_id=cid,
                layer_label=(labels[idx_map[cid]]
                             if cid in idx_map else f"cell_{cid}"),
                albedo_mean=float(grp["albedo_mean"].mean()),
                albedo_std=float(
                    np.sqrt(np.nansum(grp["albedo_std"]**2)) / max(len(grp), 1)
                ),
                n_surfaces=len(grp),
            ))
        df_cell = (pd.DataFrame(cell_rows).sort_values("cell_id")
                   if cell_rows else pd.DataFrame())
        df_cell.to_csv(
            outdir / f"cell_albedo_summary_external_{particle}_{chunk_key}.csv",
            index=False,
        )

        xpos = np.arange(len(labels))
        y    = np.full(len(labels), np.nan)
        e    = np.full(len(labels), np.nan)
        for _, row in df_cell.iterrows():
            cid = int(row["cell_id"])
            if cid in idx_map:
                y[idx_map[cid]] = row["albedo_mean"]
                e[idx_map[cid]] = row["albedo_std"]

        fig, ax = plt.subplots(figsize=(11, 6))
        ax.errorbar(xpos, y, yerr=e, fmt="none", elinewidth=1, capsize=4)
        ax.scatter(xpos, y, marker="x", s=90, linewidths=2,
                   label=f"{particle.capitalize()} cell avg")
        for r in special_rows:
            key = str(r.get("special_key", ""))
            cid = int(r["cell_id"])
            if cid not in idx_map:
                continue
            x0, y0 = idx_map[cid], float(r["albedo_mean"])
            if not np.isfinite(y0):
                continue
            col = SPECIAL_COLORS.get(key, "black")
            ax.scatter([x0], [y0], marker="x", s=120, color=col, zorder=6,
                       label=SPECIAL_LABELS.get(key, key))
        ax.set_xticks(xpos)
        ax.set_xticklabels(labels, rotation=45, fontsize=11)
        ax.set_xlabel("Layer"); ax.set_ylabel("Albedo")
        ax.set_title(f"{particle.capitalize()} Albedo: {chunk_key}")
        ax.grid(axis="y", linestyle="--", alpha=0.6)
        ax.legend(); fig.tight_layout()
        fig.savefig(outdir / f"albedo_layers_{particle}_{chunk_key}.png", dpi=300)
        plt.close(fig)


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

    is_ob      = chunk_key.startswith("OB_")
    region_tag = "OB" if is_ob else "IB"
    labels     = layer_names_for_chunk(cell_ids, region_tag)
    scaling    = scaling_for_cells(cell_ids)
    base_df    = make_chunk_base_df(
        chunk_key=chunk_key, cell_ids=cell_ids,
        labels=labels, xcent=xcentroids, xedges=xedges,
    )
    write_struct_origin_csv(outdir, chunk_key, cell_ids, labels)

    # ── Flux spectrum ─────────────────────────────────────────────────────────
    t_flux_id     = require_tally_id(sp, flux_tally, "flux")
    t_flux        = sp.get_tally(id=t_flux_id)
    cell_bins     = get_cell_bins(t_flux)
    mean_all, std_all = get_flux_spectrum_arrays(t_flux)

    nf_chunk = mean_all[:, 0, :]; nf_std = std_all[:, 0, :]
    pf_chunk = mean_all[:, 1, :]; pf_std = std_all[:, 1, :]

    nf = subset_by_cells(nf_chunk, cell_bins, cell_ids)
    ns = subset_by_cells(nf_std,   cell_bins, cell_ids)
    pf = subset_by_cells(pf_chunk, cell_bins, cell_ids)
    ps = subset_by_cells(pf_std,   cell_bins, cell_ids)

    # Neutron spectrum plot
    fig, ax = plt.subplots()
    for i, cid in enumerate(cell_ids):
        fl = nf[i].flatten() * scaling[int(cid)] / unit_lethargy
        sl = ns[i].flatten() * scaling[int(cid)] / unit_lethargy
        ax.loglog(energies[:-1], fl,
                  label=f"{labels[i]} ({xcentroids[i]:.1f} cm)",
                  color=colors[i])
        ax.fill_between(energies[:-1],
                        np.maximum(fl - sl, 1e-99), fl + sl,
                        alpha=0.25, color=colors[i], lw=0)
    ax.set(xlim=[1e-3, 1e8], ylim=[1e6, 1e17],
           xlabel="Energy [eV]",
           ylabel="Neutron flux per unit lethargy [1/cm²/s]")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.legend(fontsize=8, ncol=2)
    fig.savefig(outdir / f"n_flux_spectrum_{chunk_key}.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Photon spectrum plot
    fig, ax = plt.subplots()
    for i, cid in enumerate(cell_ids):
        fl = pf[i].flatten() * scaling[int(cid)] / unit_lethargy
        sl = ps[i].flatten() * scaling[int(cid)] / unit_lethargy
        ax.loglog(energies[:-1], fl,
                  label=f"{labels[i]} ({xcentroids[i]:.1f} cm)",
                  color=colors[i])
        ax.fill_between(energies[:-1],
                        np.maximum(fl - sl, 1e-99), fl + sl,
                        alpha=0.25, color=colors[i], lw=0)
    ax.set(xlim=[1e3, 1e8], ylim=[1e6, 1e17],
           xlabel="Energy [eV]",
           ylabel="Photon flux per unit lethargy [1/cm²/s]")
    ax.grid(True, which="both"); ax.legend(fontsize=8, ncol=2)
    fig.savefig(outdir / f"p_flux_spectrum_{chunk_key}.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ── Spectrum CSVs ─────────────────────────────────────────────────────────
    for name, fl_arr, fl_std in [("neutron", nf, ns), ("photon", pf, ps)]:
        df = pd.DataFrame({"E_mid_eV": E_mid})
        for i, lab in enumerate(labels):
            s = scaling[int(cell_ids[i])]
            df[lab]          = fl_arr[i].flatten() * s / unit_lethargy
            df[f"{lab}_std"] = fl_std[i].flatten() * s / unit_lethargy
        df.to_csv(outdir / f"{name}_spectrum_{chunk_key}.csv", index=False)

    # ── Total flux ────────────────────────────────────────────────────────────
    def _total(arr, sarr, cids):
        mn  = np.array([arr[i].sum() * scaling[int(c)]
                        for i, c in enumerate(cids)])
        std = np.array([np.sqrt((sarr[i]**2).sum()) * scaling[int(c)]
                        for i, c in enumerate(cids)])
        return mn, std

    tn, tn_s = _total(nf, ns, cell_ids)
    tp, tp_s = _total(pf, ps, cell_ids)

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    for (y, s, lab) in [(tn, tn_s, "neutron"), (tp, tp_s, "photon")]:
        lo = np.maximum(y - s, 1e-30); hi = y + s
        ax.step(xedges, np.r_[y, y[-1]], where="post", label=lab)
        ax.fill_between(xedges, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]],
                        step="post", alpha=0.15)
    ax.set(xlabel="Radial Position [cm]", ylabel="Total Flux [1/cm²/s]")
    ax.grid(True, which="both", linestyle="--"); ax.legend()
    fig.savefig(outdir / f"flux_total_{chunk_key}.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)

    save_profile(base_df, outdir, quantity="flux_total_neutron",
                 mean=tn, std=tn_s, units="1/cm2/s")
    save_profile(base_df, outdir, quantity="flux_total_photon",
                 mean=tp, std=tp_s, units="1/cm2/s")

    # ── Heating ───────────────────────────────────────────────────────────────
    t_heat_id = require_tally_id(sp, heating_tally, "heating")
    t_heat    = sp.get_tally(id=t_heat_id)
    hbins     = get_cell_bins(t_heat)
    hm_all    = t_heat.get_values(value="mean").flatten()
    hs_all    = t_heat.get_values(value="std_dev").flatten()
    hm        = subset_by_cells(hm_all, hbins, cell_ids)
    hs        = subset_by_cells(hs_all, hbins, cell_ids)

    hw  = np.array([hm[i] * scaling[int(c)] * ev_to_joule
                    for i, c in enumerate(cell_ids)])
    hws = np.array([hs[i] * scaling[int(c)] * ev_to_joule
                    for i, c in enumerate(cell_ids)])

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    lo = np.maximum(hw - hws, 1e-30); hi = hw + hws
    ax.step(xedges, np.r_[hw, hw[-1]], where="post")
    ax.fill_between(xedges, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]],
                    step="post", alpha=0.3)
    ax.set(xlabel="Radial Position [cm]", ylabel="Heating [W/cm³]")
    ax.grid(True, which="both", linestyle="--")
    fig.savefig(outdir / f"heating_{chunk_key}.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    save_profile(base_df, outdir, quantity="heating",
                 mean=hw, std=hws, units="W/cm3")

    # ── H production (H1 + H2 + H3) ──────────────────────────────────────────
    h_vals = np.zeros(len(cell_ids)); h_vstd = np.zeros(len(cell_ids))
    for i, cid in enumerate(cell_ids):
        cid = int(cid)
        t   = dpa_gas_map.get(cid)
        if t is None:
            continue
        scores_set = {str(s) for s in (t.scores or [])}
        if GAS_SCORES_EXPLICIT.issubset(scores_set):
            h1_score, h2_score, h3_score = (
                "H1-production", "H2-production", "H3-production"
            )
        elif GAS_SCORES_REACTION.issubset(scores_set):
            h1_score, h2_score, h3_score = "(n,Xp)", "(n,Xd)", "(n,Xt)"
        else:
            continue
        f_struct = cell_struct_origin_frac.get(cid, {}) or {}
        nu       = list(t.nuclides or [])
        m1, s1   = corrected_sum(t, score=h1_score, nuclides=nu, f_struct=f_struct)
        m2, s2   = corrected_sum(t, score=h2_score, nuclides=nu, f_struct=f_struct)
        m3, s3   = corrected_sum(t, score=h3_score, nuclides=nu, f_struct=f_struct)
        m        = m1 + m2 + m3
        s        = math.sqrt(s1**2 + s2**2 + s3**2)
        denom    = float(cell_total_atoms_struct.get(cid, 0.0))
        if denom > 0:
            fac       = neutron_source_rate * s_in_y / denom * 1e6
            h_vals[i] = m * fac; h_vstd[i] = s * fac

    lo = np.maximum(h_vals - h_vstd, 1e-30); hi = h_vals + h_vstd
    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.step(xedges, np.r_[h_vals, h_vals[-1]], where="post")
    ax.fill_between(xedges, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]],
                    step="post", alpha=0.3)
    ax.set(xlabel="Radial Position [cm]", ylabel="H [appm/fpy]")
    ax.grid(True, which="both", linestyle="--")
    fig.savefig(outdir / f"h_{chunk_key}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    save_profile(base_df, outdir, quantity="H_appm_fpy_struct_origin",
                 mean=h_vals, std=h_vstd, units="appm/fpy")

    # ── He production (He3 + He4) ─────────────────────────────────────────────
    he_vals = np.zeros(len(cell_ids)); he_vstd = np.zeros(len(cell_ids))
    for i, cid in enumerate(cell_ids):
        cid = int(cid)
        t   = dpa_gas_map.get(cid)
        if t is None:
            continue
        scores_set = {str(s) for s in (t.scores or [])}
        if GAS_SCORES_EXPLICIT.issubset(scores_set):
            sc3, sc4 = "He3-production", "He4-production"
        elif GAS_SCORES_REACTION.issubset(scores_set):
            sc3, sc4 = "(n,X3He)", "(n,Xa)"
        else:
            continue
        f  = cell_struct_origin_frac.get(cid, {}) or {}
        nu = list(t.nuclides or [])
        m3, s3 = corrected_sum(t, score=sc3, nuclides=nu, f_struct=f)
        m4, s4 = corrected_sum(t, score=sc4, nuclides=nu, f_struct=f)
        m  = m3 + m4; s = math.sqrt(s3**2 + s4**2)
        denom = float(cell_total_atoms_struct.get(cid, 0.0))
        if denom > 0:
            fac        = neutron_source_rate * s_in_y / denom * 1e6
            he_vals[i] = m * fac; he_vstd[i] = s * fac

    lo = np.maximum(he_vals - he_vstd, 1e-30); hi = he_vals + he_vstd
    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.step(xedges, np.r_[he_vals, he_vals[-1]], where="post")
    ax.fill_between(xedges, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]],
                    step="post", alpha=0.3)
    ax.set(xlabel="Radial Position [cm]", ylabel="He [appm/fpy]")
    ax.grid(True, which="both", linestyle="--")
    fig.savefig(outdir / f"he_{chunk_key}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    save_profile(base_df, outdir, quantity="He_appm_fpy_struct_origin",
                 mean=he_vals, std=he_vstd, units="appm/fpy")

    # ── DPA/fpy ───────────────────────────────────────────────────────────────
    dpa_vals = np.zeros(len(cell_ids)); dpa_vstd = np.zeros(len(cell_ids))
    for i, cid in enumerate(cell_ids):
        cid = int(cid)
        t   = dpa_gas_map.get(cid)
        if t is None:
            continue
        denom = float(cell_total_atoms_struct.get(cid, 0.0))
        if denom <= 0:
            continue
        f  = cell_struct_origin_frac.get(cid, {}) or {}
        sm = sv = 0.0
        for nuc in list(t.nuclides or []):
            nuc = str(nuc)
            try:
                dm, ds = tally_mean_std(t, score="damage-energy", nuclide=nuc)
            except Exception:
                continue
            fs = float(f.get(nuc, 0.0))
            if fs == 0.0:
                continue
            dm *= fs; ds *= fs
            Ed  = float(geo.materials.Ed(element_from_nuclide(nuc)))
            sc  = 0.8 / (2.0 * Ed) * neutron_source_rate * s_in_y
            sm += sc * dm / denom
            sv += (sc * ds / denom) ** 2
        dpa_vals[i] = sm; dpa_vstd[i] = math.sqrt(sv)

    lo = np.maximum(dpa_vals - dpa_vstd, 1e-30); hi = dpa_vals + dpa_vstd
    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.step(xedges, np.r_[dpa_vals, dpa_vals[-1]], where="post")
    ax.fill_between(xedges, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]],
                    step="post", alpha=0.3)
    ax.set(xlabel="Radial Position [cm]", ylabel="NRT-dpa/fpy")
    ax.grid(True, which="both", linestyle="--")
    fig.savefig(outdir / f"dpa_{chunk_key}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    save_profile(base_df, outdir, quantity="dpa_fpy_struct_origin",
                 mean=dpa_vals, std=dpa_vstd, units="DPA/fpy")

    if do_albedo:
        compute_albedo_for_chunk(
            sp, cell_ids=cell_ids, labels=labels,
            outdir=outdir, chunk_key=chunk_key,
            excluded_surfaces=set(),
        )

    print(f"[done] {chunk_key} → {outdir}")


# ──────────────────────────────────────────────────────────────────────────────
# Poloidal layer sweep  (tokamak only)
# ──────────────────────────────────────────────────────────────────────────────
def _layer_cells_region1(layer_index: int) -> Tuple[List[int], List[str]]:
    keys_ob = [f"OB_1_b{b}" for b in range(1, n_breeder + 1)]
    keys_ib = [f"IB_1_b{b}" for b in range(1, n_breeder + 1)]
    ids: List[int] = []; labs: List[str] = []
    for k in keys_ob + keys_ib:
        chunk = ob_by_key.get(k) if k.startswith("OB_") else ib_by_key.get(k)
        if not chunk or layer_index >= len(chunk):
            continue
        ids.append(int(chunk[layer_index])); labs.append(k)
    return ids, labs


def _empty_poloidal_plot(x_data, num_pts):
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    split   = n_breeder - 0.5
    ax.set_xlim([-0.5, len(x_data) - 2.5])
    ax.axvline(split, linestyle="--", linewidth=1, color="k")
    ax.text(split - 0.15, 0.95, "← outboard   inboard →",
            transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=10)
    xticks = np.arange(num_pts)
    ax.set_xticks(xticks)
    ax.set_xticklabels(
        [f"OB{i+1}" if i < n_breeder else f"IB{i+1-n_breeder}"
         for i in range(num_pts)],
        rotation=45, ha="right", fontsize=9,
    )
    ax.set_xlabel("Poloidal regions", fontsize=11)
    ax.grid(True, which="major", linestyle="--", linewidth=0.6, alpha=0.8)
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
    outdir.mkdir(parents=True, exist_ok=True)
    cell_ids_layer, _ = _layer_cells_region1(layer_index)
    n = len(cell_ids_layer)
    if n == 0:
        print(f"[warn] layer_index={layer_index} ({layer_tag}): "
              f"no cells found, skipping.")
        return

    scaling = scaling_for_cells(cell_ids_layer)
    x       = _pad(np.arange(n, dtype=float))

    def _plot_step(y, ys, ylabel, fname, yscale="auto"):
        y = np.asarray(y, float); ys = np.asarray(ys, float)
        yp  = _pad(y);  ysp = _pad(ys)
        lo  = np.clip(yp - ysp, 1e-99, None); hi = yp + ysp
        fig, ax = _empty_poloidal_plot(x, n)
        set_ylim_and_ticks(ax, yp, scale=yscale)
        ax.set_yscale("log" if yscale == "log" else ax.get_yscale())
        ax.step(x, yp, where="mid", linewidth=2)
        ax.fill_between(x, lo, hi, step="mid", alpha=0.15)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(f"{ylabel} — {layer_tag}", fontsize=12)
        fig.savefig(outdir / fname, dpi=300, bbox_inches="tight")
        plt.close(fig)

    # ── Flux ──────────────────────────────────────────────────────────────────
    t_flux    = sp.get_tally(id=require_tally_id(sp, flux_tally, "flux"))
    cell_bins = get_cell_bins(t_flux)
    ma, sa    = get_flux_spectrum_arrays(t_flux)
    nfl = subset_by_cells(ma[:, 0, :], cell_bins, cell_ids_layer)
    nfs = subset_by_cells(sa[:, 0, :], cell_bins, cell_ids_layer)
    pfl = subset_by_cells(ma[:, 1, :], cell_bins, cell_ids_layer)
    pfs = subset_by_cells(sa[:, 1, :], cell_bins, cell_ids_layer)
    neut   = np.array([nfl[i].sum() * scaling[int(c)]
                       for i, c in enumerate(cell_ids_layer)])
    neut_s = np.array([np.sqrt((nfs[i]**2).sum()) * scaling[int(c)]
                       for i, c in enumerate(cell_ids_layer)])
    phot   = np.array([pfl[i].sum() * scaling[int(c)]
                       for i, c in enumerate(cell_ids_layer)])
    phot_s = np.array([np.sqrt((pfs[i]**2).sum()) * scaling[int(c)]
                       for i, c in enumerate(cell_ids_layer)])

    fig, ax = _empty_poloidal_plot(x, n)
    set_ylim_and_ticks(ax, np.r_[_pad(neut), _pad(phot)], scale="auto")
    ax.step(x, _pad(np.clip(neut, 1e-99, None)),
            where="mid", linewidth=2, label="neutron")
    ax.step(x, _pad(np.clip(phot, 1e-99, None)),
            where="mid", linewidth=2, label="photon")
    ax.fill_between(x, _pad(np.clip(neut - neut_s, 1e-99, None)),
                    _pad(neut + neut_s), step="mid", alpha=0.15)
    ax.fill_between(x, _pad(np.clip(phot - phot_s, 1e-99, None)),
                    _pad(phot + phot_s), step="mid", alpha=0.15)
    ax.set_ylabel("Total Flux [1/cm²/s]")
    ax.set_title(f"Total Flux — {layer_tag}")
    ax.legend()
    fig.savefig(outdir / f"flux_total_{layer_tag}.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ── Heating ───────────────────────────────────────────────────────────────
    t_heat   = sp.get_tally(id=require_tally_id(sp, heating_tally, "heating"))
    hbins    = get_cell_bins(t_heat)
    hm_all   = t_heat.get_values(value="mean").flatten()
    hs_all   = t_heat.get_values(value="std_dev").flatten()
    hm_layer = subset_by_cells(hm_all, hbins, cell_ids_layer)
    hs_layer = subset_by_cells(hs_all, hbins, cell_ids_layer)
    hw  = np.array([hm_layer[i] * scaling[int(c)] * ev_to_joule
                    for i, c in enumerate(cell_ids_layer)])
    hws = np.array([hs_layer[i] * scaling[int(c)] * ev_to_joule
                    for i, c in enumerate(cell_ids_layer)])
    _plot_step(hw, hws, "Heating [W/cm³]", f"heating_{layer_tag}.png")

    # ── DPA ───────────────────────────────────────────────────────────────────
    dpa_y = np.zeros(n); dpa_s = np.zeros(n)
    for i, cid in enumerate(cell_ids_layer):
        cid = int(cid)
        t   = dpa_gas_map.get(cid)
        if t is None:
            continue
        denom = float(cell_total_atoms_struct.get(cid, 0.0))
        if denom <= 0:
            continue
        f  = cell_struct_origin_frac.get(cid, {}) or {}
        sm = sv = 0.0
        for nuc in list(t.nuclides or []):
            nuc = str(nuc)
            try:
                dm, ds = tally_mean_std(t, score="damage-energy", nuclide=nuc)
            except Exception:
                continue
            fs = float(f.get(nuc, 0.0))
            if fs == 0.0:
                continue
            dm *= fs; ds *= fs
            Ed  = float(geo.materials.Ed(element_from_nuclide(nuc)))
            sc  = 0.8 / (2.0 * Ed) * neutron_source_rate * s_in_y
            sm += sc * dm / denom; sv += (sc * ds / denom) ** 2
        dpa_y[i] = sm; dpa_s[i] = math.sqrt(max(sv, 0.0))
    _plot_step(dpa_y, dpa_s, "NRT-dpa/fpy", f"dpa_{layer_tag}.png",
               yscale="linear")

    # ── H production (H1 + H2 + H3) ──────────────────────────────────────────
    h_vals = np.zeros(n); h_vstd = np.zeros(n)
    for i, cid in enumerate(cell_ids_layer):
        cid = int(cid); t = dpa_gas_map.get(cid)
        if t is None:
            continue
        sc_set = {str(s) for s in (t.scores or [])}
        if GAS_SCORES_EXPLICIT.issubset(sc_set):
            h1_score, h2_score, h3_score = (
                "H1-production", "H2-production", "H3-production"
            )
        elif GAS_SCORES_REACTION.issubset(sc_set):
            h1_score, h2_score, h3_score = "(n,Xp)", "(n,Xd)", "(n,Xt)"
        else:
            continue
        f  = cell_struct_origin_frac.get(cid, {}) or {}
        nu = list(t.nuclides or [])
        m1, s1 = corrected_sum(t, score=h1_score, nuclides=nu, f_struct=f)
        m2, s2 = corrected_sum(t, score=h2_score, nuclides=nu, f_struct=f)
        m3, s3 = corrected_sum(t, score=h3_score, nuclides=nu, f_struct=f)
        m      = m1 + m2 + m3
        s      = math.sqrt(s1**2 + s2**2 + s3**2)
        denom  = float(cell_total_atoms_struct.get(cid, 0.0))
        if denom > 0:
            fac       = neutron_source_rate * s_in_y / denom * 1e6
            h_vals[i] = m * fac; h_vstd[i] = s * fac
    _plot_step(h_vals, h_vstd,
               "H production [appm/fpy]", f"h1_{layer_tag}.png")

    # ── He production ─────────────────────────────────────────────────────────
    he_vals = np.zeros(n); he_vstd = np.zeros(n)
    for i, cid in enumerate(cell_ids_layer):
        cid = int(cid); t = dpa_gas_map.get(cid)
        if t is None:
            continue
        sc_set = {str(s) for s in (t.scores or [])}
        if GAS_SCORES_EXPLICIT.issubset(sc_set):
            sc3, sc4 = "He3-production", "He4-production"
        elif GAS_SCORES_REACTION.issubset(sc_set):
            sc3, sc4 = "(n,X3He)", "(n,Xa)"
        else:
            continue
        f  = cell_struct_origin_frac.get(cid, {}) or {}
        nu = list(t.nuclides or [])
        m3, s3 = corrected_sum(t, score=sc3, nuclides=nu, f_struct=f)
        m4, s4 = corrected_sum(t, score=sc4, nuclides=nu, f_struct=f)
        m = m3 + m4; s = math.sqrt(s3**2 + s4**2)
        denom = float(cell_total_atoms_struct.get(cid, 0.0))
        if denom > 0:
            fac         = neutron_source_rate * s_in_y / denom * 1e6
            he_vals[i]  = m * fac; he_vstd[i] = s * fac
    _plot_step(he_vals, he_vstd,
               "He production [appm/fpy]", f"he_{layer_tag}.png")

    print(f"[done] layer {layer_tag} → {outdir}")


# ──────────────────────────────────────────────────────────────────────────────
# MAIN RUN
# ──────────────────────────────────────────────────────────────────────────────
with openmc.StatePoint(str(STATEPOINT_FILE)) as sp:
    dpa_gas_map = build_dpa_gas_map(sp)

    if cfg.SIM_TYPE == "slab":
        keys_to_process = [cfg.ALBEDO_CHUNK_KEY]
    else:
        keys_to_process = cfg.NEUTRONICS_KEYS_TO_PROCESS

    do_albedo = (cfg.DO_ALBEDO and cfg.SIM_TYPE == "tokamak")

    for key in keys_to_process:
        chunk_cells = ob_by_key.get(key) or ib_by_key.get(key)
        if chunk_cells is None:
            print(f"[warn] Unknown chunk key: {key} – skipping")
            continue
        try:
            xcent, _, xedges = radial_bins_for_key(key)
        except Exception as e:
            print(f"[warn] radial bins for {key}: {e} – skipping")
            continue
        process_chunk(
            sp, key, chunk_cells,
            dpa_gas_map=dpa_gas_map,
            xcentroids=xcent,
            xedges=xedges,
            do_albedo=do_albedo,
        )

    if cfg.SIM_TYPE == "tokamak":
        for layer_tag, layer_index in cfg.POLOIDAL_LAYERS:
            layer_outdir = RESULTS_DIR / f"layer_{layer_tag}_region1_only"
            process_layer_region1_only(
                sp, layer_index, layer_tag, layer_outdir,
                dpa_gas_map=dpa_gas_map,
            )
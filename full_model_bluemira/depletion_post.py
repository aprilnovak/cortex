#!/usr/bin/env python3
"""
depletion_post.py
=================
Post-process OpenMC depletion results (activity + decay heat) for both
tokamak and slab simulations.

Reads
-----
  cfg.DEPLETION_RUN_DIR/r2s/activation/depletion_results.h5
  cfg.DEPLETION_RUN_DIR/r2s/activation/cell_material_map.csv

Outputs per chunk key → cfg.DEPLETION_RESULTS_DIR/<chunk_key>/
  activity_cell_<region>.png
  activity_all_cells.png / .csv
  radial_activity_profiles.png
  decayheat_<region>.png
  decayheat_all_cells.png / .csv
  radial_decayheat_profiles.png
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import cycle
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import openmc
import openmc.deplete

import inputs as cfg

# ──────────────────────────────────────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────────────────────────────────────
R2S_ACTIV_DIR         = cfg.DEPLETION_RUN_DIR / "r2s" / "activation"
DEPLETION_RESULTS_H5  = R2S_ACTIV_DIR / "depletion_results.h5"
CELL_MATERIAL_MAP_CSV = R2S_ACTIV_DIR / "cell_material_map.csv"
RESULTS_OUT_DIR       = cfg.DEPLETION_RESULTS_DIR
RESULTS_OUT_DIR.mkdir(parents=True, exist_ok=True)

SECONDS_PER_YEAR = 365.0 * 24.0 * 3600.0
MAX_HALFLIFE = 60*60*24*365/12*1e6*100
# ──────────────────────────────────────────────────────────────────────────────
# Colour / style helpers
# ──────────────────────────────────────────────────────────────────────────────
def generate_colors(n):
    """Generates a smooth rainbow gradient of n RGB colors."""
    cmap = plt.get_cmap('turbo') # plasma also looks nice
    color_range = cmap(np.linspace(1, 0, n))
    return color_range

_MARKERS    = ["o", "s", "^", "D", "v", "P", "X", "<", ">", "h"]
_LINESTYLES = ["-"]

def _style_cycle(reset: bool = True):
    """Return new independent cycles so each plot starts consistent."""
    m = cycle(_MARKERS)
    ls = cycle(_LINESTYLES)
    return m, ls

def _finalize_legend(ax, *, ncol: int = 2, fontsize: int = 8, title: str = ""):
    """Sorted, multi-column, outside legend."""
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return

    order = np.argsort(labels)
    handles = [handles[i] for i in order]
    labels = [labels[i] for i in order]

    ax.legend(
        handles,
        labels,
        title=title,
        fontsize=fontsize,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        ncol=ncol,
        frameon=False,
        borderaxespad=0.0,
        handlelength=2.0,
        columnspacing=1.0,
        labelspacing=0.4,
    )

def _format_log_axes(ax):
    """Cleaner log grid and ticks."""
    ax.grid(True, which="major", linewidth=0.8, alpha=0.35)
    ax.grid(True, which="minor", linewidth=0.5, alpha=0.15)

    ax.xaxis.set_minor_locator(mticker.LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_locator(mticker.LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())

def _has_positive_finite(y) -> bool:
    y = np.asarray(y, float)
    return bool(np.any(np.isfinite(y) & (y > 0.0)))

def _ymin_from_topn_edge(global_min_topn_edge: float) -> Optional[float]:
    try:
        v = float(global_min_topn_edge)
    except Exception:
        return None
    if not np.isfinite(v) or v <= 0.0:
        return None
    return 0.1 * v

def _add_time_reference_lines(ax):
    """Useful reference lines (x-axis in years)."""

    _S  = 1.0
    _M  = 60.0
    _H  = 3_600.0
    _D  = 86_400.0
    _W  = 7.0 * 86_400.0           
    _MO = SECONDS_PER_YEAR / 12.0

    refs = [
        (_S  / SECONDS_PER_YEAR, "1 s"),
        (_M  / SECONDS_PER_YEAR, "1 m"),
        (_H  / SECONDS_PER_YEAR, "1 h"),
        (_D  / SECONDS_PER_YEAR, "1 d"),   
        (_W  / SECONDS_PER_YEAR, "1 w"),
        (_MO / SECONDS_PER_YEAR, "1 mo"),
        (1.0,                    "1 y"),
        (10.0,                   "10 y"),
        (100.0,                  "100 y"),
        (1_000.0,                "1000 y"),
    ]

    for x, txt in refs:
        ax.axvline(x, color="gray", linestyle="--", linewidth=1, alpha=0.6)
        ax.text(
            x,
            0.98,
            txt,
            transform=ax.get_xaxis_transform(),
            rotation=90,
            va="top",
            ha="right",
            fontsize=8,
            alpha=0.8,
        )

def _safe_half_life(nuclide: str) -> Optional[float]:
    """
    Return numeric half-life in seconds when available, else None.
    """
    try:
        hl = openmc.data.half_life(nuclide)
    except Exception:
        return None

    if hl is None:
        return None

    try:
        hl = float(hl)
    except Exception:
        return None

    if not np.isfinite(hl):
        return None

    return hl

def _sortable_half_life(nuclide: str) -> float:
    """
    Sort key for nuclides by half-life.
    Stable/unknown nuclides are pushed to the end.
    """
    hl = _safe_half_life(nuclide)
    return hl if hl is not None else np.inf

def display_half_life(nuclide):
  half_life = openmc.data.half_life(nuclide)

  if half_life is None:
    return " (stable/unknown)"

  # display in best units for s, h, d, y; this guarantees that each of these will be at minimum
  # 0.01 in their base units
  if (half_life < 60): # less than 1 min, display in s
    return ' ({0:.2f} s)'.format(half_life)
  elif (half_life < 60*60): # less than 1 hour, display in m
    return ' ({0:.2f} m)'.format(half_life/(60))
  elif (half_life < 60*60*24): # less than 1 day, display in h
    return ' ({0:.2f} h)'.format(half_life / (60*60))
  elif (half_life < 60*60*24*365/12.): # less than 1 month, display in d
    return ' ({0:.2f} d)'.format(half_life / (60*60*24))
  elif (half_life < 60*60*24*365/12*100000): # less than 100000 year, display in y
    return ' ({0:.2f} y)'.format(half_life / (60*60*24*365))
  elif (half_life < 60*60*24*365/12*1e6*100): # less than 100 My year, display in My
    return ' ({0:.2f} My)'.format(half_life / (60*60*24*365*1e6))
  else:
    return ' ( > 100 My)'
    #return ' ({0:.2f} Gy)'.format(half_life / (60*60*24*365/12.*1e9))

# ──────────────────────────────────────────────────────────────────────────────
# Nuclide selection
# ──────────────────────────────────────────────────────────────────────────────
def _select_topn(
    sorted_items: list[tuple[str, float]]
    ) -> list[tuple[str, float]]:
    """
    For a given step, we add nuclides to plot until the next incremental percent of
    total activity to show is below 5% of the total. This limits the number of nuclides
    plotted and is more robust than simply plotting the top 95% of nuclides because in
    short cooling times, there's dozens of nuclides which contribute to the total activity,
    which would result in way too many nuclides to plot.
    """

    total = 0.0
    for i in range(len(sorted_items)):
      total += float(sorted_items[i][1])

    # No induced activity 
    if total == 0.0:
        print("zero induced activity")
        return []

    running_total = 0.0
    prev_running_total = 0.0
    index = len(sorted_items) - 1
    for i in range(len(sorted_items)):
      prev_running_total = running_total
      running_total += float(sorted_items[i][1])

      # could have just a single very strong nuclide
      if (i == 0 and running_total / total >= 0.99):
        index = i
        break

      # otherwise, compare the gain
      if ((running_total - prev_running_total) / total <= 0.05):
        index = i
        break

    return sorted_items[:index + 1]

def _select(dict_list, n_steps):
    # selections per step; first, loop through all the time steps to find the top nuclides
    # on each given step. Then, we take the union of these and then obtain the data to
    # plot on each step by writing that nuclide for all time steps
    topn_list_by_step: list[list[tuple[str, float]]] = [[] for _ in range(n_steps)]
    top_nucs_union: set[str] = set()

    for istep in range(n_steps):
        d = dict_list[istep] or {}
        if not d:
            continue

        step_sorted = sorted(d.items(), key=lambda x: x[1], reverse=True)
        for nuc, _ in _select_topn(step_sorted):
            top_nucs_union.add(nuc)

    for istep in range(n_steps):
        d = dict_list[istep] or {}
        if not d:
            continue

        step_sorted = sorted(d.items(), key=lambda x: x[1], reverse=True)
        selected = []
        for s in step_sorted:
          if (s[0] in top_nucs_union):
            selected.append(s)

        topn_list_by_step[istep] = selected

    return topn_list_by_step, top_nucs_union

# ──────────────────────────────────────────────────────────────────────────────
# Depletion mapping
# ──────────────────────────────────────────────────────────────────────────────
@dataclass(frozen=True)
class DepletionMapping:
    cell_to_mat: Dict[int, str]            # cell_id -> mat_id (string)
    cells_by_mat: Dict[str, List[int]]     # mat_id -> [cell_id, ...]
    mat_id_to_name: Dict[str, str]         # mat_id -> name
    vol_by_mat: Dict[str, float]           # mat_id -> volume (from results[0].volume)

def load_depletion_mapping(
    results: openmc.deplete.Results,
    map_csv: Path = CELL_MATERIAL_MAP_CSV,
) -> DepletionMapping:
    """
    Loads depletion_run/r2s/activation/cell_material_map.csv and volumes from step0.
    """
    if not map_csv.is_file():
        raise FileNotFoundError(f"Cell-material map CSV not found: {map_csv}")

    map_df = pd.read_csv(str(map_csv))
    map_df["mat_id"] = map_df["mat_id"].astype(str)

    return DepletionMapping(
        cell_to_mat   = dict(zip(map_df["cell_id"], map_df["mat_id"])),
        cells_by_mat  = map_df.groupby("mat_id")["cell_id"].apply(list).to_dict(),
        mat_id_to_name= map_df.groupby("mat_id")["material_name"].first().to_dict(),
        vol_by_mat    = results[0].volume,
    )

# ──────────────────────────────────────────────────────────────────────────────
# Layer tags
# ──────────────────────────────────────────────────────────────────────────────
def make_default_layer_tags(
    n_breeder_layers: int,
    include_armor: bool = True,
    include_fw: bool = True,
    include_vv: bool = True,
    breeder_prefix: str = "OB",
    ) -> List[str]:
    """
    Build tags like:
      Armor, First Wall, OB_1..OB_N, VV
    where N is provided at runtime.

    Assumes per-chunk ordering:
      [Armor, FW, breeder_layers..., VV]
    """
    tags: List[str] = []
    if include_armor:
        tags.append("Armor")
    if include_fw:
        tags.append("First wall")

    tags.extend([f"Breeder layer {i}" for i in range(1, n_breeder_layers + 1)])

    if include_vv:
        tags.append("VV")
    return tags

def build_cell_id_to_name(
    cell_ids: Sequence[int],
    layer_tags: Optional[Sequence[str]] = None,
    chunk_key: Optional[str] = None,
) -> Dict[int, str]:
    """
    Map each cell_id to human-friendly region name.
    If tags mismatch, fall back to Layer_00-style names.
    """
    if layer_tags and len(layer_tags) == len(cell_ids):
        return {cid: str(layer_tags[i]) for i, cid in enumerate(cell_ids)}
    
    prefix = f"{chunk_key}_" if chunk_key else ""
    return {cid: f"{prefix}Layer_{i:02d}" for i, cid in enumerate(cell_ids)}

# ──────────────────────────────────────────────────────────────────────────────
# CSV helpers
# ──────────────────────────────────────────────────────────────────────────────
def _safe_csv_label(text: str) -> str:
    return str(text).strip().replace(" ", "_").replace("/", "_")

def save_total_timeseries_csv(
    total_by_cell: Dict[int, Tuple[np.ndarray, np.ndarray]],
    cell_id_to_name: Dict[int, str],
    out_csv: Path,
    value_label: str,
) -> None:
    """
    Save total quantity over cooling time for all cells into one CSV.

    Output columns:
      Cooling_time_years,
      <RegionName>__cell_<cid>, ...
    """
    if not total_by_cell:
        return
    
    # Build one dataframe per cell, then outer-merge on time
    dfs = []
    for cid, (t_rel_y, values) in total_by_cell.items():
        col = f"{_safe_csv_label(cell_id_to_name.get(cid, str(cid)))}__cell_{cid}"
        dfs.append(pd.DataFrame({"Cooling_time_years": t_rel_y, col: values}))
    df_out = dfs[0]
    for df_i in dfs[1:]:
        df_out = df_out.merge(df_i, on="Cooling_time_years", how="outer")
    df_out.sort_values("Cooling_time_years").to_csv(out_csv, index=False)

# ──────────────────────────────────────────────────────────────────────────────
# Shutdown index
# ──────────────────────────────────────────────────────────────────────────────
def shutdown_index_from_source_rates(source_rates: np.ndarray, n_steps: int) -> int:
    """
    source_rates length = number of intervals.
    results times length = n_steps.
    End-of-irradiation 'state' index = last non-zero interval + 1.
    """
    irr_intervals = np.where(source_rates > 0.0)[0]
    if irr_intervals.size == 0:
        return 0
    idx = int(irr_intervals[-1] + 1)
    return min(idx, n_steps - 1)

# ──────────────────────────────────────────────────────────────────────────────
# Activity
# ──────────────────────────────────────────────────────────────────────────────
def plot_activity_nuclides_per_cell(
    results: openmc.deplete.Results,
    mapping: DepletionMapping,
    cell_ids: Sequence[int],
    cell_id_to_name: Dict[int, str],
    out_dir: Path,
    *,
    idx_shutdown: int,
    activity_units: str = "Bq/kg",
) -> Tuple[Dict[int, Tuple[np.ndarray, np.ndarray]], float]:
    """
    Returns:
      total_activity_all[cid] = (t_rel_years, total_activity_in_units)
    """

    out_dir.mkdir(parents=True, exist_ok=True)
    total_all:  Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    global_min_topn_edge = np.inf

    for cid in cell_ids:
        if cid not in mapping.cell_to_mat:
            print(f"[Activity] cell {cid} not in map; skipping")
            continue

        mat_id = str(mapping.cell_to_mat[cid])
        if mat_id not in mapping.vol_by_mat:
            print(f"[Activity] mat {mat_id} missing volume; skipping cell {cid}")
            continue
        vol = mapping.vol_by_mat[mat_id]

        _, total_act   = results.get_activity(mat=mat_id, by_nuclide=False,
                                               units=activity_units, volume=vol)
        time_nuc, act_list = results.get_activity(mat=mat_id, by_nuclide=True,
                                                   units=activity_units, volume=vol)

        time_grid  = np.asarray(time_nuc, float)
        n_steps    = len(time_grid)

        if not (0 <= int(idx_shutdown) < n_steps):
            raise ValueError(
                f"idx_shutdown={idx_shutdown} out of bounds for n_steps={n_steps}"
            )
        if len(total_act) != n_steps:
            raise ValueError(
                f"Time grid mismatch: len(total_act)={len(total_act)} vs n_steps={n_steps} "
                f"for cell {cid}, mat {mat_id}."
            )

        total_act  = np.asarray(total_act, float)
        t_rel      = time_grid - time_grid[int(idx_shutdown)]
        mask       = t_rel > 0.0
        t_plot     = (t_rel[mask] / SECONDS_PER_YEAR).astype(float)
        act_plot   = total_act[mask]
        total_all[cid] = (t_plot, act_plot)

        region = cell_id_to_name.get(cid, str(cid))
        if not _has_positive_finite(act_plot):
            print(f"[Activity] Skipping plot for cell {cid} ({region}): zero activity.")
            continue

        topn_by_step, top_nucs_union = _select(act_list, n_steps)
        top_nucs = sorted(top_nucs_union)

        nuc_series    = {nuc: np.full(n_steps, np.nan) for nuc in top_nucs}
        others_series = np.full(n_steps, np.nan)

        for istep in range(n_steps):
            d = act_list[istep] or {}
            if not d:
                continue
            top_sum = 0.0
            for nuc, val in topn_by_step[istep]:
                nuc_series[nuc][istep] = float(val)
                top_sum += float(val)
            others_series[istep] = max(0.0, float(total_act[istep]) - top_sum)

        def _m(arr): return np.asarray(arr, float)[mask]

        others_plot = _m(others_series)
        hlists      = sorted((_sortable_half_life(n), n) for n in top_nucs)
        sorted_nucs = [n for _, n in hlists]
        colors      = generate_colors(len(top_nucs))

        fig, ax = plt.subplots(figsize=(11, 6))
        mc, lsc = _style_cycle()

        # y-axis limits
        others_pos = others_plot[others_plot > 0.0]
        total_pos  = act_plot[act_plot > 0.0]
        if others_pos.size > 0 and total_pos.size > 0:
            ax.set_ylim(
                10 ** math.floor(math.log10(np.min(others_pos))),
                10 ** math.ceil( math.log10(np.max(total_pos))),
            )

        for nuc in top_nucs:
            vp = _m(nuc_series[nuc])
            if not _has_positive_finite(vp):
                continue
            ax.loglog(t_plot, vp,
                      label=nuc + display_half_life(nuc),
                      color=colors[sorted_nucs.index(nuc) % len(colors)],
                      marker=next(mc), linestyle=next(lsc),
                      markersize=3, linewidth=1.3)

        if _has_positive_finite(others_plot):
            ax.loglog(t_plot, others_plot, label="Others",
                      linewidth=2, color="black", linestyle="--")
        if _has_positive_finite(act_plot):
            ax.loglog(t_plot, act_plot, color="black", linewidth=2,
                      label="Total", zorder=10)

        ax.set(xlabel="Time after irradiation [years]",
               ylabel=f"Activity [{activity_units}]",
               title=f"Activity — {region}")
        _add_time_reference_lines(ax); _format_log_axes(ax)
        _finalize_legend(ax, ncol=2 if len(top_nucs) <= 24 else 3,
                         title="Components")
        fig.subplots_adjust(right=0.72)
        fig.savefig(out_dir / f"activity_cell_{_safe_csv_label(region)}.png",
                    dpi=200, bbox_inches="tight")
        plt.close(fig)

    return total_all, float(global_min_topn_edge) if np.isfinite(global_min_topn_edge) else 0.0


def plot_activity_all_cells(
    total_all: Dict[int, Tuple[np.ndarray, np.ndarray]],
    cell_id_to_name: Dict[int, str],
    out_dir: Path,
    *,
    activity_units: str = "Bq/kg",
    global_min_topn_edge: float = 0.0,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(11, 6))
    mc, lsc = _style_cycle()
    colors  = plt.cm.tab20(np.linspace(0, 1, max(1, min(len(total_all), 20))))
    ymax    = 0.0

    for i, (cid, (t, y)) in enumerate(total_all.items()):
        yy = np.asarray(y, float)
        if not _has_positive_finite(yy):
            continue
        ymax = max(ymax, float(np.nanmax(yy)))
        ax.loglog(t, yy, label=cell_id_to_name.get(cid, str(cid)),
                  color=colors[i % len(colors)],
                  marker=next(mc), linestyle=next(lsc),
                  markersize=3.5, linewidth=1.1)

    ax.set(xlabel="Time after irradiation [years]",
           ylabel=f"Activity [{activity_units}]",
           title="Total Activity over cooling time")
    _add_time_reference_lines(ax); _format_log_axes(ax)
    ymin = _ymin_from_topn_edge(global_min_topn_edge)
    if ymin and ymax > 0:
        ax.set_ylim(ymin, 10 * ymax)
    _finalize_legend(ax, title="Regions")
    fig.subplots_adjust(right=0.72)
    fig.savefig(out_dir / "activity_all_cells.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_activity_radial(
    total_all: Dict[int, Tuple[np.ndarray, np.ndarray]],
    cell_ids: Sequence[int],
    xcentroids: Sequence[float],
    out_dir: Path,
    *,
    idx_to_plot: Sequence[int] = (0, 4, 8, 12, 16, 20),
    activity_units: str = "Bq/kg",
    global_min_topn_edge: float = 0.0,
) -> None:
    ref_cid = next((c for c in cell_ids if c in total_all), None)
    if ref_cid is None:
        return
    t_rel_y, _ = total_all[ref_cid]
    n_times     = len(t_rel_y)
    idx_to_plot = [int(i) for i in idx_to_plot if int(i) < n_times]
    if not idx_to_plot:
        return

    x      = np.asarray(xcentroids, float)
    radial = np.full((len(cell_ids), n_times), np.nan)
    for i, cid in enumerate(cell_ids):
        if cid in total_all:
            radial[i, :] = np.asarray(total_all[cid][1], float)

    ymax = float(np.nanmax(radial)) if np.any(
        np.isfinite(radial) & (radial > 0)
    ) else 1.0

    out_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 6))
    mc, lsc = _style_cycle()
    colors  = plt.cm.tab20(np.linspace(0, 1, max(1, len(idx_to_plot))))

    for k, j in enumerate(idx_to_plot):
        y_safe = np.where(radial[:, j] > 0, radial[:, j], np.nan)
        ax.semilogy(x, y_safe, label=f"t = {t_rel_y[j]:.2e} y",
                    color=colors[k % len(colors)],
                    marker=next(mc), linestyle=next(lsc),
                    markersize=4, linewidth=1.5)

    ax.set(xlabel="Radial position [cm]",
           ylabel=f"Total activity [{activity_units}]",
           title="Radial activity profiles")
    ax.grid(True, which="both", alpha=0.3)
    ymin = _ymin_from_topn_edge(global_min_topn_edge)
    if ymin and ymax > 0:
        ax.set_ylim(ymin, 10 * ymax)
    _finalize_legend(ax, ncol=1 if len(idx_to_plot) <= 12 else 2,
                     title="Timesteps")
    fig.subplots_adjust(right=0.76)
    fig.savefig(out_dir / "radial_activity_profiles.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────────
# Decay heat
# ──────────────────────────────────────────────────────────────────────────────
def plot_decayheat_nuclides_per_cell(
    results: openmc.deplete.Results,
    mapping: DepletionMapping,
    cell_ids: Sequence[int],
    cell_id_to_name: Dict[int, str],
    out_dir: Path,
    *,
    idx_shutdown: int,
    decayheat_units: str = "W/cm3",
) -> Tuple[Dict[int, Tuple[np.ndarray, np.ndarray]], float, float]:
    out_dir.mkdir(parents=True, exist_ok=True)
    total_all:  Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    max_heat   = 0.0
    global_min_topn_edge = np.inf

    for cid in cell_ids:
        if cid not in mapping.cell_to_mat:
            continue
        mat_id = str(mapping.cell_to_mat[cid])
        if mat_id not in mapping.vol_by_mat:
            continue
        vol = mapping.vol_by_mat[mat_id]

        _, total_h         = results.get_decay_heat(mat=mat_id, by_nuclide=False,
                                                     units=decayheat_units, volume=vol)
        time_nuc, heat_list = results.get_decay_heat(mat=mat_id, by_nuclide=True,
                                                      units=decayheat_units, volume=vol)

        time_grid = np.asarray(time_nuc, float)
        n_steps   = len(time_grid)

        if not (0 <= int(idx_shutdown) < n_steps):
            raise ValueError(
                f"idx_shutdown={idx_shutdown} out of bounds for n_steps={n_steps}"
            )
        if len(total_h) != n_steps:
            raise ValueError(
                f"Time grid mismatch: len(total_act)={len(total_h)} vs n_steps={n_steps} "
                f"for cell {cid}, mat {mat_id}."
            )

        total_h   = np.asarray(total_h, float)
        t_rel     = time_grid - time_grid[int(idx_shutdown)]
        mask      = t_rel > 0.0
        t_plot    = (t_rel[mask] / SECONDS_PER_YEAR).astype(float)
        h_plot    = total_h[mask]
        total_all[cid] = (t_plot, h_plot)

        if h_plot.size and np.nanmax(h_plot) > 0:
            max_heat = max(max_heat, float(np.nanmax(h_plot)))

        region = cell_id_to_name.get(cid, str(cid))
        if not _has_positive_finite(h_plot):
            print(f"[DecayHeat] Skipping plot for cell {cid} ({region}): zero decay heat.")
            continue

        topn_by_step, top_nucs_union = _select(heat_list, n_steps)
        top_nucs = sorted(top_nucs_union)

        nuc_series    = {nuc: np.full(n_steps, np.nan) for nuc in top_nucs}
        others_series = np.full(n_steps, np.nan)

        for istep in range(n_steps):
            d = heat_list[istep] or {}
            if not d:
                continue
            top_sum = 0.0
            for nuc, val in topn_by_step[istep]:
                nuc_series[nuc][istep] = float(val)
                top_sum += float(val)
            others_series[istep] = max(0.0, float(total_h[istep]) - top_sum)

        def _m(arr): return np.asarray(arr, float)[mask]

        others_plot = _m(others_series)
        hlists      = sorted((_sortable_half_life(n), n) for n in top_nucs)
        sorted_nucs = [n for _, n in hlists]
        colors      = generate_colors(len(top_nucs))

        fig, ax = plt.subplots(figsize=(11, 6))
        mc, lsc = _style_cycle()

        # y-axis limits
        others_pos = others_plot[others_plot > 0.0]
        total_pos  = h_plot[h_plot > 0.0]
        if others_pos.size > 0 and total_pos.size > 0:
            ax.set_ylim(
                10 ** math.floor(math.log10(np.min(others_pos))),
                10 ** math.ceil( math.log10(np.max(total_pos))),
            )

        for nuc in top_nucs:
            vp = _m(nuc_series[nuc])
            if not _has_positive_finite(vp):
                continue

            others_pos = others_plot[others_plot > 0.0]
            if others_pos.size > 0 and np.max(vp[vp > 0], initial=0.0) < 10 ** math.floor(math.log10(np.min(others_pos))):
                continue

            ax.loglog(t_plot, vp,
                      label=nuc + display_half_life(nuc),
                      color=colors[sorted_nucs.index(nuc) % len(colors)],
                      marker=next(mc), linestyle=next(lsc),
                      markersize=3, linewidth=1.3)

        if _has_positive_finite(others_plot):
            ax.loglog(t_plot, others_plot, label="Others",
                      linewidth=2, color="black", linestyle="--")
        if _has_positive_finite(h_plot):
            ax.loglog(t_plot, h_plot, color="black", linewidth=3,
                      label="Total", zorder=10)

        ax.set(xlabel="Time after irradiation [years]",
               ylabel=f"Decay heat [{decayheat_units}]",
               title=f"Decay heat — {region}")
        ax.set_xlim(1e-8, 2e3)
        _add_time_reference_lines(ax); _format_log_axes(ax)
        _finalize_legend(ax, ncol=2 if len(top_nucs) <= 24 else 3,
                         title="Components")
        fig.subplots_adjust(right=0.72)
        fig.savefig(out_dir / f"decayheat_{_safe_csv_label(region)}.png",
                    dpi=200, bbox_inches="tight")
        plt.close(fig)

    return total_all, float(max_heat), float(global_min_topn_edge) if np.isfinite(global_min_topn_edge) else 0.0


def plot_decayheat_all_cells(
    total_all: Dict[int, Tuple[np.ndarray, np.ndarray]],
    cell_id_to_name: Dict[int, str],
    max_heat: float,
    out_dir: Path,
    *,
    decayheat_units: str = "W/cm3",
    global_min_topn_edge: float = 0.0,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(11, 6))
    mc, lsc = _style_cycle()
    colors  = plt.cm.tab20(np.linspace(0, 1, max(1, min(len(total_all), 20))))

    for i, (cid, (t, y)) in enumerate(total_all.items()):
        yy = np.asarray(y, float)
        if not _has_positive_finite(yy):
            continue
        ax.loglog(t, yy, label=cell_id_to_name.get(cid, str(cid)),
                  color=colors[i % len(colors)],
                  marker=next(mc), linestyle=next(lsc),
                  markersize=3.5, linewidth=1.1)

    ax.set(xlabel="Time after irradiation [years]",
           ylabel=f"Decay heat [{decayheat_units}]",
           title="Total decay heat")
    ax.set_xlim(1e-8, 2e3)
    _add_time_reference_lines(ax); _format_log_axes(ax)
    ymin = _ymin_from_topn_edge(global_min_topn_edge)
    if max_heat > 0 and ymin:
        ax.set_ylim(ymin, 10 * max_heat)
    _finalize_legend(ax, title="Regions")
    fig.subplots_adjust(right=0.72)
    fig.savefig(out_dir / "decayheat_all_cells.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_decayheat_radial(
    total_all: Dict[int, Tuple[np.ndarray, np.ndarray]],
    cell_ids: Sequence[int],
    xcentroids: Sequence[float],
    out_dir: Path,
    *,
    idx_to_plot: Sequence[int] = (0, 4, 8, 12, 16, 20),
    decayheat_units: str = "W/cm3",
    global_min_topn_edge: float = 0.0,
    max_heat: float = 0.0,
) -> None:
    ref_cid = next((c for c in cell_ids if c in total_all), None)
    if ref_cid is None:
        return
    t_rel_y, _ = total_all[ref_cid]
    n_times     = len(t_rel_y)
    idx_to_plot = [int(i) for i in idx_to_plot if int(i) < n_times]
    if not idx_to_plot:
        return

    x      = np.asarray(xcentroids, float)
    radial = np.full((len(cell_ids), n_times), np.nan)
    for i, cid in enumerate(cell_ids):
        if cid in total_all:
            radial[i, :] = np.asarray(total_all[cid][1], float)

    out_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 6))
    mc, lsc = _style_cycle()
    colors  = plt.cm.tab20(np.linspace(0, 1, max(1, len(idx_to_plot))))

    for k, j in enumerate(idx_to_plot):
        y_safe = np.where(radial[:, j] > 0, radial[:, j], np.nan)
        ax.semilogy(x, y_safe, label=f"t = {t_rel_y[j]:.2e} y",
                    color=colors[k % len(colors)],
                    marker=next(mc), linestyle=next(lsc),
                    markersize=4, linewidth=1.5)

    ax.set(xlabel="Radial position [cm]",
           ylabel=f"Decay heat [{decayheat_units}]",
           title="Radial decay heat profiles")
    ax.grid(True, which="both", alpha=0.3)
    ymin = _ymin_from_topn_edge(global_min_topn_edge)
    if max_heat > 0 and ymin:
        ax.set_ylim(ymin, 10 * max_heat)
    _finalize_legend(ax, ncol=1 if len(idx_to_plot) <= 12 else 2,
                     title="Timesteps")
    fig.subplots_adjust(right=0.76)
    fig.savefig(out_dir / "radial_decayheat_profiles.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────────
# Chunk runner
# ──────────────────────────────────────────────────────────────────────────────
def run_chunk_postprocess(
    results: openmc.deplete.Results,
    chunks: Dict[str, Sequence[int]],
    *,
    xcentroids_by_chunk: Optional[Dict[str, Sequence[float]]] = None,
    layer_tags_by_chunk: Optional[Dict[str, Sequence[str]]]   = None,
    idx_shutdown: Optional[int] = None,
    activity_units: str = "Bq/kg",
    decayheat_units: str = "W/cm3",
    idx_to_plot: Sequence[int] = (0, 4, 8, 12, 16, 20),
) -> None:
    if xcentroids_by_chunk is None: xcentroids_by_chunk = {}
    if layer_tags_by_chunk  is None: layer_tags_by_chunk  = {}

    if idx_shutdown is None:
        sr    = results.get_source_rates()
        times = results.get_times()
        idx_shutdown = shutdown_index_from_source_rates(np.asarray(sr), len(times))

    mapping = load_depletion_mapping(results)

    for chunk_key, cell_ids in chunks.items():
        cell_ids = list(cell_ids)
        if not cell_ids:
            continue

        out_dir = RESULTS_OUT_DIR / chunk_key
        out_dir.mkdir(parents=True, exist_ok=True)

        cell_id_to_name = build_cell_id_to_name(
            cell_ids,
            layer_tags=layer_tags_by_chunk.get(chunk_key),
            chunk_key=chunk_key,
        )
        x = (xcentroids_by_chunk[chunk_key]
             if chunk_key in xcentroids_by_chunk and
             len(xcentroids_by_chunk[chunk_key]) == len(cell_ids)
             else np.arange(len(cell_ids), dtype=float))

        # Activity
        ta_all, act_min = plot_activity_nuclides_per_cell(
            results, mapping, cell_ids, cell_id_to_name, out_dir,
            idx_shutdown=idx_shutdown, activity_units=activity_units,
        )
        plot_activity_all_cells(ta_all, cell_id_to_name, out_dir,
                                activity_units=activity_units,
                                global_min_topn_edge=act_min)
        save_total_timeseries_csv(ta_all, cell_id_to_name,
                                  out_dir / "activity_all_cells.csv",
                                  "activity")
        plot_activity_radial(ta_all, cell_ids, x, out_dir,
                             idx_to_plot=idx_to_plot,
                             activity_units=activity_units,
                             global_min_topn_edge=act_min)

        # Decay heat 
        th_all, max_h, dh_min = plot_decayheat_nuclides_per_cell(
            results, mapping, cell_ids, cell_id_to_name, out_dir,
            idx_shutdown=idx_shutdown, decayheat_units=decayheat_units,
        )
        plot_decayheat_all_cells(th_all, cell_id_to_name, max_h, out_dir,
                                 decayheat_units=decayheat_units,
                                 global_min_topn_edge=dh_min)
        save_total_timeseries_csv(th_all, cell_id_to_name,
                                  out_dir / "decayheat_all_cells.csv",
                                  "decayheat")
        plot_decayheat_radial(th_all, cell_ids, x, out_dir,
                              idx_to_plot=idx_to_plot,
                              decayheat_units=decayheat_units,
                              global_min_topn_edge=dh_min, max_heat=max_h)

        print(f"[DONE] {chunk_key} → {out_dir}")


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # Chain file 
    openmc.config["chain_file"] = str(cfg.OPENMC_CHAIN_FILE)
    print(f"[chain] Using: {cfg.OPENMC_CHAIN_FILE}")

    # Verify required files 
    if not DEPLETION_RESULTS_H5.is_file():
        raise FileNotFoundError(f"Not found: {DEPLETION_RESULTS_H5}")
    if not CELL_MATERIAL_MAP_CSV.is_file():
        raise FileNotFoundError(f"Not found: {CELL_MATERIAL_MAP_CSV}")

    results = openmc.deplete.Results(str(DEPLETION_RESULTS_H5))

    # Chunk helpers from geometry (no model rebuild) 
    from geometry import build_breeder_chunks

    _chunk = build_breeder_chunks(
        cfg.INPUT_JSON, default_chunk_key=cfg.ALBEDO_CHUNK_KEY
    )
    ob_by_key           = _chunk["ob_by_key"]
    ib_by_key           = _chunk["ib_by_key"]
    radial_bins_for_key = _chunk["radial_bins_for_key"]

    # Resolve chunks from cfg.DEPLETION_CHUNKS
    CHUNKS: Dict[str, Sequence[int]] = {}
    for key, override_cells in cfg.DEPLETION_CHUNKS.items():
        if override_cells is not None:
            CHUNKS[key] = list(override_cells)
        elif key in ob_by_key:
            CHUNKS[key] = ob_by_key[key]
        elif key in ib_by_key:
            CHUNKS[key] = ib_by_key[key]
        else:
            print(f"[warn] DEPLETION_CHUNKS key {key!r} not found – skipping")

    # For slab: restrict to cells present in the depletion results
    if cfg.SIM_TYPE in ("slab", "fast_slab"):
        _present = set(
            int(k) for k in load_depletion_mapping(results).cell_to_mat.keys()
        )
        CHUNKS = {
            k: [c for c in v if c in _present]
            for k, v in CHUNKS.items()
        }

    # Layer tags and centroids 
    layer_tags_by_chunk:  Dict[str, Sequence[str]]   = {}
    xcentroids_by_chunk:  Dict[str, Sequence[float]] = {}

    for key, cell_ids in CHUNKS.items():
        if not cell_ids:
            continue
        n_breeder_layers = len(cell_ids) - 3
        side = "OB" if key.startswith("OB_") else "IB"
        layer_tags_by_chunk[key] = make_default_layer_tags(
            n_breeder_layers, breeder_prefix=side
        )
        try:
            xcentroids_by_chunk[key], _, _ = radial_bins_for_key(key)
        except Exception as e:
            print(f"[warn] centroids for {key}: {e}")

    # Shutdown index 
    sr           = results.get_source_rates()
    times        = results.get_times()
    idx_shutdown = shutdown_index_from_source_rates(np.asarray(sr), len(times))
    #print(f"[info] Shutdown index: {idx_shutdown} / {len(times)}")

    # Run 
    run_chunk_postprocess(
        results,
        CHUNKS,
        xcentroids_by_chunk = xcentroids_by_chunk,
        layer_tags_by_chunk = layer_tags_by_chunk,
        idx_shutdown        = idx_shutdown,
        activity_units      = "Bq/kg",
        decayheat_units     = "W/cm3",
        idx_to_plot         = cfg.DEPLETION_IDX_TO_PLOT,
    )
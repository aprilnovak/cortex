# ============================================================
# depletion_post.py
# Post-process OpenMC depletion results (activity + decay heat)
# for selected chunks (e.g., OB_1_b6) using the cell->material
# mapping CSV written by depletion_model.py.
# ============================================================
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple, Optional
from itertools import cycle

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import openmc
import openmc.deplete

# Some reference papers:
# https://doi.org/10.1016/j.fusengdes.2021.112428
# https://doi.org/10.1016/j.fusengdes.2021.112338
# https://iopscience.iop.org/article/10.1088/1741-4326/aca61f

# parameters
SECONDS_PER_YEAR = 365.25 * 24 * 3600.0

BASE_DIR = Path.cwd()
INPUT_JSON = (BASE_DIR / "EUDEMO_HCPB_inputs.json")

DEPLETION_RUN_DIR = (BASE_DIR / "depletion_run")
R2S_ACTIVATION_DIR = (DEPLETION_RUN_DIR / "r2s" / "activation")

DEPLETION_RESULTS_FILE = R2S_ACTIVATION_DIR / "depletion_results.h5"
CELL_MATERIAL_MAP_CSV = R2S_ACTIVATION_DIR / "cell_material_map.csv"

RESULTS_OUT_DIR = (BASE_DIR / "depletion_results")
RESULTS_OUT_DIR.mkdir(parents=True, exist_ok=True)

def generate_colors(n):
    """Generates a smooth rainbow gradient of n RGB colors."""
    cmap = plt.get_cmap('turbo') # plasma also looks nice
    color_range = cmap(np.linspace(1, 0, n))
    return color_range

def display_half_life(nuclide):
  half_life = openmc.data.half_life(nuclide)

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

# -----------------------------
# Layer tags
# -----------------------------
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
    if layer_tags is not None and len(layer_tags) == len(cell_ids):
        return {cid: str(layer_tags[i]) for i, cid in enumerate(cell_ids)}

    prefix = f"{chunk_key}_" if chunk_key else ""
    return {cid: f"{prefix}Layer_{i:02d}" for i, cid in enumerate(cell_ids)}

# -----------------------------
# Results I/O + mapping
# -----------------------------
@dataclass(frozen=True)
class DepletionMapping:
    cell_to_mat: Dict[int, str]            # cell_id -> mat_id (string)
    cells_by_mat: Dict[str, List[int]]     # mat_id -> [cell_id, ...]
    mat_id_to_name: Dict[str, str]         # mat_id -> name
    vol_by_mat: Dict[str, float]           # mat_id -> volume (from results[0].volume)

def load_depletion_mapping(
    results: openmc.deplete.Results,
    map_csv: str | Path = CELL_MATERIAL_MAP_CSV,
    ) -> DepletionMapping:
    """
    Loads depletion_run/r2s/activation/cell_material_map.csv and volumes from step0.
    """
    map_csv = Path(map_csv)
    if not map_csv.is_file():
        raise FileNotFoundError(f"Cell-material map CSV not found: {map_csv}")

    map_df = pd.read_csv(str(map_csv))
    map_df["mat_id"] = map_df["mat_id"].astype(str)

    cell_to_mat = dict(zip(map_df["cell_id"], map_df["mat_id"]))
    cells_by_mat = map_df.groupby("mat_id")["cell_id"].apply(list).to_dict()
    mat_id_to_name = map_df.groupby("mat_id")["material_name"].first().to_dict()

    step0 = results[0]
    vol_by_mat = step0.volume  # dict-like: mat_id(str) -> volume

    return DepletionMapping(
        cell_to_mat=cell_to_mat,
        cells_by_mat=cells_by_mat,
        mat_id_to_name=mat_id_to_name,
        vol_by_mat=vol_by_mat,
    )

def ensure_chunk_outdir(base_out_dir: str | Path, chunk_key: str) -> Path:
    out = Path(base_out_dir) / chunk_key
    out.mkdir(parents=True, exist_ok=True)
    return out

# -----------------------------
# CSV export helpers
# -----------------------------
def _safe_csv_label(text: str) -> str:
    return str(text).strip().replace(" ", "_").replace("/", "_")

def save_total_timeseries_csv(
    total_by_cell: Dict[int, Tuple[np.ndarray, np.ndarray]],
    cell_id_to_name: Dict[int, str],
    out_csv: str | Path,
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
        region = cell_id_to_name.get(cid, str(cid))
        col_name = f"{_safe_csv_label(region)}__cell_{cid}"

        df_i = pd.DataFrame({
            "Cooling_time_years": np.asarray(t_rel_y, dtype=float),
            col_name: np.asarray(values, dtype=float),
        })
        dfs.append(df_i)

    df_out = dfs[0]
    for df_i in dfs[1:]:
        df_out = df_out.merge(df_i, on="Cooling_time_years", how="outer")

    df_out = df_out.sort_values("Cooling_time_years").reset_index(drop=True)
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(out_csv, index=False)
    #print(f"[CSV] Wrote {value_label} totals to {out_csv}")

# -----------------------------
# Time utilities
# -----------------------------
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

# ============================================================
# Formatting helpers
# ============================================================
_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "<", ">", "h"]
#_LINESTYLES = ["-", "--", ":", "-."]
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

def _has_positive_finite(y: np.ndarray) -> bool:
    y = np.asarray(y, dtype=float)
    return bool(np.any(np.isfinite(y) & (y > 0.0)))

def _ymin_from_topn_edge(global_min_topn_edge: float) -> Optional[float]:
    try:
        v = float(global_min_topn_edge)
    except Exception:
        return None
    if not np.isfinite(v) or v <= 0.0:
        return None
    return 0.1 * v

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

    running_total = 0.0
    prev_running_total = 0.0
    index = 0
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

    return sorted_items[:i]

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

def _add_time_reference_lines(ax):
    """Useful reference lines (x-axis in years)."""

    SEC_PER_YEAR = 365.0 * 24.0 * 3600.0
    MIN_PER_YEAR = 365.0 * 24.0 * 60
    HOUR_PER_YEAR = 365.0 * 24.0
    WEEK_PER_YEAR = 52
    DAY_PER_YEAR = 365.0
    MONTH_PER_YEAR = 12.0

    refs = [
        (1.0 / SEC_PER_YEAR,   "1 s"),
        (1.0 / MIN_PER_YEAR,   "1 m"),
        (1.0 / HOUR_PER_YEAR,  "1 h"),
        (1.0 / WEEK_PER_YEAR,  "1 w"),
        (1.0 / DAY_PER_YEAR,   "1 d"),
        (1.0 / MONTH_PER_YEAR, "4 w"),
        (1.0,                 "1 y"),
        (10.0,                "10 y"),
        (100.0,               "100 y"),
        (1000.0,              "1000 y"),
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



# ============================================================
# Activity
# ============================================================
def plot_activity_nuclides_per_cell(
    results: openmc.deplete.Results,
    mapping,
    cell_ids: Sequence[int],
    cell_id_to_name: Dict[int, str],
    out_dir: Path,
    *,
    idx_shutdown: int,
    top_n: int = 5,
    activity_units: str = "Bq/kg",      # "Bq/kg" or "Bq"
    similarity_threshold: float = 0.10, # 10%
    print_similarity: bool = True,
    ) -> Tuple[Dict[int, Tuple[np.ndarray, np.ndarray]], float]:
    """
    Returns:
      total_activity_all[cid] = (t_rel_years, total_activity_in_units)
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    total_activity_all: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    global_min_topn_edge = np.inf

    for cid in cell_ids:
        if cid not in mapping.cell_to_mat:
            print(f"[Activity] WARNING: cell {cid} not found in map; skipping.")
            continue

        mat_id = str(mapping.cell_to_mat[cid])
        if mat_id not in mapping.vol_by_mat:
            print(f"[Activity] WARNING: mat_id {mat_id} missing volume; skipping cell {cid}.")
            continue

        vol = mapping.vol_by_mat[mat_id]

        # Total activity
        _, total_act = results.get_activity(mat=mat_id, by_nuclide=False, units=activity_units, volume=vol)

        # Per-nuclide activity
        time_nuc, act_dict_list = results.get_activity(mat=mat_id, by_nuclide=True, units=activity_units, volume=vol)

        time_grid = np.asarray(time_nuc, dtype=float)
        n_steps = len(time_grid)

        if not (0 <= int(idx_shutdown) < n_steps):
            raise ValueError(f"[Activity] idx_shutdown={idx_shutdown} out of bounds for n_steps={n_steps}")

        total_act = np.asarray(total_act, dtype=float)
        if len(total_act) != n_steps:
            raise ValueError(
                f"[Activity] Time grid mismatch total vs nuclide activity "
                f"(len(total_act)={len(total_act)} vs n_steps={n_steps}) for cell {cid}, mat {mat_id}."
            )

        t_rel = time_grid - time_grid[int(idx_shutdown)]
        mask = t_rel > 0.0

        t_rel_plot = (t_rel[mask] / SECONDS_PER_YEAR).astype(float)
        total_act_plot = total_act[mask]
        total_activity_all[cid] = (t_rel_plot, total_act_plot)

        topn_list_by_step, top_nucs_union = _select(act_dict_list, n_steps)
        top_nucs = sorted(top_nucs_union)

        nuc_series: Dict[str, np.ndarray] = {nuc: np.full(n_steps, np.nan, dtype=float) for nuc in top_nucs}
        others_series = np.full(n_steps, np.nan, dtype=float)

        for istep in range(n_steps):
            d = act_dict_list[istep] or {}
            if not d:
                continue

            sel = topn_list_by_step[istep]
            top_sum = 0.0
            for nuc, val in sel:
                nuc_series[nuc][istep] = float(val)
                top_sum += float(val)

            tot_val = float(total_act[istep])
            others_series[istep] = max(0.0, tot_val - top_sum)

        def _mask(arr: np.ndarray) -> np.ndarray:
            return np.asarray(arr, dtype=float)[mask]

        others_plot = _mask(others_series)

        # Plot per-cell breakdown
        fig, ax = plt.subplots(figsize=(11, 6))
        mcycle, lscycle = _style_cycle()

        # Nuclides (markers + varied linestyles)
        halflife_list = [openmc.data.half_life(i) for i in top_nucs]
        nuclide_list = top_nucs
        sorted_halflife_list, sorted_nuclide_list = zip(*sorted(zip(halflife_list, nuclide_list)))
        colors = generate_colors(len(top_nucs))

        # y-min cutoff; cut off at 1 order of magnitude below the other activity and 1 order of magnitude
        # above the total activity (rounded to powers of 10)
        min_act = 10 ** math.floor(math.log10(np.min(others_plot)))
        max_act = 10 ** math.ceil(math.log10(np.max(total_act_plot)))
        ax.set_ylim([min_act, max_act])

        for i, nuc in enumerate(top_nucs):
            vals_plot = _mask(nuc_series[nuc])
            if not _has_positive_finite(vals_plot):
                continue

            if np.max(vals_plot) < min_act:
                continue

            ax.loglog(
                t_rel_plot,
                vals_plot,
                label=nuc + display_half_life(nuc),
                color=colors[sorted_nuclide_list.index(nuc) % len(colors)],
                marker=next(mcycle),
                linestyle=next(lscycle),
                markersize=3,
                linewidth=1.3,
            )

        # Others
        if _has_positive_finite(others_plot):
            ax.loglog(t_rel_plot, others_plot, label="Others", linewidth=2.0, color='black', linestyle='--')

        # Total
        if _has_positive_finite(total_act_plot):
            ax.loglog(t_rel_plot, total_act_plot, color="black", linewidth=2.0, label="Total", zorder=10)

        region = cell_id_to_name.get(cid, str(cid))
        ax.set_xlabel("Time after irradiation [years]")
        ax.set_ylabel(f"Activity [{activity_units}]")
        ax.set_title(f"Activity — {region} (mat {mat_id})")

        _add_time_reference_lines(ax)
        _format_log_axes(ax)

        # Legend
        n_entries = len(ax.get_legend_handles_labels()[1])
        ncol = 2 if n_entries <= 24 else 3
        _finalize_legend(ax, ncol=ncol, title="Components", fontsize=8)

        fig.subplots_adjust(right=0.72)
        safe_region = str(region).replace(" ", "_").replace("/", "_")
        fig.savefig(out_dir / f"activity_cell_{safe_region}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

    if not np.isfinite(global_min_topn_edge):
        global_min_topn_edge = 0.0

    return total_activity_all, float(global_min_topn_edge)

def plot_activity_all_cells(
    total_activity_all: Dict[int, Tuple[np.ndarray, np.ndarray]],
    cell_id_to_name: Dict[int, str],
    out_dir: Path,
    *,
    activity_units: str = "Bq/kg",
    global_min_topn_edge: float = 0.0,
    ):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(11, 6))
    mcycle, lscycle = _style_cycle()

    num_cells = max(1, len(total_activity_all))
    colors = plt.cm.tab20(np.linspace(0, 1, max(1, min(num_cells, 20))))

    ymax = 0.0
    for i, (cid, (t_rel_y, y)) in enumerate(total_activity_all.items()):
        region = cell_id_to_name.get(cid, str(cid))
        yy = np.asarray(y, float)
        if yy.size:
            ymax = max(ymax, float(np.nanmax(yy)))

        ax.loglog(
            t_rel_y,
            yy,
            label=region,
            color=colors[i % len(colors)],
            marker=next(mcycle),
            linestyle=next(lscycle),
            markersize=3.5,
            linewidth=1.1,
        )

    ax.set_xlabel("Time after irradiation [years]")
    ax.set_ylabel(f"Activity [{activity_units}]")
    ax.set_title("Total Activity over cooling time")

    _add_time_reference_lines(ax)
    _format_log_axes(ax)

    ymin = _ymin_from_topn_edge(global_min_topn_edge)
    if ymin is not None and ymax > 0.0:
        ax.set_ylim(ymin, 10.0 * ymax)

    n_entries = len(ax.get_legend_handles_labels()[1])
    ncol = 2 if n_entries <= 24 else 3
    _finalize_legend(ax, ncol=ncol, title="Regions", fontsize=8)

    fig.subplots_adjust(right=0.72)
    fig.savefig(out_dir / "activity_all_cells.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

def plot_activity_radial_profiles(
    total_activity_all: Dict[int, Tuple[np.ndarray, np.ndarray]],
    cell_ids: Sequence[int],
    xcentroids: Sequence[float],
    out_dir: Path,
    *,
    idx_to_plot: Sequence[int] = (0, 3, 6, 9, 12, 15),
    activity_units: str = "Bq/kg",
    global_min_topn_edge: float = 0.0,
    ):
    if not cell_ids:
        return

    ref_cid = next((cid for cid in cell_ids if cid in total_activity_all), None)
    if ref_cid is None:
        return

    t_rel_y, _ = total_activity_all[ref_cid]
    n_times = len(t_rel_y)

    idx_to_plot = [int(i) for i in idx_to_plot if int(i) < n_times]
    if not idx_to_plot:
        return

    x = np.asarray(xcentroids, dtype=float)
    n_cells = len(cell_ids)
    radial = np.full((n_cells, n_times), np.nan, dtype=float)

    ymax = 0.0
    for i, cid in enumerate(cell_ids):
        if cid not in total_activity_all:
            continue
        _, act = total_activity_all[cid]
        radial[i, :] = np.asarray(act, dtype=float)
        if np.asarray(act).size:
            ymax = max(ymax, float(np.nanmax(act)))

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    mcycle, lscycle = _style_cycle()

    colors = plt.cm.tab20(np.linspace(0, 1, max(1, len(idx_to_plot))))

    for k, j in enumerate(idx_to_plot):
        y = radial[:, j]
        ax.semilogy(
            x,
            y,
            label=f"t = {t_rel_y[j]:.2e} y",
            color=colors[k % len(colors)],
            marker=next(mcycle),
            linestyle=next(lscycle),
            markersize=4,
            linewidth=1.5,
        )

    ax.set_xlabel("Radial position [cm]")
    ax.set_ylabel(f"Total activity [{activity_units}]")
    ax.set_title("Radial profile of Total activity per cooling time")

    ax.grid(True, which="major", linewidth=0.8, alpha=0.35)
    ax.grid(True, which="minor", linewidth=0.5, alpha=0.15)

    ymin = _ymin_from_topn_edge(global_min_topn_edge)
    if ymin is not None and ymax > 0.0:
        ax.set_ylim(ymin, 10.0 * ymax)

    n_entries = len(ax.get_legend_handles_labels()[1])
    ncol = 1 if n_entries <= 12 else 2
    _finalize_legend(ax, ncol=ncol, title="Timesteps", fontsize=8)

    fig.subplots_adjust(right=0.76)
    fig.savefig(out_dir / "radial_activity_profiles.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

# ============================================================
# Decay heat
# ============================================================
def plot_decayheat_nuclides_per_cell(
    results: openmc.deplete.Results,
    mapping,
    cell_ids: Sequence[int],
    cell_id_to_name: Dict[int, str],
    out_dir: Path,
    *,
    idx_shutdown: int,
    top_n: int = 10,
    decayheat_units: str = "W/cm3",
    similarity_threshold: float = 0.10,
    print_similarity: bool = True,
    ) -> Tuple[Dict[int, Tuple[np.ndarray, np.ndarray]], float, float]:
    """
    Returns:
      total_heat_all[cid] = (t_rel_years, total_heat_in_units)
      max_decay_comb = global max for y-max scaling
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    total_heat_all: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    max_decay_comb = 0.0
    global_min_topn_edge = np.inf

    for cid in cell_ids:
        if cid not in mapping.cell_to_mat:
            print(f"[DecayHeat] WARNING: cell {cid} not found in map; skipping.")
            continue

        mat_id = str(mapping.cell_to_mat[cid])
        if mat_id not in mapping.vol_by_mat:
            print(f"[DecayHeat] WARNING: mat_id {mat_id} missing volume; skipping cell {cid}.")
            continue

        vol = mapping.vol_by_mat[mat_id]

        _, total_h = results.get_decay_heat(mat=mat_id, by_nuclide=False, units=decayheat_units, volume=vol)
        time_nuc_h, heat_dict_list = results.get_decay_heat(mat=mat_id, by_nuclide=True, units=decayheat_units, volume=vol)

        time_grid = np.asarray(time_nuc_h, dtype=float)
        n_steps = len(time_grid)

        if not (0 <= int(idx_shutdown) < n_steps):
            raise ValueError(f"[DecayHeat] idx_shutdown={idx_shutdown} out of bounds for n_steps={n_steps}")

        total_h = np.asarray(total_h, dtype=float)
        if len(total_h) != n_steps:
            raise ValueError(
                f"[DecayHeat] Time grid mismatch total vs nuclide heat "
                f"(len(total_h)={len(total_h)} vs n_steps={n_steps}) for cell {cid}, mat {mat_id}."
            )

        t_rel = time_grid - time_grid[int(idx_shutdown)]
        mask = t_rel > 0.0

        t_rel_plot = (t_rel[mask] / SECONDS_PER_YEAR).astype(float)
        total_h_plot = total_h[mask]
        total_heat_all[cid] = (t_rel_plot, total_h_plot)

        if total_h_plot.size and np.nanmax(total_h_plot) > 0.0:
            max_decay_comb = max(max_decay_comb, float(np.nanmax(total_h_plot)))

        topn_list_by_step, top_nucs_union = _select(heat_dict_list, n_steps)
        top_nucs = sorted(top_nucs_union)

        nuc_series: Dict[str, np.ndarray] = {nuc: np.full(n_steps, np.nan, dtype=float) for nuc in top_nucs}
        others_series = np.full(n_steps, np.nan, dtype=float)

        for istep in range(n_steps):
            d = heat_dict_list[istep] or {}
            if not d:
                continue

            sel = topn_list_by_step[istep]
            top_sum = 0.0
            for nuc, val in sel:
                nuc_series[nuc][istep] = float(val)
                top_sum += float(val)

            tot_val = float(total_h[istep])
            others_series[istep] = max(0.0, tot_val - top_sum)

        def _mask(arr: np.ndarray) -> np.ndarray:
            return np.asarray(arr, dtype=float)[mask]

        others_plot = _mask(others_series)

        # Plot per-cell breakdown (formatting 1,2,3)
        fig, ax = plt.subplots(figsize=(11, 6))
        mcycle, lscycle = _style_cycle()

        # form a color scale based on the half lives
        halflife_list = [openmc.data.half_life(i) for i in top_nucs]
        nuclide_list = top_nucs
        sorted_halflife_list, sorted_nuclide_list = zip(*sorted(zip(halflife_list, nuclide_list)))
        colors = generate_colors(len(top_nucs))

        for i, nuc in enumerate(top_nucs):
            vals_plot = _mask(nuc_series[nuc])
            if not _has_positive_finite(vals_plot):
                continue

            ax.loglog(
                t_rel_plot,
                vals_plot,
                label=nuc + display_half_life(nuc),
                color=colors[sorted_nuclide_list.index(nuc) % len(colors)],
                marker=next(mcycle),
                linestyle=next(lscycle),
                markersize=3,
                linewidth=1.3,
            )

        if _has_positive_finite(others_plot):
            ax.loglog(t_rel_plot, others_plot, label="Others", linewidth=2.0, color='black', linestyle='--')

        if _has_positive_finite(total_h_plot):
            ax.loglog(t_rel_plot, total_h_plot, color="black", linewidth=3.0, label="Total", zorder=10)

        region = cell_id_to_name.get(cid, str(cid))
        ax.set_xlabel("Time after irradiation [years]")
        ax.set_ylabel(f"Decay heat [{decayheat_units}]")
        ax.set_title(f"Decay heat — {region} (mat {mat_id})")

        # y-min cutoff; cut off at 1 order of magnitude below the other decay heat and 1 order of magnitude
        # above the total decay heat (rounded to powers of 10)
        ax.set_ylim([10 ** math.floor(math.log10(np.min(others_plot))), 10 ** math.ceil(math.log10(np.max(total_h_plot)))])

        _add_time_reference_lines(ax)
        _format_log_axes(ax)

        if total_h_plot.size and np.nanmax(total_h_plot) > 0.0:
            ymax = float(np.nanmax(total_h_plot))
            ymin = _ymin_from_topn_edge(global_min_topn_edge)
            ax.set_xlim(1e-8, 2e3)
            if ymin is not None:
                ax.set_ylim(ymin, 10.0 * ymax)

        n_entries = len(ax.get_legend_handles_labels()[1])
        ncol = 2 if n_entries <= 24 else 3
        _finalize_legend(ax, ncol=ncol, title="Components", fontsize=8)

        fig.subplots_adjust(right=0.72)
        safe_region = str(region).replace(" ", "_").replace("/", "_")
        fig.savefig(out_dir / f"decayheat_{safe_region}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

    if not np.isfinite(global_min_topn_edge):
        global_min_topn_edge = 0.0

    return total_heat_all, float(max_decay_comb), float(global_min_topn_edge)

def plot_decayheat_all_cells(
    total_heat_all: Dict[int, Tuple[np.ndarray, np.ndarray]],
    cell_id_to_name: Dict[int, str],
    max_decay_comb: float,
    out_dir: Path,
    *,
    decayheat_units: str = "W/cm3",
    global_min_topn_edge: float = 0.0,
    ):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(11, 6))
    mcycle, lscycle = _style_cycle()

    num_cells = max(1, len(total_heat_all))
    colors = plt.cm.tab20(np.linspace(0, 1, max(1, min(num_cells, 20))))

    for i, (cid, (t_rel_y, y)) in enumerate(total_heat_all.items()):
        region = cell_id_to_name.get(cid, str(cid))
        yy = np.asarray(y, float)
        ax.loglog(
            t_rel_y,
            yy,
            label=region,
            color=colors[i % len(colors)],
            marker=next(mcycle),
            linestyle=next(lscycle),
            markersize=3.5,
            linewidth=1.1,
        )

    ax.set_xlabel("Time after irradiation [years]")
    ax.set_ylabel(f"Decay heat [{decayheat_units}]")
    ax.set_title("Total decay heat")

    _add_time_reference_lines(ax)
    _format_log_axes(ax)

    ax.set_xlim(1e-8, 2e3)
    ymin = _ymin_from_topn_edge(global_min_topn_edge)
    if max_decay_comb > 0 and ymin is not None:
        ax.set_ylim(ymin, 10.0 * max_decay_comb)

    n_entries = len(ax.get_legend_handles_labels()[1])
    ncol = 2 if n_entries <= 24 else 3
    _finalize_legend(ax, ncol=ncol, title="Regions", fontsize=8)

    fig.subplots_adjust(right=0.72)
    fig.savefig(out_dir / "decayheat_all_cells.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

def plot_decayheat_radial_profiles(
    total_heat_all: Dict[int, Tuple[np.ndarray, np.ndarray]],
    cell_ids: Sequence[int],
    xcentroids: Sequence[float],
    out_dir: Path,
    *,
    idx_to_plot: Sequence[int] = (0, 3, 6, 9, 12, 15),
    decayheat_units: str = "W/cm3",
    global_min_topn_edge: float = 0.0,
    max_decay_comb: float = 0.0,
    ):
    if not cell_ids:
        return

    ref_cid = next((cid for cid in cell_ids if cid in total_heat_all), None)
    if ref_cid is None:
        return

    t_rel_y, _ = total_heat_all[ref_cid]
    n_times = len(t_rel_y)

    idx_to_plot = [int(i) for i in idx_to_plot if int(i) < n_times]
    if not idx_to_plot:
        return

    x = np.asarray(xcentroids, dtype=float)
    n_cells = len(cell_ids)
    radial = np.full((n_cells, n_times), np.nan, dtype=float)

    for i, cid in enumerate(cell_ids):
        if cid not in total_heat_all:
            continue
        _, h = total_heat_all[cid]
        radial[i, :] = np.asarray(h, dtype=float)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    mcycle, lscycle = _style_cycle()
    colors = plt.cm.tab20(np.linspace(0, 1, max(1, len(idx_to_plot))))

    for k, j in enumerate(idx_to_plot):
        y = radial[:, j]
        ax.semilogy(
            x,
            y,
            label=f"t = {t_rel_y[j]:.2e} y",
            color=colors[k % len(colors)],
            marker=next(mcycle),
            linestyle=next(lscycle),
            markersize=4,
            linewidth=1.5,
        )

    ax.set_xlabel("Radial position [cm]")
    ax.set_ylabel(f"Decay heat [{decayheat_units}]")
    ax.set_title("Radial profile of decay heat per cooling time")
    ax.grid(True, which="major", linewidth=0.8, alpha=0.35)
    ax.grid(True, which="minor", linewidth=0.5, alpha=0.15)

    ymin = _ymin_from_topn_edge(global_min_topn_edge)
    if max_decay_comb > 0 and ymin is not None:
        ax.set_ylim(ymin, 10.0 * max_decay_comb)

    n_entries = len(ax.get_legend_handles_labels()[1])
    ncol = 1 if n_entries <= 12 else 2
    _finalize_legend(ax, ncol=ncol, title="Timesteps", fontsize=8)

    fig.subplots_adjust(right=0.76)
    fig.savefig(out_dir / "radial_decayheat_profiles.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

# ============================================================
# Chunk postprocessing
# ============================================================
def run_chunk_postprocess(
    results: openmc.deplete.Results,
    chunks: Dict[str, Sequence[int]],
    *,
    base_out_dir: str | Path = RESULTS_OUT_DIR,
    xcentroids_by_chunk: Optional[Dict[str, Sequence[float]]] = None,
    layer_tags_by_chunk: Optional[Dict[str, Sequence[str]]] = None,
    source_rates: Optional[np.ndarray] = None,
    idx_shutdown: Optional[int] = None,
    activity_units: str = "Bq/kg",
    decayheat_units: str = "W/cm3",
    activity_top_n: int = 5,
    decayheat_top_n: int = 5,
    similarity_threshold: float = 0.10,
    idx_to_plot: Sequence[int] = (0, 4, 8, 12, 16, 20),
    ):
    if xcentroids_by_chunk is None:
        xcentroids_by_chunk = {}
    if layer_tags_by_chunk is None:
        layer_tags_by_chunk = {}

    # Shutdown index once
    if idx_shutdown is None:
        if source_rates is None:
            source_rates = results.get_source_rates()
        times = results.get_times()
        idx_shutdown = shutdown_index_from_source_rates(source_rates, n_steps=len(times))

    mapping = load_depletion_mapping(results)

    for chunk_key, cell_ids in chunks.items():
        cell_ids = list(cell_ids)
        if not cell_ids:
            continue

        out_dir = ensure_chunk_outdir(base_out_dir, chunk_key)

        cell_id_to_name = build_cell_id_to_name(
            cell_ids,
            layer_tags=layer_tags_by_chunk.get(chunk_key, None),
            chunk_key=chunk_key,
        )

        if chunk_key in xcentroids_by_chunk and len(xcentroids_by_chunk[chunk_key]) == len(cell_ids):
            x = xcentroids_by_chunk[chunk_key]
        else:
            x = np.arange(len(cell_ids), dtype=float)

        # Activity
        total_activity_all, act_min_topn_edge = plot_activity_nuclides_per_cell(
            results, mapping, cell_ids, cell_id_to_name,
            out_dir=out_dir,
            idx_shutdown=int(idx_shutdown),
            top_n=activity_top_n,
            activity_units=activity_units,
            similarity_threshold=similarity_threshold,
        )

        plot_activity_all_cells(
            total_activity_all, cell_id_to_name,
            out_dir=out_dir,
            activity_units=activity_units,
            global_min_topn_edge=act_min_topn_edge,
        )

        save_total_timeseries_csv(
            total_activity_all,
            cell_id_to_name,
            out_dir / "activity_all_cells.csv",
            value_label="activity",
        )

        plot_activity_radial_profiles(
            total_activity_all, cell_ids, x,
            out_dir=out_dir,
            idx_to_plot=idx_to_plot,
            activity_units=activity_units,
            global_min_topn_edge=act_min_topn_edge,
        )

        # Decay heat
        total_heat_all, max_decay_comb, dh_min_topn_edge = plot_decayheat_nuclides_per_cell(
            results, mapping, cell_ids, cell_id_to_name,
            out_dir=out_dir,
            idx_shutdown=int(idx_shutdown),
            top_n=decayheat_top_n,
            decayheat_units=decayheat_units,
            similarity_threshold=similarity_threshold,
        )

        plot_decayheat_all_cells(
            total_heat_all, cell_id_to_name,
            max_decay_comb=max_decay_comb,
            out_dir=out_dir,
            decayheat_units=decayheat_units,
            global_min_topn_edge=dh_min_topn_edge,
        )

        save_total_timeseries_csv(
            total_heat_all,
            cell_id_to_name,
            out_dir / "decayheat_all_cells.csv",
            value_label="decayheat",
        )

        plot_decayheat_radial_profiles(
            total_heat_all, cell_ids, x,
            out_dir=out_dir,
            idx_to_plot=idx_to_plot,
            decayheat_units=decayheat_units,
            global_min_topn_edge=dh_min_topn_edge,
            max_decay_comb=max_decay_comb,
        )

        print(f"[DONE] {chunk_key} -> {out_dir}")


# ============================================================
# Run
# ============================================================
if __name__ == "__main__":
    print(f"[chain] openmc.config['chain_file'] = {openmc.config.get('chain_file', 'NOT SET')}")

    if not DEPLETION_RESULTS_FILE.is_file():
        raise FileNotFoundError(f"Depletion results file not found: {DEPLETION_RESULTS_FILE}")

    if not CELL_MATERIAL_MAP_CSV.is_file():
        raise FileNotFoundError(f"Cell-material map CSV not found: {CELL_MATERIAL_MAP_CSV}")

    results = openmc.deplete.Results(str(DEPLETION_RESULTS_FILE))


    from neutronics_model import build_breeder_chunks

    _cells_chunk = build_breeder_chunks(
        INPUT_JSON,
        default_chunk_key="OB_1_b6",
        gap_cm=2.0,
        start_cm=0.0,
    )

    data      = _cells_chunk["data"]
    geom      = _cells_chunk["geom"]
    cell_ids  = _cells_chunk["cell_ids_all"]
    ob_by_key = _cells_chunk["ob_by_key"]
    ib_by_key = _cells_chunk["ib_by_key"]
    OB_CHUNK_SIZE = _cells_chunk["OB_CHUNK_SIZE"]
    IB_CHUNK_SIZE = _cells_chunk["IB_CHUNK_SIZE"]
    n_breeder  = int(geom["n_breeder"])

    cell_ids_for_key    = _cells_chunk["cell_ids_for_key"]
    radial_bins_for_key = _cells_chunk["radial_bins_for_key"]

    OB_KEY = "OB_1_b6"
    IB_KEY = "IB_1_b4"

    xcentroids_ob, _, _ = radial_bins_for_key(OB_KEY)
    xcentroids_ib, _, _ = radial_bins_for_key(IB_KEY)

    CHUNKS = {
        OB_KEY: ob_by_key.get(OB_KEY, []),
        IB_KEY: ib_by_key.get(IB_KEY, []),
    }

    layer_tags_by_chunk: Dict[str, Sequence[str]] = {}
    xcentroids_by_chunk: Dict[str, Sequence[float]] = {}

    if CHUNKS[OB_KEY]:
        n_ob_layers = len(CHUNKS[OB_KEY]) - 3
        layer_tags_by_chunk[OB_KEY] = make_default_layer_tags(n_ob_layers, breeder_prefix="OB")
        xcentroids_by_chunk[OB_KEY], _, _ = radial_bins_for_key(OB_KEY)

    if CHUNKS[IB_KEY]:
        n_ib_layers = len(CHUNKS[IB_KEY]) - 3
        layer_tags_by_chunk[IB_KEY] = make_default_layer_tags(n_ib_layers, breeder_prefix="IB")
        xcentroids_by_chunk[IB_KEY], _, _ = radial_bins_for_key(IB_KEY)

    source_rates = results.get_source_rates()
    times = results.get_times()
    idx_shutdown = shutdown_index_from_source_rates(source_rates, n_steps=len(times))

    run_chunk_postprocess(
        results,
        CHUNKS,
        base_out_dir=RESULTS_OUT_DIR,
        xcentroids_by_chunk=xcentroids_by_chunk,
        layer_tags_by_chunk=layer_tags_by_chunk,
        idx_shutdown=int(idx_shutdown),
        activity_units="Bq/kg",       # "Bq" or "Bq/kg"
        decayheat_units="W/cm3",
        activity_top_n=10,
        decayheat_top_n=5,
        similarity_threshold=0.10,
        idx_to_plot=(0, 4, 8, 12, 16, 20),
    )

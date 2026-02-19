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

import time
from collections import OrderedDict

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
INPUT_JSON = Path("Tokamak_inputs.json")

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
        tags.append("First Wall")

    tags.extend([f"{breeder_prefix}_{i}" for i in range(1, n_breeder_layers + 1)])

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
    map_csv: str | Path = "r2s/activation/cell_material_map.csv",
    ) -> DepletionMapping:
    """
    Loads r2s/activation/cell_material_map.csv and volumes from step0.
    """
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


class Timer:
    def __init__(self):
        self._start = {}
        self.elapsed = OrderedDict()

    def start(self, name: str):
        self._start[name] = time.perf_counter()

    def stop(self, name: str):
        if name not in self._start:
            raise RuntimeError(f"Timer '{name}' was not started")
        dt = time.perf_counter() - self._start.pop(name)
        self.elapsed[name] = self.elapsed.get(name, 0.0) + dt

    def summary(self):
        total = sum(self.elapsed.values())
        print("\n=== Timing summary ===")
        for k, v in self.elapsed.items():
            print(f"{k:30s}: {v:8.3f} s ({100*v/total:5.1f}%)")
        print(f"{'TOTAL':30s}: {total:8.3f} s")


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

def _select_topn_with_tie(
    sorted_items: list[tuple[str, float]],
    *,
    top_n: int,
    similarity_threshold: float,
    ) -> list[tuple[str, float]]:
    """
    Always include top_n.
    If (n+1) within similarity_threshold of nth (relative to nth), include (n+1) too.
    """
    if top_n <= 0 or not sorted_items:
        return []

    base = sorted_items[:top_n]
    if len(sorted_items) >= top_n + 1:
        nth_val = float(sorted_items[top_n - 1][1])
        n1_val = float(sorted_items[top_n][1])
        if nth_val > 0.0:
            rel_diff = abs(nth_val - n1_val) / nth_val
            if rel_diff < similarity_threshold:
                base.append(sorted_items[top_n])
    return base

def _add_time_reference_lines(ax):
    """Useful reference lines (x-axis in years)."""

    SEC_PER_YEAR = 365.0 * 24.0 * 3600.0
    HOUR_PER_YEAR = 365.0 * 24.0
    DAY_PER_YEAR = 365.0
    MONTH_PER_YEAR = 12.0  

    refs = [
        (1.0 / SEC_PER_YEAR,   "1 s"),
        (1.0 / HOUR_PER_YEAR,  "1 h"),
        (1.0 / DAY_PER_YEAR,   "1 d"),
        (1.0 / MONTH_PER_YEAR, "1 m"),
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

        # selections per step
        topn_list_by_step: list[list[tuple[str, float]]] = [[] for _ in range(n_steps)]
        top_nucs_union: set[str] = set()

        for istep in range(n_steps):
            d = act_dict_list[istep] or {}
            if not d:
                continue

            step_sorted = sorted(d.items(), key=lambda x: x[1], reverse=True)

            # global min of true Nth in cooling portion
            if mask[istep] and top_n > 0 and len(step_sorted) >= top_n:
                nth_val = float(step_sorted[top_n - 1][1])
                if np.isfinite(nth_val) and nth_val > 0.0:
                    global_min_topn_edge = min(global_min_topn_edge, nth_val)

            selected = _select_topn_with_tie(
                step_sorted,
                top_n=int(top_n),
                similarity_threshold=float(similarity_threshold),
            )
            topn_list_by_step[istep] = selected

            if (
                print_similarity
                and top_n > 0
                and len(step_sorted) >= top_n + 1
                and len(selected) == top_n + 1
            ):
                nth_nuc, nth_val = step_sorted[top_n - 1]
                n1_nuc, n1_val = step_sorted[top_n]
                if float(nth_val) > 0.0:
                    rel_diff = abs(float(nth_val) - float(n1_val)) / float(nth_val)
                    region = cell_id_to_name.get(cid, str(cid))
                    t_years = float((time_grid[istep] - time_grid[int(idx_shutdown)]) / SECONDS_PER_YEAR)
                    #print(
                    #    f"[Activity] Included (n+1) due to tie in {region} (cell {cid}, mat {mat_id}) "
                    #    f"at step {istep} (t_rel={t_years:.3e} y): "
                    #    f"nth {nth_nuc}={float(nth_val):.3e} vs (n+1) {n1_nuc}={float(n1_val):.3e} "
                    #    f"(rel diff={rel_diff:.3%} < {similarity_threshold:.1%})"
                    #)

            for nuc, _ in selected:
                top_nucs_union.add(nuc)

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
        num_nucs = max(1, len(top_nucs))
        colors = plt.cm.tab20(np.linspace(0, 1, num_nucs))

        for i, nuc in enumerate(top_nucs):
            vals_plot = _mask(nuc_series[nuc])
            if not _has_positive_finite(vals_plot):
                continue

            ax.loglog(
                t_rel_plot,
                vals_plot,
                label=nuc,
                color=colors[i % len(colors)],
                marker=next(mcycle),
                linestyle=next(lscycle),
                markersize=4,
                linewidth=1.3,
            )

        # Others
        if _has_positive_finite(others_plot):
            ax.loglog(t_rel_plot, others_plot, label="Others", linewidth=2.0)

        # Total
        if _has_positive_finite(total_act_plot):
            ax.loglog(t_rel_plot, total_act_plot, color="black", linewidth=3.0, label="Total", zorder=10)

        region = cell_id_to_name.get(cid, str(cid))
        ax.set_xlabel("Time after irradiation [years]")
        ax.set_ylabel(f"Activity [{activity_units}]")
        ax.set_title(f"Activity (Top-{top_n} per timestep) — {region} (mat {mat_id})")

        _add_time_reference_lines(ax)
        _format_log_axes(ax)

        # y-min cutoff
        if total_act_plot.size and np.nanmax(total_act_plot) > 0.0:
            ymax = float(np.nanmax(total_act_plot))
            ymin = _ymin_from_topn_edge(global_min_topn_edge)
            if ymin is not None:
                ax.set_ylim(ymin, 10.0 * ymax)

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
    top_n: int = 5,
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

        topn_list_by_step: list[list[tuple[str, float]]] = [[] for _ in range(n_steps)]
        top_nucs_union: set[str] = set()

        for istep in range(n_steps):
            d = heat_dict_list[istep] or {}
            if not d:
                continue

            step_sorted = sorted(d.items(), key=lambda x: x[1], reverse=True)

            if mask[istep] and top_n > 0 and len(step_sorted) >= top_n:
                nth_val = float(step_sorted[top_n - 1][1])
                if np.isfinite(nth_val) and nth_val > 0.0:
                    global_min_topn_edge = min(global_min_topn_edge, nth_val)

            selected = _select_topn_with_tie(
                step_sorted,
                top_n=int(top_n),
                similarity_threshold=float(similarity_threshold),
            )
            topn_list_by_step[istep] = selected

            if (
                print_similarity
                and top_n > 0
                and len(step_sorted) >= top_n + 1
                and len(selected) == top_n + 1
            ):
                nth_nuc, nth_val = step_sorted[top_n - 1]
                n1_nuc, n1_val = step_sorted[top_n]
                if float(nth_val) > 0.0:
                    rel_diff = abs(float(nth_val) - float(n1_val)) / float(nth_val)
                    region = cell_id_to_name.get(cid, str(cid))
                    t_years = float((time_grid[istep] - time_grid[int(idx_shutdown)]) / SECONDS_PER_YEAR)
                    #print(
                    #    f"[DecayHeat] Included (n+1) due to tie in {region} (cell {cid}, mat {mat_id}) "
                    #    f"at step {istep} (t_rel={t_years:.3e} y): "
                    #    f"nth {nth_nuc}={float(nth_val):.3e} vs (n+1) {n1_nuc}={float(n1_val):.3e} "
                    #    f"(rel diff={rel_diff:.3%} < {similarity_threshold:.1%})"
                    #)

            for nuc, _ in selected:
                top_nucs_union.add(nuc)

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

        num_nucs = max(1, len(top_nucs))
        colors = plt.cm.tab20(np.linspace(0, 1, num_nucs))

        for i, nuc in enumerate(top_nucs):
            vals_plot = _mask(nuc_series[nuc])
            if not _has_positive_finite(vals_plot):
                continue

            ax.loglog(
                t_rel_plot,
                vals_plot,
                label=nuc,
                color=colors[i % len(colors)],
                marker=next(mcycle),
                linestyle=next(lscycle),
                markersize=4,
                linewidth=1.3,
            )

        if _has_positive_finite(others_plot):
            ax.loglog(t_rel_plot, others_plot, label="Others", linewidth=2.0)

        if _has_positive_finite(total_h_plot):
            ax.loglog(t_rel_plot, total_h_plot, color="black", linewidth=3.0, label="Total", zorder=10)

        region = cell_id_to_name.get(cid, str(cid))
        ax.set_xlabel("Time after irradiation [years]")
        ax.set_ylabel(f"Decay heat [{decayheat_units}]")
        ax.set_title(f"Decay heat (Top-{top_n} per timestep) — {region} (mat {mat_id})")

        _add_time_reference_lines(ax)
        _format_log_axes(ax)

        if total_h_plot.size and np.nanmax(total_h_plot) > 0.0:
            ymax = float(np.nanmax(total_h_plot))
            ymin = _ymin_from_topn_edge(global_min_topn_edge)
            ax.set_xlim(1e-10, 2e3)
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

    ax.set_xlim(1e-10, 2e3)
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
    base_out_dir: str | Path = "./depletion_results",
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

        #  Activity 
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
    
    #
    timer = Timer()

    timer.start("Import decay-chain")
    here = Path(__file__).resolve().parent
    # Prefer the reduced chain produced by depletion_model.py
    REDUCED_CHAIN = (here / "bluemira_chain.xml").resolve()
    # Fallback: full ENDF/B-VIII.0 chain
    FULL_CHAIN = (here.parent / "depletion_chain" / "chain_endfb80_sfr.xml").resolve()

    if REDUCED_CHAIN.exists():
        CHAIN_FILE = REDUCED_CHAIN
        print(f"[info] Using reduced depletion chain: {CHAIN_FILE}")
    elif FULL_CHAIN.exists():
        CHAIN_FILE = FULL_CHAIN
        print(f"[warn] Reduced chain not found; using full chain: {CHAIN_FILE}")
    else:
        raise FileNotFoundError(
            "No depletion chain found. Expected one of:\n"
            f"  - {REDUCED_CHAIN}\n"
            f"  - {FULL_CHAIN}"
        )

    openmc.config["chain_file"] = str(CHAIN_FILE)
    timer.stop("Import decay-chain")
    
    timer.start("import depletion resuts")
    results = openmc.deplete.Results("r2s/activation/depletion_results.h5")
    timer.stop("import depletion resuts")

    # To be updated:
    xcentroids_ob = (0.1, 1.1, 6.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 84.4, 156.0)
    xcentroids_ib = (0.1, 1.1, 6.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 73.8, 108.0)

    timer.start("Build neutronics model")
    from neutronics_model import build_breeder_chunks

    _cells_chunk = build_breeder_chunks(
    INPUT_JSON,
    default_equatorial_ob_key="OB_1_b6", 
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

    cell_ids_for_key   = _cells_chunk["cell_ids_for_key"]      
    radial_bins_for_key = _cells_chunk["radial_bins_for_key"]  

    OB_KEY = "OB_1_b6"
    IB_KEY = "IB_1_b4"

    CHUNKS = {
        OB_KEY: ob_by_key.get(OB_KEY, []),
        IB_KEY: ib_by_key.get(IB_KEY, []),
    }

    layer_tags_by_chunk: Dict[str, Sequence[str]] = {}
    xcentroids_by_chunk: Dict[str, Sequence[float]] = {}

    if CHUNKS[OB_KEY]:
        n_ob_layers = len(CHUNKS[OB_KEY]) - 3
        layer_tags_by_chunk[OB_KEY] = make_default_layer_tags(n_ob_layers, breeder_prefix="OB")
        if len(xcentroids_ob) == len(CHUNKS[OB_KEY]):
            xcentroids_by_chunk[OB_KEY] = xcentroids_ob

    if CHUNKS[IB_KEY]:
        n_ib_layers = len(CHUNKS[IB_KEY]) - 3
        layer_tags_by_chunk[IB_KEY] = make_default_layer_tags(n_ib_layers, breeder_prefix="IB")
        if len(xcentroids_ib) == len(CHUNKS[IB_KEY]):
            xcentroids_by_chunk[IB_KEY] = xcentroids_ib

    source_rates = results.get_source_rates()
    times = results.get_times()
    idx_shutdown = shutdown_index_from_source_rates(source_rates, n_steps=len(times))

    timer.stop("Build neutronics model")

    timer.start("Post-processing")
    run_chunk_postprocess(
        results,
        CHUNKS,
        base_out_dir="./depletion_results",
        xcentroids_by_chunk=xcentroids_by_chunk,
        layer_tags_by_chunk=layer_tags_by_chunk,
        idx_shutdown=int(idx_shutdown),
        activity_units="Bq/kg",       # "Bq" or "Bq/kg"
        decayheat_units="W/cm3",
        activity_top_n=15,
        decayheat_top_n=15,
        similarity_threshold=0.10,
        idx_to_plot=(0, 4, 8, 12, 16, 20),
    )
    timer.stop("Post-processing")
    timer.summary()

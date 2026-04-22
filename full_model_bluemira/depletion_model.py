#!/usr/bin/env python3
"""
depletion_model.py
==================
D1S shutdown dose rate + R2S activation depletion for the EU-DEMO blanket.

Imports from:
  inputs.py           : user configuration flags and paths
  geometry.py         : pydagmc_model, chunk helpers, neutron_source_rate
  neutronics_model.py : model (synced), all_cells

Feature toggles
---------------
  RUN_D1S        : bool  : run D1S shutdown dose calculation
  RUN_DEPLETION  : bool  : run R2S activation depletion
"""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import openmc
import openmc.data
import openmc.deplete
from openmc.deplete import d1s

import inputs as cfg
import geometry as geo

# ──────────────────────────────────────────────────────────────────────────────
# Feature toggles
# ──────────────────────────────────────────────────────────────────────────────
RUN_D1S       = True
RUN_DEPLETION = True

# ──────────────────────────────────────────────────────────────────────────────
# Import neutronics_model — gives us the fully synced model + all_cells
# This step will be removed once running neutronics+depletion together
# ──────────────────────────────────────────────────────────────────────────────
import neutronics_model as nm

model     = nm.model
all_cells = nm.all_cells

# Save the current tallies to restore later
orig_tallies = list(model.tallies)

# ──────────────────────────────────────────────────────────────────────────────
# Aliases from geometry
# ──────────────────────────────────────────────────────────────────────────────
pydagmc_model       = geo.pydagmc_model
neutron_source_rate = geo.neutron_source_rate
ob_by_key           = geo.ob_by_key
ib_by_key           = geo.ib_by_key
OB_CHUNK_SIZE       = geo.OB_CHUNK_SIZE
IB_CHUNK_SIZE       = geo.IB_CHUNK_SIZE

# ──────────────────────────────────────────────────────────────────────────────
# Output directories
# ──────────────────────────────────────────────────────────────────────────────
DEPLETION_RUN_DIR  = cfg.DEPLETION_RUN_DIR
DEPLETION_RUN_DIR.mkdir(parents=True, exist_ok=True)

D1S_DIR            = DEPLETION_RUN_DIR / "d1s"
D1S_DIR.mkdir(parents=True, exist_ok=True)

R2S_DIR            = DEPLETION_RUN_DIR / "r2s"
R2S_ACTIVATION_DIR = R2S_DIR / "activation"
R2S_ACTIVATION_DIR.mkdir(parents=True, exist_ok=True)

SDR_DIR     = D1S_DIR / "sdr"
SDR_DIR.mkdir(parents=True, exist_ok=True)

SDR_CSV_DIR = SDR_DIR / "csv"
SDR_CSV_DIR.mkdir(parents=True, exist_ok=True)

# ──────────────────────────────────────────────────────────────────────────────
# Timer
# ──────────────────────────────────────────────────────────────────────────────
class Timer:
    def __init__(self):
        self._start  = {}
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


timer = Timer()

# ──────────────────────────────────────────────────────────────────────────────
# Chunk helpers
# ──────────────────────────────────────────────────────────────────────────────
timer.start("Build chunk helpers")

OB_KEY      = cfg.ALBEDO_CHUNK_KEY
chunk_cells = geo.build_breeder_chunks(cfg.INPUT_JSON, default_chunk_key=OB_KEY)

cell_ids_all            = chunk_cells["cell_ids_all"]
cell_ids_selected_chunk = chunk_cells["selected_chunk_cell_ids"]
radial_bins_for_key     = chunk_cells["radial_bins_for_key"]

if OB_KEY not in ob_by_key:
    raise KeyError(
        f"{OB_KEY} not found in ob_by_key. "
        f"Available (sample): {list(ob_by_key)[:10]}"
    )

ob_1_b6_cells       = ob_by_key[OB_KEY]
xcentroids_ob, _, _ = radial_bins_for_key(OB_KEY)

if len(ob_1_b6_cells) != OB_CHUNK_SIZE:
    raise ValueError(
        f"len(ob_1_b6_cells)={len(ob_1_b6_cells)} "
        f"but OB_CHUNK_SIZE={OB_CHUNK_SIZE}."
    )

center_by_cell_ob6 = dict(zip(ob_1_b6_cells, xcentroids_ob))

timer.stop("Build chunk helpers")

# ──────────────────────────────────────────────────────────────────────────────
# Time grids / source rates
# ──────────────────────────────────────────────────────────────────────────────
s_to_h           = 3_600
y_to_s           = 24 * 365 * s_to_h
to_μSv           = 1e-6
to_mSv           = 1e-9
SECONDS_PER_YEAR = y_to_s

irradiation_time_y      = np.array([cfg.IRRADIATION_YEARS], dtype=float)
irradiation_intervals_s = (irradiation_time_y * y_to_s).tolist()

# Cooling milestone times since shutdown (absolute times)
_milestones_s = np.array([
    1.0,                # 1 second
    60.0,               # 1 minute
    3_600.0,            # 1 hour
    86_400.0,           # 1 day
    7 * 86_400.0,       # 1 week
    y_to_s / 12.0,      # 1 month
    1 * y_to_s,         # 1 year
    10 * y_to_s,        # 10 years
    100 * y_to_s,       # 100 years
    1_000 * y_to_s,     # 1000 years
], dtype=float)

# 2 interior log-spaced points between each consecutive milestone pair
_between_s = []
for _a, _b in zip(_milestones_s[:-1], _milestones_s[1:]):
    _between_s.extend(np.logspace(np.log10(_a), np.log10(_b), 4)[1:-1].tolist())

cooling_times_abs_s = np.unique(
    np.concatenate([_milestones_s, _between_s])
).astype(float)

cooling_intervals_s = np.diff(
    np.concatenate(([0.0], cooling_times_abs_s))
).tolist()

timesteps = irradiation_intervals_s + cooling_intervals_s

source_rates = (
    [cfg.CONSTANT_POWER_RATIO * neutron_source_rate] * len(irradiation_intervals_s)
    + [0.0] * len(cooling_intervals_s)
)

# Debug check: reconstruct absolute cooling times from intervals
_reconstructed_cooling_abs_s = np.cumsum(cooling_intervals_s)
print("[info] Cooling milestone check (years):")
for _t in [1, 60, 3600, 86400, 7*86400, y_to_s/12.0, y_to_s,
           10*y_to_s, 100*y_to_s, 1000*y_to_s]:
    _matches = np.isclose(
        _reconstructed_cooling_abs_s, _t,
        rtol=0.0, atol=1e-9 * max(1.0, _t)
    )
    print(f"  target = {_t / y_to_s:12.6e} y | found = {bool(np.any(_matches))}")

# ──────────────────────────────────────────────────────────────────────────────
# Cell volumes from all_cells (populated after DAGMC sync in neutronics_model)
# ──────────────────────────────────────────────────────────────────────────────
vol_by_cell = {
    int(cid): float(cell.volume)
    for cid, cell in all_cells.items()
    if cell.volume is not None and cell.volume > 0.0
}
print(f"[info] vol_by_cell populated for {len(vol_by_cell)} cells")

# ──────────────────────────────────────────────────────────────────────────────
# Identify plasma and VV port-fill cells (last two DAGMC volumes)
# ──────────────────────────────────────────────────────────────────────────────
vols = list(pydagmc_model.volumes)
if len(vols) < 2:
    raise RuntimeError(
        f"PyDAGMC model has only {len(vols)} volumes; expected >= 2."
    )

plasma_vol_id     = vols[-2].id
vvportfill_vol_id = vols[-1].id
print(f"[info] SDR cells: plasma={plasma_vol_id}, vv_port_fill={vvportfill_vol_id}")

for _cid in [plasma_vol_id, vvportfill_vol_id]:
    if _cid not in vol_by_cell:
        raise KeyError(
            f"Cell {_cid} has no volume. "
            f"Available (sample): {list(vol_by_cell)[:10]}"
        )
    print(f"[info] Cell {_cid} volume = {vol_by_cell[_cid]:.4e} cm³")

# ──────────────────────────────────────────────────────────────────────────────
# Plot helpers
# ──────────────────────────────────────────────────────────────────────────────
def _add_time_reference_lines(ax):
    """Vertical reference lines at human-readable time milestones (x in years).

    All x-positions are derived from SECONDS_PER_YEAR so the lines land
    exactly on data points from the matching time grid. Labels are in strictly
    ascending x order (seconds -> millennia).
    """
    refs = [
        (1.0           / SECONDS_PER_YEAR, "1 s"),
        (60.0          / SECONDS_PER_YEAR, "1 m"),
        (3_600.0       / SECONDS_PER_YEAR, "1 h"),
        (86_400.0      / SECONDS_PER_YEAR, "1 d"),
        (7 * 86_400.0  / SECONDS_PER_YEAR, "1 w"),
        (1.0 / 12.0,                       "1 mo"),
        (1.0,                              "1 y"),
        (10.0,                             "10 y"),
        (100.0,                            "100 y"),
        (1_000.0,                          "1 ky"),
    ]
    for x, txt in refs:
        ax.axvline(x, color="gray", linestyle="--", linewidth=1, alpha=0.6)
        ax.text(
            x, 0.98, txt,
            transform=ax.get_xaxis_transform(),
            rotation=90, va="top", ha="right",
            fontsize=8, alpha=0.8,
        )

def _format_log_axes(ax):
    """Apply log-log scale with clean minor ticks and grid."""
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="major", linewidth=0.8, alpha=0.35)
    ax.grid(True, which="minor", linewidth=0.5, alpha=0.15)
    ax.xaxis.set_minor_locator(
        mticker.LogLocator(base=10, subs=np.arange(2, 10) * 0.1)
    )
    ax.yaxis.set_minor_locator(
        mticker.LogLocator(base=10, subs=np.arange(2, 10) * 0.1)
    )
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())

# ──────────────────────────────────────────────────────────────────────────────
# D1S shutdown dose rate
# ──────────────────────────────────────────────────────────────────────────────
if RUN_D1S:
    print("\n==========================")
    print("D1S Shutdown Dose Rate")
    print("==========================")
    timer.start("Build D1S model")

    dose_cells_time = [plasma_vol_id, vvportfill_vol_id]
    dose_cells_prof = list(ob_1_b6_cells)
    dose_cells_all  = dose_cells_time + dose_cells_prof

    cell_filter2 = openmc.CellFilter(dose_cells_all)

    vol_by_cell_all: dict = {}
    for cid in dose_cells_all:
        v = vol_by_cell.get(cid)
        if v is None:
            raise KeyError(f"Missing volume for cell {cid}.")
        vol_by_cell_all[cid] = float(v)

    # Dose coefficients
    energies, pSv_cm2 = openmc.data.dose_coefficients(
        particle="photon", geometry="AP"
    )
    dose_filter   = openmc.EnergyFunctionFilter(
        energies, pSv_cm2, interpolation="cubic"
    )
    photon_filter = openmc.ParticleFilter("photon")

    dose_tally         = openmc.Tally(name="dose tally")
    dose_tally.filters = [dose_filter, photon_filter, cell_filter2]
    dose_tally.scores  = ["flux"]

    model.tallies                            = openmc.Tallies([dose_tally])
    model.settings.photon_transport          = True
    model.settings.use_decay_photons         = True

    nuclides = d1s.prepare_tallies(model)
    factors  = d1s.time_correction_factors(nuclides, timesteps, source_rates)

    timer.stop("Build D1S model")

    print("\n---------------------------")
    print("Performing D1S run")
    print("---------------------------")
    timer.start("D1S run")

    statepoint = model.run(
        cwd=str(D1S_DIR),
        output=False,
        threads=cfg.OPENMP_THREADS,
    )

    with openmc.StatePoint(statepoint) as sp:
        tally = sp.get_tally(name="dose tally")

    corrected_tallies = [
        d1s.apply_time_correction(tally, factors, i + 1)
        for i in range(len(timesteps))
    ]

    # Parse D1S results
    time_s            = []
    dose_time_by_cell = {plasma_vol_id: [], vvportfill_vol_id: []}
    profiles          = []
    cell_col          = "cell"  # default; overwritten inside loop

    for t_cool, ctally in zip(cooling_times_abs_s, corrected_tallies[1:]):
        d1s_df = ctally.get_pandas_dataframe()

        if "cell" in d1s_df.columns:
            cell_col = "cell"
        elif "cell_id" in d1s_df.columns:
            cell_col = "cell_id"
        else:
            raise KeyError("Could not find a cell column in tally dataframe.")

        d1s_df[cell_col]      = pd.to_numeric(d1s_df[cell_col]).astype(int)
        d1s_df["mean"]        = pd.to_numeric(d1s_df["mean"])
        d1s_df["cell_volume"] = d1s_df[cell_col].map(vol_by_cell_all)

        if d1s_df["cell_volume"].isna().any():
            missing = d1s_df.loc[
                d1s_df["cell_volume"].isna(), cell_col
            ].unique().tolist()
            raise KeyError(f"Missing volumes for cells: {missing}")

        d1s_df["μSv/h"] = (
            d1s_df["mean"] * (s_to_h * to_μSv) / d1s_df["cell_volume"]
        )
        d1s_df["mSv/h"] = (
            d1s_df["mean"] * (s_to_h * to_mSv) / d1s_df["cell_volume"]
        )

        # Time series: plasma & VV port-fill
        df_time = d1s_df[d1s_df[cell_col].isin(dose_cells_time)].copy()
        time_s.append(float(t_cool))

        for cid in dose_cells_time:
            cid = int(cid)
            sub = df_time.loc[df_time[cell_col] == cid, "μSv/h"]
            if len(sub) == 0:
                raise RuntimeError(
                    f"Missing μSv/h row for cid={cid} at t={t_cool:.3e}s"
                )
            dose_time_by_cell[cid].append(float(sub.iloc[0]))

        # Spatial profile: OB chunk
        df_prof            = d1s_df[d1s_df[cell_col].isin(ob_1_b6_cells)].copy()
        df_prof["centers"] = df_prof[cell_col].map(center_by_cell_ob6)

        if df_prof["centers"].isna().any():
            missing = df_prof.loc[
                df_prof["centers"].isna(), cell_col
            ].unique().tolist()
            raise KeyError(f"Missing centroids for {OB_KEY} cells: {missing}")

        df_prof = df_prof.sort_values("centers")
        profiles.append({"t_s": float(t_cool), "df": df_prof.copy()})

        print(
            f"Cooling {t_cool:.3e} s | "
            f"plasma={dose_time_by_cell[int(plasma_vol_id)][-1]:.3e} μSv/h | "
            f"vvpf={dose_time_by_cell[int(vvportfill_vol_id)][-1]:.3e} μSv/h | "
            f"{OB_KEY} rows={len(df_prof)}"
        )

    # Plot 1A: plasma time series
    t_rel_plot   = [t / SECONDS_PER_YEAR for t in time_s]
    plasma_doses = dose_time_by_cell[int(plasma_vol_id)]

    fig, ax = plt.subplots()
    ax.set_xscale("log")
    _add_time_reference_lines(ax)

    has_positive = any(v > 0 for v in plasma_doses)
    if has_positive:
        ax.set_yscale("log")

    ax.plot(t_rel_plot, plasma_doses)
    ax.axhline(0.1,   linestyle="--", color="k", linewidth=1.0)
    ax.axhline(10,    linestyle="--", color="k", linewidth=1.0)
    ax.axhline(10000, linestyle="--", color="k", linewidth=1.0)
    ax.text(1e-9, 0.1   * 1.5, "Natural background")
    ax.text(1e-9, 10    * 1.5, "Hands-on limit")
    ax.text(1e-9, 10000 * 1.5, "Remote recycling limit")
    ax.set_ylabel("Shutdown Dose [μSv/h]")
    ax.set_xlabel("Cooling Time [y]")

    if not has_positive:
        ax.set_title(
            f"Plasma SDR — zero values "
            f"({'expected for slab' if cfg.SIM_TYPE == 'slab' else 'check D1S result'})"
        )
        print(
            f"[warn] Plasma SDR has no positive values — "
            f"{'expected for slab SIM_TYPE' if cfg.SIM_TYPE == 'slab' else 'check D1S tally'}"
        )

    fig.savefig(SDR_DIR / "sdr_time_plasma.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Plot 1B: VV port-fill (optional) 
    show_VV = False
    if show_VV:
        fig, ax = plt.subplots()
        _add_time_reference_lines(ax)
        _format_log_axes(ax)
        ax.plot(
            t_rel_plot,
            dose_time_by_cell[int(vvportfill_vol_id)],
            label="VV port-fill (D1S)",
        )
        ax.set_ylabel("Shutdown Dose [μSv/h]")
        ax.set_xlabel("Cooling Time [y]")
        ax.legend()
        fig.savefig(SDR_DIR / "sdr_time_vvportfill.png",
                    dpi=300, bbox_inches="tight")
        plt.close(fig)

    # Plot 2: OB chunk spatial profiles
    if not profiles:
        raise RuntimeError(
            f"profiles is empty: {OB_KEY} rows were never captured."
        )

    ncurves = min(10, len(profiles))
    idxs    = np.linspace(0, len(profiles) - 1, ncurves, dtype=int)

    fig, ax = plt.subplots()
    for i in idxs:
        t_s_i = profiles[i]["t_s"]
        dfp   = profiles[i]["df"]
        ax.plot(
            dfp["centers"].to_numpy(float),
            dfp["μSv/h"].to_numpy(float),
            label=f"{t_s_i:.1e} s",
        )
    ax.set_yscale("log")
    ax.set_ylim(1e-6, 1e12)
    ax.set_ylabel("Shutdown Dose [μSv/h]")
    ax.set_xlabel("Radial Position [cm]")
    ax.set_title(f"D1S spatial profile: {OB_KEY}")
    ax.grid(True, which="both")
    ax.legend()
    fig.savefig(SDR_DIR / f"sdr_profile_{OB_KEY}.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Save SDR CSV outputs
    for prof in profiles:
        t_s    = prof["t_s"]
        t_y    = t_s / y_to_s
        df_out = prof["df"][[cell_col, "centers", "μSv/h", "mSv/h"]].copy()
        df_out = df_out.rename(columns={
            cell_col:  "cell_id",
            "centers": "radial_center_cm",
        })
        df_out["t_s"] = t_s
        df_out["t_y"] = t_y
        df_out.to_csv(
            SDR_CSV_DIR / f"sdr_profile_t{t_s:.4e}s.csv", index=False
        )

    pd.DataFrame({
        "t_s": time_s,
        "t_y": [t / y_to_s for t in time_s],
        f"plasma_cell{int(plasma_vol_id)}_uSvh":
            list(dose_time_by_cell[int(plasma_vol_id)]),
        f"vvpf_cell{int(vvportfill_vol_id)}_uSvh":
            list(dose_time_by_cell[int(vvportfill_vol_id)]),
    }).to_csv(SDR_DIR / "sdr_timeseries_plasma_vvpf.csv", index=False)

    summary_rows = []
    for prof in profiles:
        t_s = prof["t_s"]
        row = {"t_s": t_s, "t_y": t_s / y_to_s}
        for _, r in prof["df"].iterrows():
            cid = int(r[cell_col])
            cx  = float(r["centers"])
            row[f"cell_{cid}_x{cx:.1f}cm_uSvh"] = float(r["μSv/h"])
            row[f"cell_{cid}_x{cx:.1f}cm_mSvh"] = float(r["mSv/h"])
        summary_rows.append(row)

    pd.DataFrame(summary_rows).sort_values("t_s").reset_index(drop=True).to_csv(
        SDR_DIR / f"sdr_summary_{OB_KEY}.csv", index=False
    )

    print(f"[info] Wrote {len(profiles)} per-timestep SDR CSVs to {SDR_CSV_DIR}")
    print(f"[info] Wrote plasma/vvpf time series → "
          f"{SDR_DIR / 'sdr_timeseries_plasma_vvpf.csv'}")
    print(f"[info] Wrote OB summary → {SDR_DIR / f'sdr_summary_{OB_KEY}.csv'}")

    timer.stop("D1S run")

    # Restore model state after D1S
    model.tallies                    = orig_tallies
    model.settings.use_decay_photons = False

else:
    print("\n[info] RUN_D1S=False — skipping D1S shutdown dose calculation")
    model.tallies = orig_tallies

# ──────────────────────────────────────────────────────────────────────────────
# Depletion
# ──────────────────────────────────────────────────────────────────────────────
if RUN_DEPLETION:
    print("\n--------------------------------")
    print(f"Performing depletion ({cfg.SIM_TYPE})")
    print("--------------------------------")
    timer.start("Set up depletion model")

    openmc.deplete.pool.NUM_PROCESSES = cfg.DEPLETION_PROCESSES

    model.settings.use_decay_photons = False
    model.tallies                    = orig_tallies

    # Sync materials and determine geometry paths BEFORE differentiate_mats (safety)
    model.materials = openmc.Materials(
        list(model.geometry.get_all_materials().values())
    )
    model.geometry.determine_paths()

    # Step 1: mark depletable materials BEFORE differentiate_mats 
    # (Figure it out which materials will be depleted) -> matters more for slab
    if cfg.SIM_TYPE == "tokamak":
        _pre_source_cells = all_cells
        _pre_target_ids   = sorted(_pre_source_cells.keys())
        model.settings.photon_transport = False
    else:
        _pre_source_cells = model.geometry.get_all_cells()
        _all_json_ids     = set(int(cid) for cid in ob_by_key.get(OB_KEY, []))
        _present_ids      = set(int(c) for c in _pre_source_cells.keys())
        _pre_target_ids   = sorted(_all_json_ids & _present_ids)

    _n_marked = 0
    for _cid in _pre_target_ids:
        _cell = _pre_source_cells.get(int(_cid))
        if _cell is None:
            continue
        _mat = _cell.fill
        if not isinstance(_mat, openmc.Material):
            continue
        _vol = vol_by_cell.get(int(_cid))
        if _vol is None or float(_vol) <= 0.0:
            continue
        if _mat.volume is None:
            _mat.volume = float(_vol)
        _mat.depletable = True
        _n_marked += 1

    print(f"[depletion] Marked {_n_marked} materials as depletable")

    # Step 2: differentiate — one unique material per depletable cell 
    model.differentiate_mats("match cell", depletable_only=True)

    # Sync materials again after differentiation
    model.materials = openmc.Materials(
        list(model.geometry.get_all_materials().values())
    )
    model.geometry.determine_paths()

    # Step 3: collect the now-unique (cell, material) pairs
    dagmc_cell_ids: list = []
    deplete_cells:  list = []
    deplete_mats:   list = []

    if cfg.SIM_TYPE == "tokamak":
        _source_cells = all_cells
        target_ids    = sorted(_source_cells.keys())
        print(
            f"[depletion] tokamak: scanning "
            f"{len(target_ids)} cells for depletable materials"
        )
    else:
        _source_cells = model.geometry.get_all_cells()
        _all_json_ids = set(int(cid) for cid in ob_by_key.get(OB_KEY, []))
        _present_ids  = set(int(c) for c in _source_cells.keys())
        target_ids    = sorted(_all_json_ids & _present_ids)
        print(
            f"[depletion] slab: {len(target_ids)} cells from "
            f"{OB_KEY} present in wrapper geometry"
        )

    if not target_ids:
        raise RuntimeError(
            f"[depletion] No target cells found for SIM_TYPE={cfg.SIM_TYPE}. "
            f"Check geometry and chunk key {OB_KEY}."
        )

    for cid in target_ids:
        cell = _source_cells.get(int(cid))
        if cell is None:
            continue
        mat = cell.fill
        if not isinstance(mat, openmc.Material):
            continue
        if not mat.depletable:
            continue
        dagmc_cell_ids.append(int(cid))
        deplete_cells.append(cell)
        deplete_mats.append(mat)

    print(
        f"[depletion] Found {len(deplete_mats)} depletable materials "
        f"across {len(dagmc_cell_ids)} cells"
    )

    if not deplete_mats:
        raise RuntimeError(
            "deplete_mats is empty. Check material assignments in "
            "geometry.py / eudemo_materials.py."
        )

    # Mappings for post-processing
    cell_to_mat    = {cid: str(mat.id) for cid, mat in zip(dagmc_cell_ids, deplete_mats)}
    mat_to_cell    = {str(mat.id): cid for cid, mat in zip(dagmc_cell_ids, deplete_mats)}
    mat_id_to_name = {
        str(mat.id): (mat.name or f"material_{mat.id}")
        for mat in deplete_mats
    }

    # Build reduced chain
    initial_nuclides = model.geometry.get_all_nuclides()
    chain            = openmc.deplete.Chain.from_xml(str(cfg.OPENMC_CHAIN_FILE))
    print(
        f"[chain] Loaded {len(chain.nuclides)} nuclides "
        f"from: {cfg.OPENMC_CHAIN_FILE}"
    )

    reduced_chain = chain.reduce(initial_nuclides, level=cfg.REDUCED_CHAIN_LEVEL)
    print(f"[chain] Reduced to {len(reduced_chain.nuclides)} nuclides")

    if len(reduced_chain.nuclides) == 0:
        raise RuntimeError(
            "Reduced chain has 0 nuclides — nuclide names may not match chain. "
            f"Sample initial nuclides: {list(initial_nuclides)[:10]}"
        )

    bluemira_chain = DEPLETION_RUN_DIR / "bluemira_chain.xml"
    reduced_chain.export_to_xml(str(bluemira_chain))
    print(f"[info] Wrote reduced chain: {bluemira_chain}")

    openmc.config["chain_file"] = str(bluemira_chain)
    model.settings.depletion    = {"chain_file": str(bluemira_chain)}

    print(
        f"[depletion] Calling get_microxs_and_flux with "
        f"{len(deplete_mats)} materials ..."
    )

    # Generate fluxes and microscopic cross-sections
    fluxes, micros = openmc.deplete.get_microxs_and_flux(
        model,
        deplete_mats,
        chain_file=str(bluemira_chain),
        run_kwargs={
            "cwd":     str(DEPLETION_RUN_DIR),
            "output":  True,
            "threads": cfg.OPENMP_THREADS,
        },
    )

    # Define depletion operator (coupled/independent)
    operator = openmc.deplete.IndependentOperator(
        deplete_mats,
        fluxes,
        micros,
        chain_file=str(bluemira_chain),
        normalization_mode="source-rate",
    )

    # Set output directory
    operator.output_dir = str(R2S_ACTIVATION_DIR)

    # Set predictor integrator
    integrator = openmc.deplete.PredictorIntegrator(
        operator,
        timesteps,
        source_rates=source_rates,
    )

    timer.stop("Set up depletion model")

    timer.start("Depletion run")
    
    # Run depletion
    integrator.integrate()

    results = openmc.deplete.Results(
        str(R2S_ACTIVATION_DIR / "depletion_results.h5")
    )

    # Save cell_id and material_id mapping for depletion_post.py
    pd.DataFrame([
        {"cell_id": cid, "mat_id": int(mat.id), "material_name": mat.name}
        for cid, mat in zip(dagmc_cell_ids, deplete_mats)
    ]).to_csv(R2S_ACTIVATION_DIR / "cell_material_map.csv", index=False)

    print(f"[info] Written cell_material_map.csv to {R2S_ACTIVATION_DIR}")

    timer.stop("Depletion run")

else:
    print("\n[info] RUN_DEPLETION=False — skipping R2S activation depletion")

# ──────────────────────────────────────────────────────────────────────────────
# Timing summary
# ──────────────────────────────────────────────────────────────────────────────
timer.summary()

# ──────────────────────────────────────────────────────────────────────────────
# R2S decay-gamma run (placeholder)
# ──────────────────────────────────────────────────────────────────────────────
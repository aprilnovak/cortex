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

import json
import h5py
import os

import inputs as cfg
import geometry as geo

# Changes the random number seed to initiate the depletion run; used in conjunction
# with run_statistics_depletion.py
#model.settings.seed = int(os.getenv('OPENMC_RNG_SEED'))

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

# ──────────────────────────────────────────────────────────────────────────────
# Read neutronics run metadata (batches, particles) for matching transport runs
# ──────────────────────────────────────────────────────────────────────────────
_run_meta_path = cfg.NEUTRONICS_RESULTS_DIR / "run_meta.json"
with open(_run_meta_path) as _f:
    _run_meta = json.load(_f)

_neutronics_batches = int(_run_meta.get("batches_completed") or cfg.FIXED_BATCHES)
_neutronics_ppb     = int(_run_meta.get("particles_per_batch") or cfg.PARTICLES_PER_BATCH)
print(f"[Loading] Matching neutronics: {_neutronics_batches} batches × {_neutronics_ppb} particles")

# ──────────────────────────────────────────────────────────────────────────────
# Transport overrides for D1S and depletion microXS runs
# ──────────────────────────────────────────────────────────────────────────────
# Set to None to inherit from the neutronics run metadata above.
D1S_BATCHES          = None          # None → _neutronics_batches
D1S_PARTICLES        = None    # None → _neutronics_ppb
DEPLETION_BATCHES    = None          # None → _neutronics_batches
DEPLETION_PARTICLES  = None    # None → _neutronics_ppb

_d1s_batches    = D1S_BATCHES       if D1S_BATCHES       is not None else _neutronics_batches
_d1s_particles  = D1S_PARTICLES     if D1S_PARTICLES     is not None else _neutronics_ppb
_dep_batches    = DEPLETION_BATCHES  if DEPLETION_BATCHES  is not None else _neutronics_batches
_dep_particles  = DEPLETION_PARTICLES if DEPLETION_PARTICLES is not None else _neutronics_ppb

print(f"[D1S]       transport: {_d1s_batches} batches - {_d1s_particles} particles/batch")
print(f"[depletion] transport: {_dep_batches} batches - {_dep_particles} particles/batch")

# ──────────────────────────────────────────────────────────────────────────────
# Aliases from geometry
# ──────────────────────────────────────────────────────────────────────────────
pydagmc_model       = geo.pydagmc_model
neutron_source_rate = geo.neutron_source_rate
ob_by_key           = geo.ob_by_key
ib_by_key           = geo.ib_by_key
OB_CHUNK_SIZE       = geo.OB_CHUNK_SIZE
IB_CHUNK_SIZE       = geo.IB_CHUNK_SIZE

TRANSPORT_THREADS   = cfg.OPENMP_THREADS
DEPLETION_PROCESSES = cfg.DEPLETION_PROCESSES

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
# Neutron source rate
# ──────────────────────────────────────────────────────────────────────────────
if cfg.SIM_TYPE in ("slab", "fast_slab"):

    # Neutron ratio from surface source file (mirrors neutronics_model.py)
    if not cfg.SURFACE_SOURCE_FILE.is_file():
        raise FileNotFoundError(
            f"Surface source file not found: {cfg.SURFACE_SOURCE_FILE}"
        )
    with h5py.File(str(cfg.SURFACE_SOURCE_FILE), "r") as _f:
        _particles = _f["source_bank"]["particle"][:]
    _total                = len(_particles)
    _n_neutrons           = int((_particles == 0).sum())
    _neutron_ratio_source = _total / _n_neutrons if _n_neutrons > 0 else 1.0

    # J_in from armor current JSON (written by tokamak neutronics_post.py)
    if not cfg.ARMOR_CURRENT_JSON.is_file():
        raise FileNotFoundError(
            f"Armor current JSON not found: {cfg.ARMOR_CURRENT_JSON}\n"
            f"Run tokamak neutronics_post.py first to generate it."
        )
    with open(cfg.ARMOR_CURRENT_JSON) as _f:
        _records = json.load(_f)
    _J_in = float(_records[0]["J_in_mean"])

    neutron_source_rate: float = (
        _neutron_ratio_source
        * _J_in
        * (cfg.TOTAL_FUSION_POWER_W / cfg.NUMBER_OF_SECTORS)
        / (geo.ev_to_joule * cfg.EV_PER_FUSION)
    )
    print(f"[source] SIM_TYPE={cfg.SIM_TYPE}: neutron_source_rate (surface-scaled) = {neutron_source_rate:.6e} n/s")
    print(f"  neutron_ratio_source = {_neutron_ratio_source:.6f}")
    print(f"  J_in (ratio of neutrons from full tokamak source) = {_J_in:.6e}")

else:
    neutron_source_rate: float = geo.neutron_source_rate
    print(f"[source] SIM_TYPE={cfg.SIM_TYPE}: neutron_source_rate = {neutron_source_rate:.6e} n/s")

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

# 3 interior log-spaced points between each consecutive milestone pair
# except for the last three, where we add 8
_between_s = []
cnt = 0
for _a, _b in zip(_milestones_s[:-1], _milestones_s[1:]):
    if (cnt >= 6):
      _between_s.extend(np.logspace(np.log10(_a), np.log10(_b), 10)[1:-1].tolist())
    else:
      _between_s.extend(np.logspace(np.log10(_a), np.log10(_b), 5)[1:-1].tolist())
    cnt += 1

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

# ──────────────────────────────────────────────────────────────────────────────
# Surface-source helper for slab depletion
# ──────────────────────────────────────────────────────────────────────────────
def write_neutron_only_surface_source(src_path: Path, dst_path: Path) -> tuple[Path, float]:
    """
    Copy an OpenMC surface_source.h5 file but keep only neutron source particles.

    OpenMC particle IDs in source_bank:
      0 = neutron
      1 = photon
      2 = electron
      3 = positron

    Returns
    -------
    dst_path : Path
        Path to the neutron-only surface source.
    """
    src_path = Path(src_path)
    dst_path = Path(dst_path)

    if not src_path.exists():
        raise FileNotFoundError(f"Surface source file not found: {src_path}")

    dst_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(src_path, "r") as src:
        if "source_bank" not in src:
            raise KeyError(f"{src_path} does not contain /source_bank")

        bank = src["source_bank"][:]

        if "particle" not in bank.dtype.names:
            raise KeyError("source_bank does not contain a 'particle' field")

        neutron_bank = bank[bank["particle"] == 0]

        n_total = len(bank)
        n_neutrons = len(neutron_bank)

        if n_total == 0:
            raise RuntimeError(f"{src_path} contains an empty source_bank")

        if n_neutrons == 0:
            raise RuntimeError(
                f"No neutron particles found in {src_path}. "
                "Cannot use this source for neutron depletion."
            )

        with h5py.File(dst_path, "w") as dst:
            # Copy root attributes
            for key, val in src.attrs.items():
                dst.attrs[key] = val

            # Copy every dataset/group except source_bank.
            # Most surface source files only need source_bank + attributes,
            # but this keeps the copy safer if extra metadata exists.
            for key in src.keys():
                if key == "source_bank":
                    continue
                src.copy(key, dst)

            dst.create_dataset("source_bank", data=neutron_bank)

    print(
        "\n[surface source] Wrote neutron-only surface source"
        f"\n  input file          : {src_path}"
        f"\n  output file         : {dst_path}"
    )

    return dst_path

# ──────────────────────────────────────────────────────────────────────────────
# Cell volumes from all_cells (populated after DAGMC sync in neutronics_model)
# ──────────────────────────────────────────────────────────────────────────────
vol_by_cell = {
    int(cid): float(cell.volume)
    for cid, cell in all_cells.items()
    if cell.volume is not None and cell.volume > 0.0
}
print(f"[Volume] vol_by_cell populated for {len(vol_by_cell)} cells")

# ──────────────────────────────────────────────────────────────────────────────
# Identify plasma and VV port-fill cells (last two DAGMC volumes)
# ──────────────────────────────────────────────────────────────────────────────
if cfg.SIM_TYPE == "tokamak":
    vols = list(pydagmc_model.volumes)
    plasma_vol_id     = vols[-2].id
    vvportfill_vol_id = vols[-1].id
    #print(f"[info] SDR cells: plasma={plasma_vol_id}, vv_port_fill={vvportfill_vol_id}")

    for _cid in [plasma_vol_id, vvportfill_vol_id]:
        if _cid not in vol_by_cell:
            raise KeyError(
                f"Cell {_cid} has no volume. "
                f"Available (sample): {list(vol_by_cell)[:10]}"
            )
        print(f"[info] Cell {_cid} volume = {vol_by_cell[_cid]:.4e} cm³")
else:
    plasma_vol_id     = None
    vvportfill_vol_id = None
    #print(f"[D1S] SIM_TYPE={cfg.SIM_TYPE}: plasma/vvportfill cells not present — skipping SDR time-series")

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

    if cfg.SIM_TYPE == "tokamak":
        dose_cells_time = [plasma_vol_id, vvportfill_vol_id]
    else:
        dose_cells_time = []   # no plasma or port-fill in slab/fast_slab geometry

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

    model.tallies = openmc.Tallies([])
    print(f"[D1S] tallies after clearing: {len(model.tallies)}")

    dose_tally         = openmc.Tally(name="dose tally")
    dose_tally.filters = [dose_filter, photon_filter, cell_filter2]
    dose_tally.scores  = ["flux"]

    model.tallies                            = openmc.Tallies([dose_tally])
    model.settings.photon_transport          = True
    model.settings.use_decay_photons         = True

    print(f"[D1S] tallies for D1S: {len(model.tallies)}")

    # Set same number of batches to neutronics_model.py
    model.settings.trigger_active = False
    model.settings.batches        = _d1s_batches
    model.settings.particles      = _d1s_particles

    nuclides = d1s.prepare_tallies(model)
    factors  = d1s.time_correction_factors(nuclides, timesteps, source_rates)

    timer.stop("Build D1S model")

    print("\n---------------------------")
    print("Performing D1S run")
    print("---------------------------")
    timer.start("D1S run")

    statepoint = model.run(
        cwd=str(D1S_DIR),
        output=True,
        threads=TRANSPORT_THREADS,
    )

    with openmc.StatePoint(statepoint) as sp:
        tally = sp.get_tally(name="dose tally")

    # Apply time correction
    corrected_tallies = [
        d1s.apply_time_correction(tally, factors, i + 1)
        for i in range(len(timesteps))
    ]

    # Storage
    time_s            = []
    dose_time_by_cell = (
        {plasma_vol_id: [], vvportfill_vol_id: []}
        if cfg.SIM_TYPE == "tokamak"
        else {}
    )
    dose_std_dev_time_by_cell = (
        {plasma_vol_id: [], vvportfill_vol_id: []}
        if cfg.SIM_TYPE == "tokamak"
        else {}
    )
    profiles          = []
    cell_col          = "cell"  # default; overwritten inside loop

    for t_cool, ctally in zip(cooling_times_abs_s, corrected_tallies[1:]):
        d1s_df = ctally.get_pandas_dataframe()

        # Identify cell column
        if "cell" in d1s_df.columns:
            cell_col = "cell"
        elif "cell_id" in d1s_df.columns:
            cell_col = "cell_id"
        else:
            raise KeyError("Could not find a cell column in tally dataframe.")

        d1s_df[cell_col]      = pd.to_numeric(d1s_df[cell_col]).astype(int)
        d1s_df["mean"]        = pd.to_numeric(d1s_df["mean"])
        d1s_df["std. dev."]   = pd.to_numeric(d1s_df["std. dev."])

        # Map volumes
        d1s_df["cell_volume"] = d1s_df[cell_col].map(vol_by_cell_all)
        if d1s_df["cell_volume"].isna().any():
            missing = d1s_df.loc[
                d1s_df["cell_volume"].isna(), cell_col
            ].unique().tolist()
            raise KeyError(f"Missing volumes for cells: {missing}")

        # Convert to dose rate
        d1s_df["μSv/h"] = (
            d1s_df["mean"] * (s_to_h * to_μSv) / d1s_df["cell_volume"]
        )
        d1s_df["mSv/h"] = (
            d1s_df["mean"] * (s_to_h * to_mSv) / d1s_df["cell_volume"]
        )
        d1s_df["std. dev. μSv/h"] = (
            d1s_df["std. dev."] * (s_to_h * to_μSv) / d1s_df["cell_volume"]
        )

        # --------------------------------------------------
        # (A) Time series: plasma & VV port-fill (tokamak only)
        # --------------------------------------------------
        if cfg.SIM_TYPE == "tokamak":

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
                sub = df_time.loc[df_time[cell_col] == cid, "std. dev. μSv/h"]
                if len(sub) == 0:
                    raise RuntimeError(
                        f"Missing std. dev. μSv/h row for cid={cid} at t={t_cool:.3e}s"
                    )
                dose_std_dev_time_by_cell[cid].append(float(sub.iloc[0]))

        # --------------------------------------------------
        # (B) Spatial profile: OB chunk
        # --------------------------------------------------
        df_prof            = d1s_df[d1s_df[cell_col].isin(ob_1_b6_cells)].copy()
        df_prof["centers"] = df_prof[cell_col].map(center_by_cell_ob6)

        if df_prof["centers"].isna().any():
            missing = df_prof.loc[
                df_prof["centers"].isna(), cell_col
            ].unique().tolist()
            raise KeyError(f"Missing centroids for {OB_KEY} cells: {missing}")

        df_prof = df_prof.sort_values("centers")
        profiles.append({"t_s": float(t_cool), "df": df_prof.copy()})

        if cfg.SIM_TYPE == "tokamak":
            print(
                f"Cooling {t_cool:.3e} s | "
                f"plasma={dose_time_by_cell[int(plasma_vol_id)][-1]:.3e} μSv/h | "
                f"vvpf={dose_time_by_cell[int(vvportfill_vol_id)][-1]:.3e} μSv/h | "
                f"{OB_KEY} rows={len(df_prof)}"
            )
        else:
            print(f"Cooling {t_cool:.3e} s | {OB_KEY} rows={len(df_prof)}")

    # ------------------------------------------------------------------
    # Plot 1A: plasma time series
    # ------------------------------------------------------------------
    if cfg.SIM_TYPE == "tokamak":
        t_rel_plot   = [t / SECONDS_PER_YEAR for t in time_s]
        plasma_doses = dose_time_by_cell[int(plasma_vol_id)]
        hi = [plasma_doses[i] + dose_std_dev_time_by_cell[int(plasma_vol_id)][i] for i in range(len(plasma_doses))]
        lo = [plasma_doses[i] - dose_std_dev_time_by_cell[int(plasma_vol_id)][i] for i in range(len(plasma_doses))]

        fig, ax = plt.subplots()
        ax.set_xscale("log")
        _add_time_reference_lines(ax)

        has_positive = any(v > 0 for v in plasma_doses)
        if has_positive:
            ax.set_yscale("log")

        ax.plot(t_rel_plot, plasma_doses, marker='o', color='b', markersize=3)
        ax.fill_between(t_rel_plot, hi, lo,
                        alpha=0.25, lw=0, color='b')
        ax.axhline(0.1,   linestyle="--", color="k", linewidth=1.0)
        ax.axhline(10,    linestyle="--", color="k", linewidth=1.0)
        ax.axhline(10000, linestyle="--", color="k", linewidth=1.0)
        ax.text(5e-8, 0.1   * 1.5, "Natural background")
        ax.text(5e-8, 10    * 1.5, "Hands-on limit")
        ax.text(5e-8, 10000 * 1.5, "Remote recycling limit")
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

        # ------------------------------------------------------------------
        # Plot 1B: VV port-fill time series (OPTIONAL)
        # ------------------------------------------------------------------
        # TODO: SDR currently shows zero here. Checks are scheduled for the following weeks
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

        # ------------------------------------------------------------------
        # Per-timestep CSV: one row per cell (tokamak only — slab below)
        # ------------------------------------------------------------------
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

        # Time-series CSV: one row per cooling time, per point-of-interest cell
        pd.DataFrame({
            "t_s": time_s,
            "t_y": [t / y_to_s for t in time_s],
            f"plasma_cell{int(plasma_vol_id)}_uSvh":
                list(dose_time_by_cell[int(plasma_vol_id)]),
            f"plasma_cell{int(plasma_vol_id)}_uSvh_std_dev":
                list(dose_std_dev_time_by_cell[int(plasma_vol_id)]),
            f"vvpf_cell{int(vvportfill_vol_id)}_uSvh":
                list(dose_time_by_cell[int(vvportfill_vol_id)]),
        }).to_csv(SDR_DIR / "sdr_timeseries_plasma_vvpf.csv", index=False)

    # ------------------------------------------------------------------
    # Plot 2A: OB chunk spatial profiles (all sim types)
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # Plot 2B: Per-layer SDR time series for OB chunk (all sim types)
    # ------------------------------------------------------------------
    _n_breeder_layers = OB_CHUNK_SIZE - 3  # Armor + FW + VV = 3 fixed
    _layer_tags = (
        ["Armor", "First Wall"]
        + [f"Breeder layer {i}" for i in range(1, _n_breeder_layers + 1)]
        + ["VV"]
    )

    # Accumulate dose timeseries per cell, in cell order
    _dose_by_layer: dict[int, list[float]] = {int(cid): [] for cid in ob_1_b6_cells}
    for prof in profiles:
        dfp = prof["df"].set_index(cell_col)
        for cid in ob_1_b6_cells:
            cid_i = int(cid)
            val   = dfp.at[cid_i, "μSv/h"] if cid_i in dfp.index else float("nan")
            _dose_by_layer[cid_i].append(float(val))

    t_rel_plot_layers = [p["t_s"] / SECONDS_PER_YEAR for p in profiles]

    _cmap = plt.get_cmap("tab20", OB_CHUNK_SIZE)

    # 2B-i: Combined overview — all layers on one axes
    fig, ax = plt.subplots(figsize=(9, 5))
    _format_log_axes(ax)
    _add_time_reference_lines(ax)

    for idx, (cid, tag) in enumerate(zip(ob_1_b6_cells, _layer_tags)):
        doses = _dose_by_layer[int(cid)]
        if any(v > 0 and not np.isnan(v) for v in doses):
            ax.plot(
                t_rel_plot_layers,
                doses,
                label=tag,
                color=_cmap(idx),
                linewidth=1.4,
            )

    ax.axhline(0.1,   linestyle="--", color="k", linewidth=0.9, alpha=0.7)
    ax.axhline(10,    linestyle="--", color="k", linewidth=0.9, alpha=0.7)
    ax.axhline(10000, linestyle="--", color="k", linewidth=0.9, alpha=0.7)
    ax.text(t_rel_plot_layers[-1], 0.1   * 1.5, "Natural background", ha="right", fontsize=8)
    ax.text(t_rel_plot_layers[-1], 10    * 1.5, "Hands-on limit",      ha="right", fontsize=8)
    ax.text(t_rel_plot_layers[-1], 10000 * 1.5, "Remote recycling",    ha="right", fontsize=8)
    ax.set_ylim(bottom=0.01)
    ax.set_ylabel("Shutdown Dose [μSv/h]")
    ax.set_xlabel("Cooling Time [y]")
    ax.set_title(f"D1S SDR per OB layer — {OB_KEY}")
    ax.legend(fontsize=7, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(SDR_DIR / f"sdr_time_layers_{OB_KEY}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[info] Wrote combined per-layer SDR → {SDR_DIR / f'sdr_time_layers_{OB_KEY}.png'}")

    # 2B-ii: Individual plot per layer
    SDR_LAYERS_DIR = SDR_DIR / "layers"
    SDR_LAYERS_DIR.mkdir(parents=True, exist_ok=True)

    for idx, (cid, tag) in enumerate(zip(ob_1_b6_cells, _layer_tags)):
        doses = _dose_by_layer[int(cid)]
        has_positive = any(v > 0 and not np.isnan(v) for v in doses)

        fig, ax = plt.subplots(figsize=(8, 4.5))
        _format_log_axes(ax)
        _add_time_reference_lines(ax)

        if has_positive:
            ax.plot(
                t_rel_plot_layers,
                doses,
                color=_cmap(idx),
                linewidth=1.6,
            )
        else:
            ax.text(
                0.5, 0.5, "No positive SDR values",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=11, color="gray",
            )

        ax.axhline(0.1,   linestyle="--", color="k", linewidth=0.9, alpha=0.7)
        ax.axhline(10,    linestyle="--", color="k", linewidth=0.9, alpha=0.7)
        ax.axhline(10000, linestyle="--", color="k", linewidth=0.9, alpha=0.7)
        ax.text(t_rel_plot_layers[-1], 0.1   * 1.5, "Natural background", ha="right", fontsize=8)
        ax.text(t_rel_plot_layers[-1], 10    * 1.5, "Hands-on limit",      ha="right", fontsize=8)
        ax.text(t_rel_plot_layers[-1], 10000 * 1.5, "Remote recycling",    ha="right", fontsize=8)
        ax.set_ylim(bottom=0.01)
        ax.set_ylabel("Shutdown Dose [μSv/h]")
        ax.set_xlabel("Cooling Time [y]")
        ax.set_title(f"D1S SDR — {tag} ({OB_KEY})")
        fig.tight_layout()

        # filename: zero-padded layer index so files sort correctly
        _safe_tag = tag.replace(" ", "_").replace("/", "_")
        _fname    = f"sdr_layer_{idx:02d}_{_safe_tag}_{OB_KEY}.png"
        fig.savefig(SDR_LAYERS_DIR / _fname, dpi=300, bbox_inches="tight")
        plt.close(fig)

    print(
        f"[info] Wrote {OB_CHUNK_SIZE} individual layer SDR plots → {SDR_LAYERS_DIR}"
    )

    # ------------------------------------------------------------------
    # Per-timestep CSV for slab
    # ------------------------------------------------------------------
    if cfg.SIM_TYPE in ("slab", "fast_slab"):
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

    # Summary CSV: one row per cooling time, columns = OB_1_b6 cells
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
    if cfg.SIM_TYPE == "tokamak":
        print(f"[info] Wrote plasma/vvpf time series → "
              f"{SDR_DIR / 'sdr_timeseries_plasma_vvpf.csv'}")
    print(f"[D1S] Wrote OB summary → {SDR_DIR / f'sdr_summary_{OB_KEY}.csv'}")

    timer.stop("D1S run")

else:
    print("\n[D1S] RUN_D1S=False — skipping D1S shutdown dose calculation")

# ──────────────────────────────────────────────────────────────────────────────
# Depletion
# ──────────────────────────────────────────────────────────────────────────────
if RUN_DEPLETION:
    print("\n--------------------------------")
    print(f"Performing depletion ({cfg.SIM_TYPE})")
    print("--------------------------------")
    timer.start("Set up depletion model")

    import warnings
    try:
        from openmc.exceptions import IDWarning
    except ImportError:
        IDWarning = UserWarning

    # Control depletion multiprocessing
    openmc.deplete.pool.NUM_PROCESSES = DEPLETION_PROCESSES

    # -------------------------------------------------------------------------
    # Reset model state after D1S
    # -------------------------------------------------------------------------
    model.settings.use_decay_photons = False
    model.settings.photon_transport  = False
    
    # Ensure no surface source is written during depletion transport
    model.settings.surf_source_write = {}
    
    # Apply neutronics batch count to model settings before microXS transport
    # or other predefined number of batches and particles
    model.settings.trigger_active = False
    model.settings.batches   = _dep_batches
    model.settings.particles = _dep_particles

    # Clear tallies
    model.tallies = openmc.Tallies([])
    print(f"[depletion] tallies after clearing: {len(model.tallies)}")

    # -------------------------------------------------------------------------
    # Slab/fast_slab: neutron-only surface source + clear plasma source
    # -------------------------------------------------------------------------
    if cfg.SIM_TYPE in ("slab", "fast_slab"):
        print("\n[depletion] Preparing neutron-only slab surface source")

        neutron_only_surface_source = DEPLETION_RUN_DIR / "surface_source_neutrons_only.h5"
        write_neutron_only_surface_source(
            Path(cfg.SURFACE_SOURCE_FILE),
            neutron_only_surface_source,
        )
        model.settings.surf_source_read = {"path": str(neutron_only_surface_source)}
        model.settings.source           = []

    print("\n[depletion] Source-rate normalization")
    print(f"  SIM_TYPE                 = {cfg.SIM_TYPE}")
    print(f"  neutron_source_rate      = {neutron_source_rate:.8e} n/s")
    print(f"  cfg.CONSTANT_POWER_RATIO = {cfg.CONSTANT_POWER_RATIO:.8e}")
    print(f"  depletion source rate    = {source_rates[0]:.8e} n/s")

    # -------------------------------------------------------------------------
    # Collect target cells: gather nuclides for chain + mark depletable
    # -------------------------------------------------------------------------
    if cfg.SIM_TYPE == "tokamak":
        _pre_source_cells = all_cells
        _pre_target_ids   = sorted(_pre_source_cells.keys())
    else:
        _pre_source_cells = model.geometry.get_all_cells()
        _all_json_ids     = set(int(cid) for cid in ob_by_key.get(OB_KEY, []))
        _present_ids      = set(int(c) for c in _pre_source_cells.keys())
        _pre_target_ids   = sorted(_all_json_ids & _present_ids)

    initial_nuclides: set[str] = set()
    _n_marked = 0

    for _cid in _pre_target_ids:
        _cell = _pre_source_cells.get(int(_cid))
        if _cell is None:
            continue
        _mat = _cell.fill
        if not isinstance(_mat, openmc.Material):
            continue

        # Collect nuclides for chain reduction
        for nuc, _, _ in _mat.nuclides:
            initial_nuclides.add(nuc)

        # Mark depletable + assign volume
        _vol = vol_by_cell.get(int(_cid))
        if _vol is None or float(_vol) <= 0.0:
            continue
        if _mat.volume is None:
            _mat.volume = float(_vol)
        _mat.depletable = True
        _n_marked += 1

    print(f"[depletion] Marked {_n_marked} materials as depletable")
    print(f"[chain] Initial nuclides from deplete_mats: {len(initial_nuclides)}")

    if _n_marked == 0:
        raise RuntimeError(
            "[depletion] No materials were marked as depletable. "
            "Check target cell IDs, material assignments, and cell volumes."
        )

    # -------------------------------------------------------------------------
    # Build reduced depletion chain
    # -------------------------------------------------------------------------

    chain = openmc.deplete.Chain.from_xml(str(cfg.OPENMC_CHAIN_FILE))
    print(f"[chain] Loaded {len(chain.nuclides)} nuclides from: {cfg.OPENMC_CHAIN_FILE}")

    reduced_chain  = chain.reduce(initial_nuclides, level=cfg.REDUCED_CHAIN_LEVEL)
    bluemira_chain = DEPLETION_RUN_DIR / "bluemira_chain.xml"
    reduced_chain.export_to_xml(str(bluemira_chain))
    print(f"[chain] Reduced to {len(reduced_chain.nuclides)} nuclides")

    openmc.config["chain_file"] = str(bluemira_chain)
    model.settings.depletion    = {"chain_file": str(bluemira_chain)}

    # -------------------------------------------------------------------------
    # Sync materials and geometry paths before differentiate_mats
    # (tokamak only — for slab/fast_slab this collapses the universe hierarchy)
    # -------------------------------------------------------------------------
    
    if cfg.SIM_TYPE == "tokamak":
        model.materials = openmc.Materials(
            list(model.geometry.get_all_materials().values())
        )
        model.geometry.determine_paths()

    # -------------------------------------------------------------------------
    # Differentiate materials (tokamak only)
    # For slab/fast_slab: each cell already has a (defined) unique material — skip.
    # -------------------------------------------------------------------------

    if cfg.SIM_TYPE == "tokamak":
        model.differentiate_mats("match cell", depletable_only=True)
        model.materials = openmc.Materials(
            list(model.geometry.get_all_materials().values())
        )
        model.geometry.determine_paths()

    _microxs_model = model

    # -------------------------------------------------------------------------
    # Collect depletable materials
    # -------------------------------------------------------------------------
    deplete_mats = sorted(
        [m for m in _microxs_model.materials if m.depletable],
        key=lambda m: int(m.id),
    )

    print(f"[depletion] {len(deplete_mats)} depletable materials")

    if not deplete_mats:
        raise RuntimeError(
            "deplete_mats is empty. Check material assignments in "
            "geometry.py / eudemo_materials.py."
        )

    # -------------------------------------------------------------------------
    # Collect cell/material pairs for post-processing mappings
    # -------------------------------------------------------------------------
    if cfg.SIM_TYPE == "tokamak":
        _source_cells = all_cells
        target_ids    = sorted(_source_cells.keys())
        #print(f"[depletion] tokamak: scanning {len(target_ids)} cells")
    else:
        _source_cells = model.geometry.get_all_cells()
        _all_json_ids = set(int(cid) for cid in ob_by_key.get(OB_KEY, []))
        _present_ids  = set(int(c) for c in _source_cells.keys())
        target_ids    = sorted(_all_json_ids & _present_ids)
        #print(f"[depletion] slab: {len(target_ids)} cells from {OB_KEY}")

    if not target_ids:
        raise RuntimeError(
            f"[depletion] No target cells found for SIM_TYPE={cfg.SIM_TYPE}."
        )

    dagmc_cell_ids:       list[int]         = []
    deplete_cells:        list[openmc.Cell]  = []
    _orig_mat_id_to_cell: dict[int, int]     = {}

    if cfg.SIM_TYPE == "tokamak":
        for _cid in target_ids:
            _cell = _source_cells.get(int(_cid))
            if _cell is None:
                continue
            _mat = _cell.fill
            if not isinstance(_mat, openmc.Material) or not _mat.depletable:
                continue
            _orig_mat_id_to_cell[int(_mat.id)] = int(_cid)
            dagmc_cell_ids.append(int(_cid))
            deplete_cells.append(_cell)
    else:
        # slab/fast_slab: mat IDs unchanged
        for _cid in target_ids:
            _cell = _source_cells.get(int(_cid))
            if _cell is None:
                continue
            _mat = _cell.fill
            if not isinstance(_mat, openmc.Material) or not _mat.depletable:
                continue
            _orig_mat_id_to_cell[int(_mat.id)] = int(_cid)
            dagmc_cell_ids.append(int(_cid))
            deplete_cells.append(_cell)

    #print(f"[depletion] Found {len(dagmc_cell_ids)} depletable cells")

    if not dagmc_cell_ids:
        raise RuntimeError(
            "dagmc_cell_ids is empty. Check material assignments in "
            "geometry.py / eudemo_materials.py."
        )

    # -------------------------------------------------------------------------
    # Mappings for post-processing
    # -------------------------------------------------------------------------
    cell_to_mat = {
        _orig_mat_id_to_cell[int(mat.id)]: str(mat.id)
        for mat in deplete_mats if int(mat.id) in _orig_mat_id_to_cell
    }
    mat_to_cell = {
        str(mat.id): _orig_mat_id_to_cell[int(mat.id)]
        for mat in deplete_mats if int(mat.id) in _orig_mat_id_to_cell
    }
    mat_id_to_name = {
        str(mat.id): (mat.name or f"material_{mat.id}")
        for mat in deplete_mats
    }

    # -------------------------------------------------------------------------
    # Generate fluxes and microscopic cross sections
    # -------------------------------------------------------------------------
    fluxes, micros = openmc.deplete.get_microxs_and_flux(
        _microxs_model,
        deplete_mats,
        chain_file=str(bluemira_chain),
        run_kwargs={
            "cwd":     str(DEPLETION_RUN_DIR),
            "output":  True,
            "threads": TRANSPORT_THREADS,
        },
    )

    # -------------------------------------------------------------------------
    # Independent depletion operator + predictor integrator
    # -------------------------------------------------------------------------
    operator = openmc.deplete.IndependentOperator(
        deplete_mats,
        fluxes,
        micros,
        chain_file=str(bluemira_chain),
        normalization_mode="source-rate",
    )
    operator.output_dir = str(R2S_ACTIVATION_DIR)

    integrator = openmc.deplete.PredictorIntegrator(
        operator,
        timesteps,
        source_rates=source_rates,
    )

    timer.stop("Set up depletion model")

    # -------------------------------------------------------------------------
    # Run depletion
    # -------------------------------------------------------------------------
    timer.start("Depletion run")

    integrator.integrate()

    results = openmc.deplete.Results(
        str(R2S_ACTIVATION_DIR / "depletion_results.h5")
    )

    # -------------------------------------------------------------------------
    # Save cell/material map for depletion_post.py
    # -------------------------------------------------------------------------
    def _strip_cell_suffix(name: str, cid: int) -> str:
        suffix = f"_{cid}"
        return name[:-len(suffix)] if name and name.endswith(suffix) else name

    _cell_to_chunk: dict[int, str] = {}
    for _key, _cids in ob_by_key.items():
        for _cid in _cids:
            _cell_to_chunk[int(_cid)] = _key
    for _key, _cids in ib_by_key.items():
        for _cid in _cids:
            _cell_to_chunk[int(_cid)] = _key

    pd.DataFrame([
        {
            "cell_id":       _orig_mat_id_to_cell.get(int(mat.id), -1),
            "mat_id":        int(mat.id),
            "material_name": _strip_cell_suffix(
                mat.name or f"material_{mat.id}",
                _orig_mat_id_to_cell.get(int(mat.id), -1),
            ),
            "chunk_key":     _cell_to_chunk.get(
                _orig_mat_id_to_cell.get(int(mat.id), -1), ""
            ),
            "volume_cm3":    float(mat.volume) if mat.volume is not None else np.nan,
        }
        for mat in deplete_mats
    ]).sort_values("cell_id").to_csv(R2S_ACTIVATION_DIR / "cell_material_map.csv", index=False)

    print(f"[depletion] Written cell_material_map.csv to {R2S_ACTIVATION_DIR}")

    timer.stop("Depletion run")

else:
    print("\n[depletion] RUN_DEPLETION=False — skipping R2S activation depletion")

# ──────────────────────────────────────────────────────────────────────────────
# Timing summary
# ──────────────────────────────────────────────────────────────────────────────
timer.summary()

# ──────────────────────────────────────────────────────────────────────────────
# R2S decay-gamma run (placeholder)
# ──────────────────────────────────────────────────────────────────────────────

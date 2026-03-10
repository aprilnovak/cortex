#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple, Optional, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------
# User input
# ----------------------------
OB_KEY = "OB_1_b6" 

# ----------------------------
# Paths (script is in cortex/model_comparison/)
# ----------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
CORTEX_DIR = SCRIPT_DIR.parent

TOK_DIR = (CORTEX_DIR / "detailed_bluemira" / "neutronics_results" / OB_KEY).resolve()
SLB_DIR = (CORTEX_DIR / "slab"             / "neutronics_results" / OB_KEY).resolve()

TOK_DEP_DIR = (CORTEX_DIR / "detailed_bluemira" / "depletion_results" / OB_KEY).resolve()
SLB_DEP_DIR = (CORTEX_DIR / "slab"             / "depletion_results" / OB_KEY).resolve()

OUTDIR = (SCRIPT_DIR / "plots" / OB_KEY).resolve()
OUTDIR.mkdir(parents=True, exist_ok=True)

OUTDIR_NEUT_PROFILES = OUTDIR / "neutronics_profiles"
OUTDIR_NEUT_SPECTRA  = OUTDIR / "neutronics_spectra"
OUTDIR_DEP_ACTIVITY  = OUTDIR / "depletion_activity"
OUTDIR_DEP_DECAY     = OUTDIR / "depletion_decayheat"

OUTDIR_NEUT_PROFILES.mkdir(parents=True, exist_ok=True)
OUTDIR_NEUT_SPECTRA.mkdir(parents=True, exist_ok=True)
OUTDIR_DEP_ACTIVITY.mkdir(parents=True, exist_ok=True)
OUTDIR_DEP_DECAY.mkdir(parents=True, exist_ok=True)

# ----------------------------
# File discovery
# ----------------------------
def find_profile_csvs(folder: Path, chunk_key: str) -> Dict[str, Path]:
    """
    Map: quantity -> path
    expects: profile_<quantity>_<chunk_key>.csv
    """
    out: Dict[str, Path] = {}
    if not folder.is_dir():
        return out

    suffix = f"_{chunk_key}.csv"
    for p in sorted(folder.glob(f"profile_*_{chunk_key}.csv")):
        name = p.name
        if not name.startswith("profile_") or not name.endswith(suffix):
            continue
        quantity = name[len("profile_") : -len(suffix)]
        out[quantity] = p
    return out

def find_spectrum_csvs(folder: Path, chunk_key: str) -> Dict[str, Path]:
    """
    Map: particle -> path
    expects:
      neutron_spectrum_<chunk_key>.csv
      photon_spectrum_<chunk_key>.csv
    """
    out: Dict[str, Path] = {}
    if not folder.is_dir():
        return out

    for particle in ["neutron", "photon"]:
        p = folder / f"{particle}_spectrum_{chunk_key}.csv"
        if p.is_file():
            out[particle] = p

    return out

def find_depletion_csvs(folder: Path) -> Dict[str, Path]:
    """
    expects in folder:
      activity_all_cells.csv
      decayheat_all_cells.csv
    """
    out: Dict[str, Path] = {}
    if not folder.is_dir():
        return out

    p1 = folder / "activity_all_cells.csv"
    p2 = folder / "decayheat_all_cells.csv"

    if p1.is_file():
        out["activity"] = p1
    if p2.is_file():
        out["decayheat"] = p2

    return out

# ----------------------------
# Loading
# ----------------------------
def load_profile(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    need = {"x_left_cm", "x_right_cm", "mean", "std"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"{csv_path}: missing columns {sorted(missing)}")

    for c in ["x_left_cm", "x_right_cm", "mean", "std"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["x_left_cm", "x_right_cm", "mean", "std"]).copy()
    df = df.sort_values("x_left_cm").reset_index(drop=True)
    return df

def profile_to_step_arrays(df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_left = df["x_left_cm"].to_numpy(float)
    x_right = df["x_right_cm"].to_numpy(float)
    mean = df["mean"].to_numpy(float)
    std = df["std"].to_numpy(float)

    xedges = np.empty(len(mean) + 1, dtype=float)
    xedges[0] = x_left[0]
    xedges[1:] = x_right
    return xedges, mean, std

def load_spectrum(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)

    if "E_mid_eV" not in df.columns:
        raise ValueError(f"{csv_path}: missing column 'E_mid_eV'")

    df["E_mid_eV"] = pd.to_numeric(df["E_mid_eV"], errors="coerce")
    df = df.dropna(subset=["E_mid_eV"]).copy()
    df = df.sort_values("E_mid_eV").reset_index(drop=True)
    return df

def load_depletion_csv(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    if "Cooling_time_years" not in df.columns:
        raise ValueError(f"{csv_path}: missing column 'Cooling_time_years'")

    df["Cooling_time_years"] = pd.to_numeric(df["Cooling_time_years"], errors="coerce")
    df = df.dropna(subset=["Cooling_time_years"]).copy()

    for c in df.columns:
        if c == "Cooling_time_years":
            continue
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.sort_values("Cooling_time_years").reset_index(drop=True)
    return df

# ----------------------------
# Helpers
# ----------------------------
def sanitize_name(text: str) -> str:
    return str(text).replace(" ", "_").replace("/", "_")

def depletion_column_to_layer(col: str) -> str:
    if "__cell_" in col:
        return col.split("__cell_")[0]
    return col

# ----------------------------
# Plotting: neutronics profiles
# ----------------------------
def plot_neutronics_profile_tok_vs_slab(
    *,
    quantity: str,
    df_tok: pd.DataFrame,
    df_slb: pd.DataFrame,
    outpath: Path,
    ):
    xedges_t, y_t, s_t = profile_to_step_arrays(df_tok)
    xedges_s, y_s, s_s = profile_to_step_arrays(df_slb)

    t_low = np.maximum(y_t - s_t, 1e-30)
    t_up  = y_t + s_t
    s_low = np.maximum(y_s - s_s, 1e-30)
    s_up  = y_s + s_s

    units = None
    for df in (df_tok, df_slb):
        if "units" in df.columns:
            u = df["units"].dropna().astype(str)
            if len(u):
                units = u.iloc[0]
                break

    xcent_t = 0.5 * (xedges_t[:-1] + xedges_t[1:])
    xcent_s = 0.5 * (xedges_s[:-1] + xedges_s[1:])
    y_s_on_t = np.interp(xcent_t, xcent_s, y_s, left=np.nan, right=np.nan)

    with np.errstate(divide="ignore", invalid="ignore"):
        rel = (y_s_on_t - y_t) / y_t
    rel[~np.isfinite(rel)] = np.nan

    fig, ax = plt.subplots()

    ax.step(xedges_t, np.r_[y_t, y_t[-1]], where="post", label="tokamak")
    ax.fill_between(
        xedges_t,
        np.r_[t_low, t_low[-1]],
        np.r_[t_up,  t_up[-1]],
        step="post",
        alpha=0.30,
    )

    ax.step(xedges_s, np.r_[y_s, y_s[-1]], where="post", label="slab")
    ax.fill_between(
        xedges_s,
        np.r_[s_low, s_low[-1]],
        np.r_[s_up,  s_up[-1]],
        step="post",
        alpha=0.30,
    )

    ax.set_yscale("log")
    ylabel = quantity
    if units:
        ylabel += f" [{units}]"
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Radial Position [cm]")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)

    ax2 = ax.twinx()
    rel_percent = rel * 100.0

    if rel_percent.size and np.any(np.isfinite(rel_percent)):
        rel_step = np.r_[rel_percent, rel_percent[-1]]
        ax2.step(
            xedges_t,
            rel_step,
            where="post",
            color="red",
            linestyle="--",
            linewidth=1.8,
            label="(slab - tok)/tok [%]",
        )
        ax2.axhline(0.0, color="red", linestyle="--", linewidth=0.8)

    ax2.set_ylabel("Relative difference [%]", color="red")
    ax2.tick_params(axis="y", colors="red")

    rel_f = rel_percent[np.isfinite(rel_percent)]
    if rel_f.size > 0:
        rmin = float(rel_f.min())
        rmax = float(rel_f.max())
        if rmin == rmax:
            pad = 5.0 if rmin == 0 else 0.1 * abs(rmin)
            ax2.set_ylim(rmin - pad, rmax + pad)
        else:
            span = rmax - rmin
            ax2.set_ylim(rmin - 0.1 * span, rmax + 0.1 * span)

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=9)

    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)

# ----------------------------
# Plotting: neutronics spectra
# ----------------------------
def plot_neutronics_spectrum_tok_vs_slab_per_layer(
    *,
    particle: str,
    df_tok: pd.DataFrame,
    df_slb: pd.DataFrame,
    outdir: Path,
    chunk_key: str,
    pick_layers: Optional[List[str]] = None,
    rel_ylim: Tuple[float, float] = (-100.0, 100.0),
    rel_clip: float = 100.0,
    tok_floor: float = 1e-30,
    ):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    layers_tok = [c for c in df_tok.columns if c != "E_mid_eV"]
    layers_slb = [c for c in df_slb.columns if c != "E_mid_eV"]

    if pick_layers is None:
        layers = [L for L in layers_tok if L in layers_slb]
    else:
        layers = [L for L in pick_layers if (L in layers_tok and L in layers_slb)]

    if not layers:
        raise RuntimeError(f"No common layers to plot for spectrum {particle}.")

    E_tok = pd.to_numeric(df_tok["E_mid_eV"], errors="coerce").to_numpy(float)
    E_slb = pd.to_numeric(df_slb["E_mid_eV"], errors="coerce").to_numpy(float)

    for layer in layers:
        y_tok = pd.to_numeric(df_tok[layer], errors="coerce").to_numpy(float)
        y_slb = pd.to_numeric(df_slb[layer], errors="coerce").to_numpy(float)

        mt = np.isfinite(E_tok) & np.isfinite(y_tok) & (E_tok > 0.0) & (y_tok > 0.0)
        ms = np.isfinite(E_slb) & np.isfinite(y_slb) & (E_slb > 0.0) & (y_slb > 0.0)

        if not np.any(mt) or not np.any(ms):
            print(f"[warn] skipping {particle} {layer}: no positive valid data")
            continue

        y_slb_on_tok = np.interp(E_tok, E_slb[ms], y_slb[ms], left=np.nan, right=np.nan)

        with np.errstate(divide="ignore", invalid="ignore"):
            rel_percent = 100.0 * (y_slb_on_tok - y_tok) / y_tok

        mr = (
            np.isfinite(E_tok)
            & np.isfinite(y_tok)
            & np.isfinite(rel_percent)
            & (E_tok > 0.0)
            & (y_tok > tok_floor)
        )

        rel_plot = np.full_like(rel_percent, np.nan, dtype=float)
        rel_plot[mr] = np.clip(rel_percent[mr], -rel_clip, rel_clip)

        fig, ax = plt.subplots(figsize=(9, 6))

        ax.loglog(E_tok[mt], y_tok[mt], label="tokamak", linewidth=2)
        ax.loglog(E_slb[ms], y_slb[ms], linestyle="--", label="slab", linewidth=2)

        ax.set_xlabel("Energy [eV]")
        if particle.lower() == "neutron":
            ax.set_ylabel("Neutron flux per unit lethargy [1/cm$^2$/s]")
            ax.set_title(f"Neutron Spectrum — {layer}")
        else:
            ax.set_ylabel("Photon flux per unit lethargy [1/cm$^2$/s]")
            ax.set_title(f"Photon Spectrum — {layer}")

        ax.set_xlim([1, 100e6])
        ax.grid(True, which="both")

        ax2 = ax.twinx()
        if np.any(np.isfinite(rel_plot)):
            ax2.semilogx(
                E_tok[np.isfinite(rel_plot)],
                rel_plot[np.isfinite(rel_plot)],
                color="red",
                linestyle="--",
                linewidth=1.8,
                label="(slab - tok)/tok [%]",
            )

        ax2.axhline(0.0, color="red", linestyle="--", linewidth=0.8)
        ax2.set_ylim(rel_ylim)
        ax2.set_ylabel("Relative difference [%]", color="red")
        ax2.tick_params(axis="y", colors="red")

        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="best", fontsize=9)

        outpath = outdir / f"{particle}_{sanitize_name(layer)}_{chunk_key}.png"
        fig.savefig(outpath, dpi=300, bbox_inches="tight")
        plt.close(fig)

        n_clipped = int(np.sum(np.isfinite(rel_percent[mr]) & (np.abs(rel_percent[mr]) > rel_clip)))
        if n_clipped > 0:
            print(f"[info] {particle} {layer}: clipped {n_clipped} relative-difference points to plot range {rel_ylim}")

        print(f"[ok] wrote {outpath}")

# ----------------------------
# Plotting: depletion
# ----------------------------
def plot_depletion_tok_vs_slab(
    *,
    quantity: str,
    layer_label: str,
    time_tok_y: np.ndarray,
    values_tok: np.ndarray,
    time_slb_y: np.ndarray,
    values_slb: np.ndarray,
    outpath: Path,
    ):
    time_tok_y = np.asarray(time_tok_y, dtype=float)
    values_tok = np.asarray(values_tok, dtype=float)
    time_slb_y = np.asarray(time_slb_y, dtype=float)
    values_slb = np.asarray(values_slb, dtype=float)

    mt = np.isfinite(time_tok_y) & np.isfinite(values_tok) & (time_tok_y > 0.0) & (values_tok > 0.0)
    ms = np.isfinite(time_slb_y) & np.isfinite(values_slb) & (time_slb_y > 0.0) & (values_slb > 0.0)

    time_tok_y = time_tok_y[mt]
    values_tok = values_tok[mt]
    time_slb_y = time_slb_y[ms]
    values_slb = values_slb[ms]

    if len(time_tok_y) == 0 or len(time_slb_y) == 0:
        print(f"[warn] skipping depletion plot for {quantity} / {layer_label}: no positive finite data")
        return

    values_slb_on_tok = np.interp(time_tok_y, time_slb_y, values_slb, left=np.nan, right=np.nan)

    with np.errstate(divide="ignore", invalid="ignore"):
        rel = (values_slb_on_tok - values_tok) / values_tok
    rel[~np.isfinite(rel)] = np.nan
    rel_percent = rel * 100.0

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.loglog(time_tok_y, values_tok, label="tokamak", linewidth=2)
    ax.loglog(time_slb_y, values_slb, label="slab", linestyle="--", linewidth=2)

    if quantity.lower() == "activity":
        ax.set_ylabel("Activity")
        ax.set_title(f"Activity — {layer_label}")
    else:
        ax.set_ylabel("Decay heat")
        ax.set_title(f"Decay heat — {layer_label}")

    ax.set_xlabel("Cooling Time [years]")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)

    ax2 = ax.twinx()
    if np.any(np.isfinite(rel_percent)):
        ax2.semilogx(
            time_tok_y[np.isfinite(rel_percent)],
            rel_percent[np.isfinite(rel_percent)],
            color="red",
            linestyle="--",
            linewidth=1.8,
            label="(slab - tok)/tok [%]",
        )

    ax2.axhline(0.0, color="red", linestyle="--", linewidth=0.8)
    ax2.set_ylabel("Relative difference [%]", color="red")
    ax2.tick_params(axis="y", colors="red")

    rel_f = rel_percent[np.isfinite(rel_percent)]
    if rel_f.size > 0:
        rmin = float(rel_f.min())
        rmax = float(rel_f.max())
        if rmin == rmax:
            pad = 5.0 if rmin == 0 else 0.1 * abs(rmin)
            ax2.set_ylim(rmin - pad, rmax + pad)
        else:
            span = rmax - rmin
            ax2.set_ylim(rmin - 0.1 * span, rmax + 0.1 * span)

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="best", fontsize=9)

    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)

def make_depletion_plots_for_quantity(
    *,
    quantity: str,
    csv_tok: Path,
    csv_slb: Path,
    outdir: Path,
    ):
    df_tok = load_depletion_csv(csv_tok)
    df_slb = load_depletion_csv(csv_slb)

    tok_cols = [c for c in df_tok.columns if c != "Cooling_time_years"]
    slb_cols = [c for c in df_slb.columns if c != "Cooling_time_years"]

    tok_map = {depletion_column_to_layer(c): c for c in tok_cols}
    slb_map = {depletion_column_to_layer(c): c for c in slb_cols}

    layers = sorted(set(tok_map) & set(slb_map))
    if not layers:
        print(f"[warn] no common layers found in {quantity} csvs")
        return

    for layer in layers:
        col_t = tok_map[layer]
        col_s = slb_map[layer]

        outpath = outdir / f"{quantity}_{sanitize_name(layer)}_{OB_KEY}.png"
        plot_depletion_tok_vs_slab(
            quantity=quantity,
            layer_label=layer,
            time_tok_y=df_tok["Cooling_time_years"].to_numpy(float),
            values_tok=df_tok[col_t].to_numpy(float),
            time_slb_y=df_slb["Cooling_time_years"].to_numpy(float),
            values_slb=df_slb[col_s].to_numpy(float),
            outpath=outpath,
        )
        print(f"[ok] wrote {outpath}")

# ----------------------------
# Main
# ----------------------------
def main() -> int:
    print(f"[info] OB_KEY = {OB_KEY}")
    print(f"[info] tokamak neutronics: {TOK_DIR}")
    print(f"[info] slab    neutronics: {SLB_DIR}")
    print(f"[info] tokamak depletion:  {TOK_DEP_DIR}")
    print(f"[info] slab    depletion:  {SLB_DEP_DIR}")
    print(f"[info] out:                {OUTDIR}")

    # ---------------------------------
    # 1) Neutronics profiles
    # ---------------------------------
    tok_profiles = find_profile_csvs(TOK_DIR, OB_KEY)
    slb_profiles = find_profile_csvs(SLB_DIR, OB_KEY)

    if tok_profiles and slb_profiles:
        common_q = sorted(set(tok_profiles) & set(slb_profiles))
        for q in common_q:
            df_t = load_profile(tok_profiles[q])
            df_s = load_profile(slb_profiles[q])
            outpath = OUTDIR_NEUT_PROFILES / f"{q}_{OB_KEY}.png"
            plot_neutronics_profile_tok_vs_slab(
                quantity=q,
                df_tok=df_t,
                df_slb=df_s,
                outpath=outpath,
            )
            print(f"[ok] wrote {outpath}")
    else:
        print("[warn] profile csvs not found for one or both models")

    # ---------------------------------
    # 2) Neutronics spectra
    # ---------------------------------
    tok_spec = find_spectrum_csvs(TOK_DIR, OB_KEY)
    slb_spec = find_spectrum_csvs(SLB_DIR, OB_KEY)

    for particle in ["neutron", "photon"]:
        if particle in tok_spec and particle in slb_spec:
            df_t = load_spectrum(tok_spec[particle])
            df_s = load_spectrum(slb_spec[particle])

            plot_neutronics_spectrum_tok_vs_slab_per_layer(
                particle=particle,
                df_tok=df_t,
                df_slb=df_s,
                outdir=OUTDIR_NEUT_SPECTRA,
                chunk_key=OB_KEY,
                pick_layers=None,
                rel_ylim=(-100.0, 100.0),
                rel_clip=100.0,
                tok_floor=1e-30,
            )
        else:
            print(f"[warn] missing {particle} spectrum csv for one or both models")

    # ---------------------------------
    # 3) Depletion
    # ---------------------------------
    tok_dep = find_depletion_csvs(TOK_DEP_DIR)
    slb_dep = find_depletion_csvs(SLB_DEP_DIR)

    if "activity" in tok_dep and "activity" in slb_dep:
        make_depletion_plots_for_quantity(
            quantity="activity",
            csv_tok=tok_dep["activity"],
            csv_slb=slb_dep["activity"],
            outdir=OUTDIR_DEP_ACTIVITY,
        )
    else:
        print("[warn] activity_all_cells.csv not found for one or both models")

    if "decayheat" in tok_dep and "decayheat" in slb_dep:
        make_depletion_plots_for_quantity(
            quantity="decayheat",
            csv_tok=tok_dep["decayheat"],
            csv_slb=slb_dep["decayheat"],
            outdir=OUTDIR_DEP_DECAY,
        )
    else:
        print("[warn] decayheat_all_cells.csv not found for one or both models")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple, Optional, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OB_KEY = "OB_1_b6"
 
# Geometry case: one of "default", "HCLL", "WCLL", "HCPB"
GEOMETRY = "default" #"HCPB" #"default"

# Set to False if depletion results are not yet available
RUN_DEPLETION = False


SCRIPT_DIR = Path.cwd()
CORTEX_DIR = SCRIPT_DIR.parent
 
_GEOMETRY_DIRS: Dict[str, Tuple[str, str]] = {
    "default": ("detailed_bluemira",      "slab"),
    "HCLL":    ("detailed_bluemira_hcll", "slab_hcll"),
    "WCLL":    ("detailed_bluemira_wcll", "slab_wcll"),
    "HCPB":    ("detailed_bluemira_hcpb", "slab_hcpb"),
}
 
if GEOMETRY not in _GEOMETRY_DIRS:
    raise ValueError(f"Unknown GEOMETRY {GEOMETRY!r}. Choose from: {sorted(_GEOMETRY_DIRS)}")
 
_TOK_DIRNAME, _SLB_DIRNAME = _GEOMETRY_DIRS[GEOMETRY]
 
TOK_DIR     = (CORTEX_DIR / _TOK_DIRNAME / "neutronics_results" / OB_KEY)
SLB_DIR     = (CORTEX_DIR / _SLB_DIRNAME / "neutronics_results" / OB_KEY)
TOK_DEP_DIR = (CORTEX_DIR / _TOK_DIRNAME / "depletion_results"  / OB_KEY)
SLB_DEP_DIR = (CORTEX_DIR / _SLB_DIRNAME / "depletion_results"  / OB_KEY)
 
OUTDIR = (SCRIPT_DIR / "plots" / OB_KEY)
OUTDIR.mkdir(parents=True, exist_ok=True)

OUTDIR_NEUT_PROFILES = OUTDIR / "neutronics_profiles"
OUTDIR_NEUT_SPECTRA  = OUTDIR / "neutronics_spectra"
OUTDIR_DEP_ACTIVITY  = OUTDIR / "depletion_activity"
OUTDIR_DEP_DECAY     = OUTDIR / "depletion_decayheat"

for d in (OUTDIR_NEUT_PROFILES, OUTDIR_NEUT_SPECTRA, OUTDIR_DEP_ACTIVITY, OUTDIR_DEP_DECAY):
    d.mkdir(parents=True, exist_ok=True)


def find_profile_csvs(folder: Path, chunk_key: str) -> Dict[str, Path]:
    out: Dict[str, Path] = {}
    if not folder.is_dir():
        return out
    suffix = f"_{chunk_key}.csv"
    for p in sorted(folder.glob(f"profile_*_{chunk_key}.csv")):
        name = p.name
        if name.startswith("profile_") and name.endswith(suffix):
            out[name[len("profile_"):-len(suffix)]] = p
    return out


def find_spectrum_csvs(folder: Path, chunk_key: str) -> Dict[str, Path]:
    out: Dict[str, Path] = {}
    if not folder.is_dir():
        return out
    for particle in ("neutron", "photon"):
        p = folder / f"{particle}_spectrum_{chunk_key}.csv"
        if p.is_file():
            out[particle] = p
    return out


def find_depletion_csvs(folder: Path) -> Dict[str, Path]:
    out: Dict[str, Path] = {}
    if not folder.is_dir():
        return out
    for key, name in (("activity", "activity_all_cells.csv"), ("decayheat", "decayheat_all_cells.csv")):
        p = folder / name
        if p.is_file():
            out[key] = p
    return out


def load_profile(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    missing = {"x_left_cm", "x_right_cm", "mean", "std"} - set(df.columns)
    if missing:
        raise ValueError(f"{csv_path}: missing columns {sorted(missing)}")
    for c in ("x_left_cm", "x_right_cm", "mean", "std"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna(subset=["x_left_cm", "x_right_cm", "mean", "std"]).sort_values("x_left_cm").reset_index(drop=True)


def profile_to_step_arrays(df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_left  = df["x_left_cm"].to_numpy(float)
    x_right = df["x_right_cm"].to_numpy(float)
    mean    = df["mean"].to_numpy(float)
    std     = df["std"].to_numpy(float)
    xedges = np.empty(len(mean) + 1, dtype=float)
    xedges[0]  = x_left[0]
    xedges[1:] = x_right
    return xedges, mean, std


def load_spectrum(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    if "E_mid_eV" not in df.columns:
        raise ValueError(f"{csv_path}: missing column 'E_mid_eV'")
    df["E_mid_eV"] = pd.to_numeric(df["E_mid_eV"], errors="coerce")
    return df.dropna(subset=["E_mid_eV"]).sort_values("E_mid_eV").reset_index(drop=True)


def load_depletion_csv(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    if "Cooling_time_years" not in df.columns:
        raise ValueError(f"{csv_path}: missing column 'Cooling_time_years'")
    df["Cooling_time_years"] = pd.to_numeric(df["Cooling_time_years"], errors="coerce")
    df = df.dropna(subset=["Cooling_time_years"])
    value_cols = [c for c in df.columns if c != "Cooling_time_years"]
    df[value_cols] = df[value_cols].apply(pd.to_numeric, errors="coerce")
    return df.sort_values("Cooling_time_years").reset_index(drop=True)


def _sanitize(text: str) -> str:
    return str(text).replace(" ", "_").replace("/", "_")


def _col_to_layer(col: str) -> str:
    return col.split("__cell_")[0] if "__cell_" in col else col


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

    units = next(
        (df["units"].dropna().astype(str).iloc[0]
         for df in (df_tok, df_slb) if "units" in df.columns and len(df["units"].dropna())),
        None,
    )

    xcent_t = 0.5 * (xedges_t[:-1] + xedges_t[1:])
    xcent_s = 0.5 * (xedges_s[:-1] + xedges_s[1:])
    y_s_on_t = np.interp(xcent_t, xcent_s, y_s, left=np.nan, right=np.nan)

    with np.errstate(divide="ignore", invalid="ignore"):
        rel = (y_s_on_t - y_t) / y_t
    rel[~np.isfinite(rel)] = np.nan
    rel_percent = rel * 100.0

    fig, ax = plt.subplots()

    ax.step(xedges_t, np.r_[y_t, y_t[-1]], where="post", label="tokamak")
    ax.fill_between(xedges_t, np.r_[t_low, t_low[-1]], np.r_[t_up, t_up[-1]], step="post", alpha=0.30)
    ax.step(xedges_s, np.r_[y_s, y_s[-1]], where="post", label="slab")
    ax.fill_between(xedges_s, np.r_[s_low, s_low[-1]], np.r_[s_up, s_up[-1]], step="post", alpha=0.30)

    ax.set_yscale("log")
    ax.set_ylabel(f"{quantity} [{units}]" if units else quantity)
    ax.set_xlabel("Radial Position [cm]")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)

    ax2 = ax.twinx()
    if np.any(np.isfinite(rel_percent)):
        ax2.step(xedges_t, np.r_[rel_percent, rel_percent[-1]], where="post",
                 color="red", linestyle="--", linewidth=1.8, label="(slab - tok)/tok [%]")
    ax2.axhline(0.0, color="red", linestyle="--", linewidth=0.8)
    ax2.set_ylabel("Relative difference [%]", color="red")
    ax2.tick_params(axis="y", colors="red")

    rel_f = rel_percent[np.isfinite(rel_percent)]
    if rel_f.size > 0:
        rmin, rmax = float(rel_f.min()), float(rel_f.max())
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


def plot_neutronics_spectrum_tok_vs_slab_per_layer(
    *,
    particle: str,
    df_tok: pd.DataFrame,
    df_slb: pd.DataFrame,
    outdir: Path,
    pick_layers: Optional[List[str]] = None,
    rel_ylim: Tuple[float, float] = (-100.0, 100.0),
    rel_clip: float = 100.0,
    tok_floor: float = 1e-30,
):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Mean columns only (exclude _std columns)
    layers_tok = [c for c in df_tok.columns if c != "E_mid_eV" and not c.endswith("_std")]
    layers_slb = [c for c in df_slb.columns if c != "E_mid_eV" and not c.endswith("_std")]

    if pick_layers is None:
        layers = [L for L in layers_tok if L in layers_slb]
    else:
        layers = [L for L in pick_layers if L in layers_tok and L in layers_slb]

    if not layers:
        raise RuntimeError(f"No common layers to plot for spectrum {particle}.")

    E_tok = pd.to_numeric(df_tok["E_mid_eV"], errors="coerce").to_numpy(float)
    E_slb = pd.to_numeric(df_slb["E_mid_eV"], errors="coerce").to_numpy(float)

    for i, layer in enumerate(layers, start=1):
        layer_tag = f"layer_{i}"

        y_tok = pd.to_numeric(df_tok[layer], errors="coerce").to_numpy(float)
        y_slb = pd.to_numeric(df_slb[layer], errors="coerce").to_numpy(float)

        # Load std if available
        std_col = f"{layer}_std"
        s_tok = pd.to_numeric(df_tok[std_col], errors="coerce").to_numpy(float) if std_col in df_tok.columns else np.zeros_like(y_tok)
        s_slb = pd.to_numeric(df_slb[std_col], errors="coerce").to_numpy(float) if std_col in df_slb.columns else np.zeros_like(y_slb)

        mt = np.isfinite(E_tok) & np.isfinite(y_tok) & (E_tok > 0.0) & (y_tok > 0.0)
        ms = np.isfinite(E_slb) & np.isfinite(y_slb) & (E_slb > 0.0) & (y_slb > 0.0)

        if not np.any(mt) or not np.any(ms):
            print(f"No positive valid data for {particle} {layer_tag}, skipping")
            continue

        y_slb_on_tok = np.interp(E_tok, E_slb[ms], y_slb[ms], left=np.nan, right=np.nan)

        with np.errstate(divide="ignore", invalid="ignore"):
            rel_percent = 100.0 * (y_slb_on_tok - y_tok) / y_tok

        mr = (
            np.isfinite(E_tok) & np.isfinite(y_tok) & np.isfinite(rel_percent)
            & (E_tok > 0.0) & (y_tok > tok_floor)
        )
        rel_plot = np.full_like(rel_percent, np.nan, dtype=float)
        rel_plot[mr] = np.clip(rel_percent[mr], -rel_clip, rel_clip)

        fig, ax = plt.subplots(figsize=(9, 6))

        # Tokamak mean + std band
        ax.loglog(E_tok[mt], y_tok[mt], label="tokamak", linewidth=2)
        tok_low = np.maximum(y_tok - s_tok, tok_floor)
        tok_up  = y_tok + s_tok
        ax.fill_between(E_tok[mt], tok_low[mt], tok_up[mt], alpha=0.20)

        # Slab mean + std band
        ax.loglog(E_slb[ms], y_slb[ms], linestyle="--", label="slab", linewidth=2)
        slb_low = np.maximum(y_slb - s_slb, tok_floor)
        slb_up  = y_slb + s_slb
        ax.fill_between(E_slb[ms], slb_low[ms], slb_up[ms], alpha=0.20)

        ax.set_xlabel("Energy [eV]")
        ax.set_xlim([1, 100e6])
        ax.grid(True, which="both")

        if particle.lower() == "neutron":
            ax.set_ylabel("Neutron flux per unit lethargy [1/cm$^2$/s]")
            ax.set_title(f"Neutron Spectrum — {layer_tag}")
        else:
            ax.set_ylabel("Photon flux per unit lethargy [1/cm$^2$/s]")
            ax.set_title(f"Photon Spectrum — {layer_tag}")

        ax2 = ax.twinx()
        if np.any(np.isfinite(rel_plot)):
            ax2.semilogx(E_tok[np.isfinite(rel_plot)], rel_plot[np.isfinite(rel_plot)],
                         color="red", linestyle="--", linewidth=1.8, label="(slab - tok)/tok [%]")
        ax2.axhline(0.0, color="red", linestyle="--", linewidth=0.8)
        ax2.set_ylim(rel_ylim)
        ax2.set_ylabel("Relative difference [%]", color="red")
        ax2.tick_params(axis="y", colors="red")

        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="best", fontsize=9)

        outpath = outdir / f"{particle}_{layer_tag}.png"
        fig.savefig(outpath, dpi=300, bbox_inches="tight")
        plt.close(fig)


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

    time_tok_y, values_tok = time_tok_y[mt], values_tok[mt]
    time_slb_y, values_slb = time_slb_y[ms], values_slb[ms]

    if len(time_tok_y) == 0 or len(time_slb_y) == 0:
        print(f"No positive finite data for {quantity} / {layer_label}, skipping")
        return

    values_slb_on_tok = np.interp(time_tok_y, time_slb_y, values_slb, left=np.nan, right=np.nan)

    with np.errstate(divide="ignore", invalid="ignore"):
        rel_percent = 100.0 * (values_slb_on_tok - values_tok) / values_tok
    rel_percent[~np.isfinite(rel_percent)] = np.nan

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.loglog(time_tok_y, values_tok, label="tokamak", linewidth=2)
    ax.loglog(time_slb_y, values_slb, label="slab", linestyle="--", linewidth=2)
    ax.set_xlabel("Cooling Time [years]")
    ax.set_ylabel("Activity" if quantity.lower() == "activity" else "Decay heat")
    ax.set_title(f"{'Activity' if quantity.lower() == 'activity' else 'Decay heat'} — {layer_label}")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)

    ax2 = ax.twinx()
    if np.any(np.isfinite(rel_percent)):
        ax2.semilogx(time_tok_y[np.isfinite(rel_percent)], rel_percent[np.isfinite(rel_percent)],
                     color="red", linestyle="--", linewidth=1.8, label="(slab - tok)/tok [%]")
    ax2.axhline(0.0, color="red", linestyle="--", linewidth=0.8)
    ax2.set_ylabel("Relative difference [%]", color="red")
    ax2.tick_params(axis="y", colors="red")

    rel_f = rel_percent[np.isfinite(rel_percent)]
    if rel_f.size > 0:
        rmin, rmax = float(rel_f.min()), float(rel_f.max())
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

    tok_map = {_col_to_layer(c): c for c in df_tok.columns if c != "Cooling_time_years"}
    slb_map = {_col_to_layer(c): c for c in df_slb.columns if c != "Cooling_time_years"}

    layers = sorted(set(tok_map) & set(slb_map))
    if not layers:
        print(f"No common layers in {quantity} CSVs")
        return

    for i, layer in enumerate(layers, start=1):
        layer_tag = f"layer_{i}"
        plot_depletion_tok_vs_slab(
            quantity=quantity,
            layer_label=layer_tag,
            time_tok_y=df_tok["Cooling_time_years"].to_numpy(float),
            values_tok=df_tok[tok_map[layer]].to_numpy(float),
            time_slb_y=df_slb["Cooling_time_years"].to_numpy(float),
            values_slb=df_slb[slb_map[layer]].to_numpy(float),
            outpath=outdir / f"{quantity}_{layer_tag}.png",
        )


def main() -> int:
    tok_profiles = find_profile_csvs(TOK_DIR, OB_KEY)
    slb_profiles = find_profile_csvs(SLB_DIR, OB_KEY)

    if tok_profiles and slb_profiles:
        for q in sorted(set(tok_profiles) & set(slb_profiles)):
            plot_neutronics_profile_tok_vs_slab(
                quantity=q,
                df_tok=load_profile(tok_profiles[q]),
                df_slb=load_profile(slb_profiles[q]),
                outpath=OUTDIR_NEUT_PROFILES / f"{q}_{OB_KEY}.png",
            )
    else:
        print("Profile CSVs not found for one or both models")

    tok_spec = find_spectrum_csvs(TOK_DIR, OB_KEY)
    slb_spec = find_spectrum_csvs(SLB_DIR, OB_KEY)

    for particle in ("neutron", "photon"):
        if particle in tok_spec and particle in slb_spec:
            plot_neutronics_spectrum_tok_vs_slab_per_layer(
                particle=particle,
                df_tok=load_spectrum(tok_spec[particle]),
                df_slb=load_spectrum(slb_spec[particle]),
                outdir=OUTDIR_NEUT_SPECTRA,
                pick_layers=None,
                rel_ylim=(-100.0, 100.0),
                rel_clip=100.0,
                tok_floor=1e-30,
            )
        else:
            print(f"No {particle} spectrum CSV for one or both models")

    if RUN_DEPLETION:
            tok_dep = find_depletion_csvs(TOK_DEP_DIR)
            slb_dep = find_depletion_csvs(SLB_DEP_DIR)

            for quantity, outdir in (("activity", OUTDIR_DEP_ACTIVITY), ("decayheat", OUTDIR_DEP_DECAY)):
                if quantity in tok_dep and quantity in slb_dep:
                    make_depletion_plots_for_quantity(
                        quantity=quantity,
                        csv_tok=tok_dep[quantity],
                        csv_slb=slb_dep[quantity],
                        outdir=outdir,
                    )
                else:
                    print(f"{quantity}_all_cells.csv not found for one or both models")
    else:
        print("Depletion comparison skipped (RUN_DEPLETION = False)")


if __name__ == "__main__":
    raise SystemExit(main())
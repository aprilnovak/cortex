#!/usr/bin/env python3
"""
Post-process OpenMC statepoint and save plots/CSVs for equatorial model

Outputs:
  neutronics_results_dagmc_sector/
    summary.csv
    n_flux_spectrum.png
    p_flux_spectrum.png
    flux_total.png
    heating.png
    h_appm_fpy.png
    he_appm_fpy.png
    dpa_fpy.png

Edit:
  - STATEPOINT_FILE
  - MODEL_MODULE (the module where your tallies/cell metadata live)

  - Apply correction to:
      * H, He appm/y: weight the production tallies by f_struct_origin for each nuclide
      * DPA/y: weight damage-energy contribution by f_struct_origin for each nuclide

"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import openmc
import openmc.mgxs

# -----------------------------------------------------------------------------
# Import your model construction module
# -----------------------------------------------------------------------------
import neutronics_model as bm

# =============================================================================
# User settings
# =============================================================================
STATEPOINT_FILE = "statepoint.10.h5"

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "neutronics_results_dagmc_sector"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Fallback (only used if bm.neutron_source_rate missing)
TOTAL_POWER_W = 2.0e9
N_SECTORS = 16
EV_TO_J = 1.60218e-19
EV_FUSION = 17.6e6
S_IN_Y = 365.0 * 24.0 * 60.0 * 60.0
SURFACE_SOURCE_POWER_RATIO = 7.171062e-02

# MG energy structure
ENERGIES = np.asarray(openmc.mgxs.GROUP_STRUCTURES["CCFE-709"], dtype=float)
UNIT_LETHARGY = np.asarray(
    [math.log(ENERGIES[i + 1] / ENERGIES[i]) for i in range(len(ENERGIES) - 1)],
    dtype=float,
)

EPS = 1e-16

# =============================================================================
# Helpers
# =============================================================================
def require_attr(mod, name: str):
    if not hasattr(mod, name):
        raise RuntimeError(f"Missing `{name}` in neutronics_model.py module.")
    val = getattr(mod, name)
    if val is None:
        raise RuntimeError(f"`{name}` exists in neutronics_model.py but is None (did build() populate it?)")
    return val

def require_tally_id(sp: openmc.StatePoint, bm_tally_obj: openmc.Tally, desc: str) -> int:
    tid = int(bm_tally_obj.id)
    if tid not in sp.tallies:
        raise RuntimeError(
            f"[fatal] {desc}: bm id={tid} not found in statepoint. "
            "Your postprocessing and statepoint were produced by different tally sets."
        )
    return tid

def get_cell_bins_from_tally(tally: openmc.Tally) -> list[int]:
    cf = next(f for f in tally.filters if isinstance(f, openmc.CellFilter))
    return [int(x) for x in cf.bins]

def subset_by_cells(data: np.ndarray, cell_bins: list[int], desired_cells: list[int]) -> np.ndarray:
    idx = [cell_bins.index(int(cid)) for cid in desired_cells]
    return data[idx]

def element_from_nuclide(nuc: str) -> str:
    m = re.match(r"[A-Za-z]+", nuc)
    return m.group(0) if m else nuc

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
    """
    Σ_n T(score,n) * f_struct_origin(n), assume uncorrelated std.
    """
    mean_tot = 0.0
    var_tot = 0.0
    for nuc in nuclides:
        w = float(f_struct.get(str(nuc), 0.0))
        if w == 0.0:
            continue
        try:
            m, sd = tally_mean_std_for_nuclide(t, score=score, nuclide=str(nuc))
        except Exception:
            continue
        mean_tot += w * m
        var_tot += (w * sd) ** 2
    return mean_tot, math.sqrt(max(var_tot, 0.0))

def scaling_for_cells(cell_ids: list[int], neutron_source_rate: float) -> dict[int, float]:
    """
    scaling[cid] = neutron_source_rate / volume_cm3
    Uses bm.model geometry (lazy-build method ensures bm.model exists).
    """
    scaling: dict[int, float] = {}
    all_cells = bm.model.geometry.get_all_cells()
    for cid in cell_ids:
        cid = int(cid)
        vol = all_cells[cid].volume
        if vol is None or vol <= 0.0:
            raise ValueError(
                f"Cell {cid} has no valid volume (got {vol}).\n"
                "If this is a DAGMC volume-id mismatch, fix ID reservation/reset logic in neutronics_model.py."
            )
        scaling[cid] = SURFACE_SOURCE_POWER_RATIO * float(neutron_source_rate) / float(vol)
    return scaling

# =============================================================================
# DPA/Gas tally discovery (per-cell)
# =============================================================================
GAS_SCORES_EXPLICIT = {"H1-production", "He3-production", "He4-production"}
GAS_SCORES_REACTION = {"(n,Xp)", "(n,X3He)", "(n,Xa)"}

def _get_cell_id_from_tally(t: openmc.Tally) -> Optional[int]:
    cf = next((f for f in t.filters if isinstance(f, openmc.CellFilter)), None)
    if cf is None or not cf.bins:
        return None
    return int(cf.bins[0])

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

def save_struct_origin_csv(
    outdir: Path,
    filename: str,
    cell_ids: list[int],
    cell_struct_origin_frac: Dict[int, Dict[str, float]],
    *,
    min_f: float = 0.0,
) -> Path:
    """
    Save structural-origin fractions per cell/nuclide to CSV.
    """

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []

    for cid in cell_ids:
        cid = int(cid)
        frac_map = cell_struct_origin_frac.get(cid, {}) or {}

        if not frac_map:
            rows.append(
                dict(
                    cell_id=cid,
                    nuclide=None,
                    f_struct_origin=0.0,
                )
            )
            continue

        for nuc, f in frac_map.items():
            f = float(f)

            rows.append(
                dict(
                    cell_id=cid,
                    nuclide=str(nuc),
                    f_struct_origin=f,
                )
            )

    df = pd.DataFrame(rows)

    outpath = outdir / filename
    df.to_csv(outpath, index=False)

    return outpath

# =============================================================================
# X-bins: use provided centroids (and build edges from midpoints)
# =============================================================================
def build_x_edges_from_widths(
    widths: np.ndarray,
    *,
    x0: float = 0.0,
    ) -> np.ndarray:
    """
    Build contiguous bin edges from widths only:
      edges[0] = x0
      edges[i+1] = edges[i] + widths[i]
    """
    widths = np.asarray(widths, dtype=float).ravel()
    if widths.size == 0:
        raise ValueError("widths array is empty.")
    if np.any(widths <= 0.0):
        raise ValueError("All widths must be positive.")
    edges = np.empty(widths.size + 1, dtype=float)
    edges[0] = float(x0)
    edges[1:] = edges[0] + np.cumsum(widths)
    return edges

def centroids_from_edges(edges: np.ndarray) -> np.ndarray:
    edges = np.asarray(edges, dtype=float).ravel()
    if edges.size < 2:
        raise ValueError("edges must have length >= 2.")
    return 0.5 * (edges[:-1] + edges[1:])

# =============================================================================
# Main processor
# =============================================================================
def process(sp: openmc.StatePoint):

    # cell ids / centroids
    if hasattr(bm, "cell_ids") and getattr(bm, "cell_ids") is not None:
        cell_ids = [int(c) for c in getattr(bm, "cell_ids")]
    else:
        # fallback: hard-code here if your model doesn't export cell_ids yet
        cell_ids = [56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66]

    if hasattr(bm, "xcentroids_ob") and getattr(bm, "xcentroids_ob") is not None:
        xcent = np.asarray(getattr(bm, "xcentroids_ob"), dtype=float).ravel()
    else:
        xcent = np.asarray((0.1, 1.1, 6.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 84.4, 156.0), dtype=float)

    if len(xcent) != len(cell_ids):
        raise RuntimeError(
            f"xcentroids length ({len(xcent)}) != cell_ids length ({len(cell_ids)}). "
            "Export bm.xcentroids_ob aligned to bm.cell_ids, or update the fallback."
        )


    x0 = float(xcent[0] - 0.5 * bm.widths[0])
    xedges = build_x_edges_from_widths(bm.widths, x0=x0)
    xcent_effective = centroids_from_edges(xedges)

    # labels from bm.model geometry
    all_cells = bm.model.geometry.get_all_cells()
    labels = [str(all_cells[c].name) for c in cell_ids]

    # scaling (source_rate / volume)
    scaling = scaling_for_cells(cell_ids, bm.neutron_source_rate)

    # structural-origin maps
    #cell_struct_origin_frac: Dict[int, Dict[str, float]] = require_attr(bm, "cell_struct_origin_frac")
    cell_struct_origin_frac = bm.cell_struct_origin_frac
    cell_total_atoms_struct= bm.cell_total_atoms_struct#: Dict[int, float] = require_attr(bm, "cell_total_atoms_struct")

      # dict[cid] -> dict[nuc] -> f_origin
    save_struct_origin_csv(
        outdir=RESULTS_DIR,
        filename="struct_origin_debug.csv",
        cell_ids=cell_ids,   # or any list of cell_ids
        cell_struct_origin_frac=bm.cell_struct_origin_frac,
        min_f=1e-6              # optional filter
    )

    # tallies (by ID)
    t_flux_id = require_tally_id(sp, bm.flux_tally, "flux spectrum tally")
    t_flux = sp.get_tally(id=t_flux_id)
    
    t_flux_tot_id = require_tally_id(sp, bm.flux_tally_total, "total flux tally")
    t_flux_tot = sp.get_tally(id=t_flux_tot_id)
    
    t_heat_id = require_tally_id(sp, bm.heating_tally, "heating tally")
    t_heat = sp.get_tally(id=t_heat_id)

    # per-cell DPA/gas tallies found in statepoint
    dpa_gas_map = build_dpa_gas_map(sp)

    # ----------------------------
    # Flux spectrum (neutron+photon)
    # ----------------------------
    cell_bins_flux = get_cell_bins_from_tally(t_flux)
    flux_mean = t_flux.get_reshaped_data(value="mean")

    # shape: (nCellBins, nParticles(2), nGroups)
    neutron_flux = flux_mean[:, 0, :]
    photon_flux = flux_mean[:, 1, :]

    neutron_flux_chunk = subset_by_cells(neutron_flux, cell_bins_flux, cell_ids)
    photon_flux_chunk = subset_by_cells(photon_flux, cell_bins_flux, cell_ids)

    # Neutron spectrum plot
    plt.figure()
    for i, cid in enumerate(cell_ids):
        spec = neutron_flux_chunk[i].ravel() * scaling[int(cid)] / UNIT_LETHARGY
        plt.loglog(ENERGIES[:-1], spec, label=f"{labels[i]} (x={xcent[i]:.2f} cm)")
    plt.grid(True, which="both")
    plt.ylabel("Neutron flux per unit lethargy [1/cm$^2$/s]")
    plt.xlabel("Energy [eV]")
    plt.xlim([1.0, 1.0e8])
    plt.legend(fontsize=8, ncol=2)
    plt.savefig(RESULTS_DIR / "n_flux_spectrum.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Photon spectrum plot
    plt.figure()
    for i, cid in enumerate(cell_ids):
        spec = photon_flux_chunk[i].ravel() * scaling[int(cid)] / UNIT_LETHARGY
        plt.loglog(ENERGIES[:-1], spec, label=f"{labels[i]} (x={xcent[i]:.2f} cm)")
    plt.grid(True, which="both")
    plt.ylabel("Photon flux per unit lethargy [1/cm$^2$/s]")
    plt.xlabel("Energy [eV]")
    plt.xlim([1.0, 1.0e8])
    plt.legend(fontsize=8, ncol=2)
    plt.savefig(RESULTS_DIR / "p_flux_spectrum.png", dpi=300, bbox_inches="tight")
    plt.close()

    # ----------------------------
    # Total flux vs x (neutron + photon)
    # ----------------------------
    cell_bins_tot = get_cell_bins_from_tally(t_flux_tot)
    tot_mean = t_flux_tot.get_reshaped_data(value="mean").squeeze()
    tot_std = t_flux_tot.get_reshaped_data(value="std_dev").squeeze()
    # tot arrays are (nCellBins, nParticles)

    tot_mean_chunk = subset_by_cells(tot_mean, cell_bins_tot, cell_ids)
    tot_std_chunk = subset_by_cells(tot_std, cell_bins_tot, cell_ids)

    total_neut = np.array([tot_mean_chunk[i, 0] * scaling[int(cid)] for i, cid in enumerate(cell_ids)], float)
    total_phot = np.array([tot_mean_chunk[i, 1] * scaling[int(cid)] for i, cid in enumerate(cell_ids)], float)
    total_neut_std = np.array([tot_std_chunk[i, 0] * scaling[int(cid)] for i, cid in enumerate(cell_ids)], float)
    total_phot_std = np.array([tot_std_chunk[i, 1] * scaling[int(cid)] for i, cid in enumerate(cell_ids)], float)

    fig, ax = plt.subplots()
    ax.set_yscale("log")

    n_lo = np.maximum(total_neut - total_neut_std, EPS)
    n_hi = np.maximum(total_neut + total_neut_std, EPS)
    p_lo = np.maximum(total_phot - total_phot_std, EPS)
    p_hi = np.maximum(total_phot + total_phot_std, EPS)

    ax.step(xedges, np.r_[total_neut, total_neut[-1]], where="post", label="neutron")
    ax.fill_between(xedges, np.r_[n_lo, n_lo[-1]], np.r_[n_hi, n_hi[-1]], step="post", alpha=0.15)

    ax.step(xedges, np.r_[total_phot, total_phot[-1]], where="post", label="photon")
    ax.fill_between(xedges, np.r_[p_lo, p_lo[-1]], np.r_[p_hi, p_hi[-1]], step="post", alpha=0.15)

    ax.set_xlabel("Radial position x [cm]")
    ax.set_ylabel("Total flux [1/cm$^2$/s]")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.legend()
    fig.savefig(RESULTS_DIR / "flux_total.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ----------------------------
    # Heating vs x  (W/cm^3)
    # ----------------------------
    cell_bins_heat = get_cell_bins_from_tally(t_heat)
    heat_mean = t_heat.get_values(value="mean").ravel()
    heat_std = t_heat.get_values(value="std_dev").ravel()

    heat_mean_chunk = subset_by_cells(heat_mean, cell_bins_heat, cell_ids)
    heat_std_chunk = subset_by_cells(heat_std, cell_bins_heat, cell_ids)

    # (eV/source) * (source/s / cm^3) -> eV/(cm^3*s) -> W/cm^3
    heat_w = np.array([heat_mean_chunk[i] * scaling[int(cid)] * EV_TO_J for i, cid in enumerate(cell_ids)], float)
    heat_w_std = np.array([heat_std_chunk[i] * scaling[int(cid)] * EV_TO_J for i, cid in enumerate(cell_ids)], float)

    fig, ax = plt.subplots()
    ax.set_yscale("log")

    h_lo = np.maximum(heat_w - heat_w_std, EPS)
    h_hi = np.maximum(heat_w + heat_w_std, EPS)

    ax.step(xedges, np.r_[heat_w, heat_w[-1]], where="post")
    ax.fill_between(xedges, np.r_[h_lo, h_lo[-1]], np.r_[h_hi, h_hi[-1]], step="post", alpha=0.25)

    ax.set_xlabel("Radial position x [cm]")
    ax.set_ylabel("Heating [W/cm$^3$]")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    fig.savefig(RESULTS_DIR / "heating.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # =======
    # source scaling
    # ========
    source_per_y = SURFACE_SOURCE_POWER_RATIO * float(bm.neutron_source_rate) * float(S_IN_Y) # n per fpy 

    # ----------------------------
    # H / He appm/fpy (struct-origin corrected)
    # ----------------------------
    h_appm_y = np.zeros(len(cell_ids), dtype=float)
    h_appm_y_std = np.zeros(len(cell_ids), dtype=float)
    he_appm_y = np.zeros(len(cell_ids), float)
    he_appm_y_std = np.zeros(len(cell_ids), float)

    for i, cid in enumerate(cell_ids):
        cid = int(cid)
        t = dpa_gas_map.get(cid)
        if t is None:
            continue

        scores = {str(s) for s in (t.scores or [])}
        f_struct = cell_struct_origin_frac.get(cid, {}) or {}
        nuclides = list(t.nuclides or [])

        denom = float(cell_total_atoms_struct.get(cid, 0.0))
        if denom <= 0.0:
            continue

        
        source_scale = source_per_y / denom * 1e6  # appm/fpy

        # H
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

        # Apply here the structural_nuclide_fraction
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
            h_appm_y[i] = h_mean_corr * source_scale
            h_appm_y_std[i] = h_std_corr * source_scale

        # He
        if GAS_SCORES_EXPLICIT.issubset(scores):
            he3_score, he4_score = "He3-production", "He4-production"
        elif GAS_SCORES_REACTION.issubset(scores):
            he3_score, he4_score = "(n,X3He)", "(n,Xa)"
        else:
            he3_score = he4_score = None

        if he3_score is not None:
            m3, sd3 = corrected_sum_mean_std_getvalues(t, score=he3_score, nuclides=nuclides, f_struct=f_struct)
            m4, sd4 = corrected_sum_mean_std_getvalues(t, score=he4_score, nuclides=nuclides, f_struct=f_struct)
            he_appm_y[i] = (m3 + m4) * source_scale
            he_appm_y_std[i] = math.sqrt(sd3**2 + sd4**2) * source_scale

    # Plot H
    fig, ax = plt.subplots()
    ax.set_yscale("log")
    lo = np.maximum(h_appm_y - h_appm_y_std, 1e-30)
    hi = np.maximum(h_appm_y + h_appm_y_std, 1e-30)
    ax.step(xedges, np.r_[h_appm_y, h_appm_y[-1]], where="post")
    ax.fill_between(xedges, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]], step="post", alpha=0.25)
    ax.set_xlabel("Radial position x [cm]")
    ax.set_ylabel("H [appm/fpy] (struct-origin)")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    fig.savefig(RESULTS_DIR / "h_appm_fpy.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Plot He
    fig, ax = plt.subplots()
    ax.set_yscale("log")
    lo = np.maximum(he_appm_y - he_appm_y_std, 1e-30)
    hi = np.maximum(he_appm_y + he_appm_y_std, 1e-30)
    ax.step(xedges, np.r_[he_appm_y, he_appm_y[-1]], where="post")
    ax.fill_between(xedges, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]], step="post", alpha=0.25)
    ax.set_xlabel("Radial position x [cm]")
    ax.set_ylabel("He [appm/fpy] (struct-origin)")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    fig.savefig(RESULTS_DIR / "he_appm_fpy.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ----------------------------
    # DPA/fpy (struct-origin corrected; NRT using Ed per element)
    # ----------------------------
    if not hasattr(bm, "materials") or not hasattr(bm.materials, "Ed"):
        raise RuntimeError("Expected bm.materials.Ed(element) to exist for DPA calculation.")

    dpa = np.zeros(len(cell_ids), float)
    dpa_std = np.zeros(len(cell_ids), float)

    for i, cid in enumerate(cell_ids):
        cid = int(cid)
        t = dpa_gas_map.get(cid)
        if t is None:
            continue

        denom_atoms = float(cell_total_atoms_struct.get(cid, 0.0))
        if denom_atoms <= 0.0:
            continue

        f_struct = cell_struct_origin_frac.get(cid, {}) or {}
        nuclides = list(t.nuclides or [])

        sum_mean = 0.0
        sum_var = 0.0

        for nuc in nuclides:
            nuc = str(nuc)
            fs = float(f_struct.get(nuc, 0.0))
            if fs == 0.0:
                continue

            try:
                dmg_mean, dmg_sd = tally_mean_std_for_nuclide(t, score="damage-energy", nuclide=nuc)
            except Exception:
                continue

            dmg_mean *= fs
            dmg_sd *= fs

            el = element_from_nuclide(nuc)
            Ed = float(bm.materials.Ed(el))  # eV
            disp_scale = 0.8 / (2.0 * Ed)    # NRT: displacements = 0.8*E_dmg/(2*Ed)

            disp_mean = disp_scale * source_per_y * dmg_mean
            disp_sd = disp_scale * source_per_y * dmg_sd

            dpa_mean = disp_mean / denom_atoms
            dpa_sd = disp_sd / denom_atoms

            sum_mean += dpa_mean
            sum_var += dpa_sd**2

        dpa[i] = sum_mean
        dpa_std[i] = math.sqrt(max(sum_var, 0.0))

    fig, ax = plt.subplots()
    ax.set_yscale("log")
    lo = np.maximum(dpa - dpa_std, 1e-30)
    hi = np.maximum(dpa + dpa_std, 1e-30)
    ax.step(xedges, np.r_[dpa, dpa[-1]], where="post")
    ax.fill_between(xedges, np.r_[lo, lo[-1]], np.r_[hi, hi[-1]], step="post", alpha=0.25)
    ax.set_xlabel("Radial position x [cm]")
    ax.set_ylabel("NRT dpa/fpy (struct-origin)")
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    fig.savefig(RESULTS_DIR / "dpa_fpy.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ----------------------------
    # CSV summary
    # ----------------------------
    rows = []
    for i, cid in enumerate(cell_ids):
        rows.append(
            dict(
                cell_id=int(cid),
                label=str(labels[i]),
                x_centroid_cm=float(xcent[i]),
                x_left_cm=float(xedges[i]),
                x_right_cm=float(xedges[i + 1]),
                total_flux_neutron=float(total_neut[i]),
                total_flux_photon=float(total_phot[i]),
                heating_W_per_cm3=float(heat_w[i]),
                H_appm_fpy=float(h_appm_y[i]),
                He_appm_fpy=float(he_appm_y[i]),
                dpa_fpy=float(dpa[i]),
            )
        )
    pd.DataFrame(rows).to_csv(RESULTS_DIR / "summary.csv", index=False)

    print(f"[ok] Wrote results to: {RESULTS_DIR}")


# =============================================================================
# Main
# =============================================================================
if __name__ == "__main__":
    with openmc.StatePoint(STATEPOINT_FILE) as sp:
        process(sp)

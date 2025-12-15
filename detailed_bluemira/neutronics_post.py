#!/usr/bin/env python3
"""
Post-processing script for the Bluemira / EU-DEMO OB sector.

This script:
  - reads the OpenMC statepoint file
  - loads geometry and materials from bluremira_model
  - post-processing and plotting results

Required files in the working directory:
  - statepoint.<batches>.h5
  - bluemira_model.py
  - EUDEMO_11_10d.h5m
  - materials/ (Python module with material definitions)
"""

import math
import os
import re
from collections import defaultdict
from enum import Enum

import matplotlib.pyplot as plt
import numpy as np
import openmc
import pandas as pd
import pydagmc

# easier to just reload full_setup
import neutronics_model as bm

import os
import sys
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'materials'))
sys.path.append(module_path)
import materials

# ----------------------------------------------------------------------
# User / model settings
# ----------------------------------------------------------------------

STATEPOINT_FILE = "statepoint.10.h5"
_DAGMC_MODEL_FILE = "EUDEMO_11_10d.h5m"

# Using this function to properly  match the plot of 
# the steps/error shades to the radial position of that layer
def make_bin_edges(centroids, widths):
    """
    Given layer centroids and widths, return bin edges for step plots.

    For N layers, returns N+1 edges:
      edge[0]   = centroid[0] - width[0]/2
      edge[1]   = centroid[0] + width[0]/2  (≈ start of 2nd layer)
      ...
      edge[N-1] = centroid[N-1] - width[N-1]/2
      edge[N]   = centroid[N-1] + width[N-1]/2
    """
    left = centroids - widths / 2.0
    right = centroids + widths / 2.0

    edges = np.empty(len(centroids) + 1, dtype=float)
    edges[:-1] = left
    edges[-1] = right[-1]
    return edges

# honestly I could have called some of these with bm.cell_ids etc
# OB equatorial cells (same as in your main script)
cell_ids = [49, 45, 79, 66, 82, 94, 103, 112, 124, 132, 2]

# Radial centroids [cm] for those layers, in same order as CELL_IDS_OB
xcentroids_ob = np.array((0.1, 1.1, 6.0, 15.0, 25.0,
                          35.0, 45.0, 55.0, 70.0, 90.0, 157.2))


step_lengths_ob = np.array([0.2, 1.8, 8.0, 10.0, 10.0,
                            10.0, 10.0, 10.0, 20.0, 20.0, 110.0])

# Source intensity -> scaling
total_power = 2e9 # 2000 MW
number_sectors = 16
section_power = total_power / number_sectors
ev_to_joule = 1.60218e-19
ev_fusion = 17.6e6 
convert_e = ev_to_joule * ev_fusion
neutron_source_rate = section_power / convert_e
s_in_y = (365 * 24 * 60 * 60) 

energies = openmc.mgxs.GROUP_STRUCTURES['CCFE-709']
unit_lethargy = [np.log(energies[i+1]/energies[i]) for i in range(len(energies)-1)]

all_cells = bm.all_cells
cell_ids = bm.cell_ids

# ---------------------------------------------------------------------
# Albedo calculation settings
# ---------------------------------------------------------------------

cell_id_to_name = {
    49:  "Armor",
    45:  "First Wall",
    79:  "OB_1",
    66:  "OB_2",
    82:  "OB_3",
    94:  "OB_4",
    103: "OB_5",
    112: "OB_6",
    124: "OB_7",
    132: "OB_8",
    2:   "VV",
}

# Surfaces to exclude from mean/uncertainty calculations
# (fill this with the IDs you want to drop)
EXCLUDED_SURFACES = [
    22,   # VV x+
    43,   # VV x-
    507,  # ARMOR x-
    896,  # OB_8 x+
    1075, #VV x+
    # The following surfaces are y,z faces
    # They are very small and originated from some CAD clean-up
    # Further model should clean those micro-surfaces
    # removed them so they dont bias the calculation
    1076, # (tiny surface)
    1077, # (tiny surface)
]

# ----------------------------------------------------------------------
# Main postprocessing
# ----------------------------------------------------------------------

with openmc.StatePoint(STATEPOINT_FILE) as sp:

    # bin edges for OB-style step plots
    x_edges_ob = make_bin_edges(xcentroids_ob, step_lengths_ob)
    x_ob = np.asarray(xcentroids_ob, dtype=float)

    # slice [0 : N_ob] -> kept this in case we want IB tallies
    idx_ob_start = 0
    idx_ob_end   = len(cell_ids)

    # -----------------------------------------
    # Scaling per cell
    # -----------------------------------------
    scaling_by_cell = {}
    for cid in cell_ids:
        vol = all_cells[cid].volume
        scaling_by_cell[cid] = (neutron_source_rate / vol)

    # -------------------------
    # Neutron & photon spectra
    # -------------------------
    particle_tally = sp.get_tally(id=bm.flux_tally.id)
    flux_mean = particle_tally.get_reshaped_data()
    flux_std  = particle_tally.get_reshaped_data(value='std_dev')

    # get only the neutrons
    neutron_flux      = flux_mean[:, 0, :]
    neutron_flux_std  = flux_std[:, 0, :]

    neutron_flux_ob     = neutron_flux[idx_ob_start:idx_ob_end, :]
    neutron_flux_ob_std = neutron_flux_std[idx_ob_start:idx_ob_end, :]

    # prepare csv flux data
    E_mid = 0.5 * (energies[:-1] + energies[1:])
    flux_data = {'Energy [eV]': E_mid}

    plt.figure()
    for i, cid in enumerate(cell_ids):
        scaling = scaling_by_cell[cid]
        flux_scaled = neutron_flux_ob[i].flatten() * scaling / unit_lethargy
        flux_data[f'Flux_cell_{cid} [n/cm^2/s/lethargy]'] = flux_scaled
        plt.loglog(
            energies[:-1],
            flux_scaled,
            label=f"Depth = {xcentroids_ob[i]:.2f} cm"
        )

    plt.legend()
    plt.grid()
    plt.ylabel('Neutron Flux Per Unit Lethargy [1/cm$^2$/s]')
    plt.xlabel('Energy [eV]')
    plt.xlim([1, 100e6])
    plt.savefig('n_flux_spectrum_ob.png', dpi=300)
    plt.close()

    # Neutron flux csv file for students at W armor (your current behavior: first flux column)
    df_flux = pd.DataFrame(flux_data)
    flux_col = df_flux.iloc[:, 1]
    df_first = pd.DataFrame({
        "Energy_mid[eV]": E_mid,
        "Flux_cell_armor[n/cm2/s/lethargy]": flux_col
    })
    df_first.to_csv("neutron_flux_spectrum.csv", index=False)

    # get only the photons
    photon_flux     = flux_mean[:, 1, :]
    photon_flux_std = flux_std[:, 1, :]

    photon_flux_ob     = photon_flux[idx_ob_start:idx_ob_end, :]
    photon_flux_ob_std = photon_flux_std[idx_ob_start:idx_ob_end, :]

    plt.figure()
    for i, cid in enumerate(cell_ids):
        scaling = scaling_by_cell[cid]
        plt.loglog(
            energies[:-1],
            photon_flux_ob[i].flatten() * scaling / unit_lethargy,
            label=f"Depth = {xcentroids_ob[i]:.2f} cm"
        )

    plt.legend()
    plt.grid()
    plt.ylabel('Photon Flux Per Unit Lethargy [1/cm$^2$/s/eV]')
    plt.xlabel('Energy [eV]')
    plt.xlim([1, 100e6])
    plt.savefig('p_flux_spectrum_ob.png', dpi=300)
    plt.close()

    # ------------------------------------------------------
    # Total flux directly from energy-integrated tally
    # ------------------------------------------------------
    # Fetch tally (by id or by name)
    t_flux_tot = sp.get_tally(id=bm.flux_tally_total.id)
    ncells_ob = len(cell_ids)
    n_cells   = len(cell_ids) 

    # Result shape depends on filters; easiest is to use get_reshaped_data()
    # Since filters are [CellFilter, ParticleFilter], we expect (cells, particles)
    flux_tot_mean = t_flux_tot.get_reshaped_data(value='mean').squeeze()
    flux_tot_std  = t_flux_tot.get_reshaped_data(value='std_dev').squeeze()

    # Ensure shape is (n_cells, n_particles)
    # If your cell filter includes exactly the same cells in the same order as cell_ids:
    #   n_particles should be 2 (neutron, photon)
    n_particles = flux_tot_mean.shape[1]

    # Slice OB cells (same slicing convention as your spectral tally)
    flux_tot_mean_ob = flux_tot_mean[idx_ob_start:idx_ob_end, :]
    flux_tot_std_ob  = flux_tot_std[idx_ob_start:idx_ob_end, :]

    # Convert per-source -> per-second using the same scaling_by_cell
    direct_total_neut = np.zeros(ncells_ob)
    direct_total_neut_std = np.zeros(ncells_ob)
    direct_total_phot = np.zeros(ncells_ob)
    direct_total_phot_std = np.zeros(ncells_ob)

    for i, cid in enumerate(cell_ids):
        sc = scaling_by_cell[cid]
        direct_total_neut[i]     = flux_tot_mean_ob[i, 0] * sc
        direct_total_neut_std[i] = flux_tot_std_ob[i, 0]  * sc
        direct_total_phot[i]     = flux_tot_mean_ob[i, 1] * sc
        direct_total_phot_std[i] = flux_tot_std_ob[i, 1]  * sc

    print('Maximum neutron flux (OB): ', f"{np.max(direct_total_neut):.4e}")
    print('Maximum photon flux  (OB): ', f"{np.max(direct_total_phot):.4e}")
    # ------------------------------------------------------
    # Compare: direct-total vs integrated-spectrum
    # ------------------------------------------------------
    fig, ax = plt.subplots()

    # Direct totals from flux_tally_total
    neut_lower_dir = np.maximum(direct_total_neut - direct_total_neut_std, 1e-2)
    neut_upper_dir = direct_total_neut + direct_total_neut_std
    phot_lower_dir = np.maximum(direct_total_phot - direct_total_phot_std, 1e-2)
    phot_upper_dir = direct_total_phot + direct_total_phot_std

    ax.step(x_edges_ob, np.r_[direct_total_neut, direct_total_neut[-1]],
            where='post', linestyle='--', label='Neutron')
    ax.fill_between(x_edges_ob, np.r_[neut_lower_dir, neut_lower_dir[-1]],
                    np.r_[neut_upper_dir, neut_upper_dir[-1]],
                    step='post', alpha=0.15, label='_nolegend_')

    ax.step(x_edges_ob, np.r_[direct_total_phot, direct_total_phot[-1]],
            where='post', linestyle='--', label='Photon')
    ax.fill_between(x_edges_ob, np.r_[phot_lower_dir, phot_lower_dir[-1]],
                    np.r_[phot_upper_dir, phot_upper_dir[-1]],
                    step='post', alpha=0.15, label='_nolegend_')

    ax.set_yscale('log')
    ax.set_ylabel('Total Flux [1/cm$^2$/s]')
    ax.set_xlabel('Radial Position [cm]')
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.legend()
    fig.savefig('flux_total_OB.png', bbox_inches='tight', dpi=300)
    plt.close(fig)

    # ------------------------------------------------------
    # Heating (OB)  
    # ------------------------------------------------------
    h_tally = sp.get_tally(id=bm.heating_tally.id)
    heating     = h_tally.get_values().flatten()
    heating_std = h_tally.get_values(value='std_dev').flatten()

    heating_ob     = heating[:ncells_ob]
    heating_std_ob = heating_std[:ncells_ob]

    def convert_heating(heating_group, heating_std_group, cell_ids_group):
        n = len(cell_ids_group)
        heat = np.zeros(n)
        heat_std = np.zeros(n)
        for i, cid in enumerate(cell_ids_group):
            scaling = scaling_by_cell[cid]
            heat[i]     = heating_group[i] * scaling * ev_to_joule
            heat_std[i] = heating_std_group[i] * scaling * ev_to_joule
        return heat, heat_std

    heat_ob, heat_std_ob = convert_heating(heating_ob, heating_std_ob, cell_ids)
    print('Maximum heating (OB): ', np.max(heat_ob))

    lower_h = np.maximum(heat_ob - heat_std_ob, 1e-16)
    upper_h = heat_ob + heat_std_ob

    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.step(
        x_edges_ob,
        np.r_[heat_ob, heat_ob[-1]],
        where='post',
        color='k',
        label='Heating (OB)'
    )
    ax.fill_between(
        x_edges_ob,
        np.r_[lower_h, lower_h[-1]],
        np.r_[upper_h, upper_h[-1]],
        step='post',
        alpha=0.3,
        edgecolor='none',
        facecolor='gray',
        label='Std. Dev.'
    )
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.set_ylabel('Heating [W/cm³]')
    ax.set_xlabel('Radial Position [cm]')
    ax.legend()
    fig.savefig("heating_uncertainty_OB.png", dpi=300, bbox_inches='tight')
    plt.close(fig)

    # ------------------------------------------------------
    # Hydrogen production (OB) 
    # ------------------------------------------------------
    hydrogen1_tally = sp.get_tally(id=bm.h1_tally.id)
    n_nuclides = len(bm.solid_nuclides)

    # mean/std arrays: shape -> (cells, nuclides)
    h1_mean = hydrogen1_tally.mean.squeeze().reshape(n_cells, n_nuclides)
    h1_std  = hydrogen1_tally.std_dev.squeeze().reshape(n_cells, n_nuclides)

    # Sum over parent nuclides -> per-cell production per source
    h1_per_cell_mean = h1_mean.sum(axis=1)
    h1_per_cell_std  = np.sqrt((h1_std**2).sum(axis=1))   # RSS across nuclides

    # OB slice
    h1_ob_mean = h1_per_cell_mean[:ncells_ob]
    h1_ob_std  = h1_per_cell_std[:ncells_ob]

    # Convert: per source -> per second -> per year -> appm/y
    h_ob = np.zeros(ncells_ob)
    h_ob_std = np.zeros(ncells_ob)

    for i, cid in enumerate(cell_ids):
        factor = neutron_source_rate * s_in_y / bm.cell_total_atoms_solid[cid] * 1e6
        h_ob[i]     = h1_ob_mean[i] * factor
        h_ob_std[i] = h1_ob_std[i]  * factor

    print('Maximum proton appm/y (OB): ', np.max(h_ob))

    lower_h1 = np.maximum(h_ob - h_ob_std, 1e-30)
    upper_h1 = h_ob + h_ob_std

    fig, ax = plt.subplots()
    ax.set_yscale('log')

    ax.step(
        x_edges_ob,
        np.r_[h_ob, h_ob[-1]],
        where='post',
        color='k',
        label='Proton [appm/y]'
    )
    ax.fill_between(
        x_edges_ob,
        np.r_[lower_h1, lower_h1[-1]],
        np.r_[upper_h1, upper_h1[-1]],
        step='post',
        alpha=0.3,
        edgecolor='none',
        facecolor='gray',
        label='Std. Dev.'
    )

    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.set_ylabel('Proton [appm/y]')
    ax.set_xlabel('Radial Position [cm]')
    ax.legend()
    fig.savefig('h1_OB.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    # ------------------------------------------------------
    # Helium production (OB)
    # ------------------------------------------------------
    helium3_tally = sp.get_tally(id=bm.he3_tally.id)
    helium4_tally = sp.get_tally(id=bm.he4_tally.id)

    he3_mean = helium3_tally.mean.squeeze().reshape(n_cells, n_nuclides)
    he3_std  = helium3_tally.std_dev.squeeze().reshape(n_cells, n_nuclides)

    he4_mean = helium4_tally.mean.squeeze().reshape(n_cells, n_nuclides)
    he4_std  = helium4_tally.std_dev.squeeze().reshape(n_cells, n_nuclides)

    # Sum over nuclides -> per cell (per source)
    he3_per_cell_mean = he3_mean.sum(axis=1)
    he3_per_cell_std  = np.sqrt((he3_std**2).sum(axis=1))

    he4_per_cell_mean = he4_mean.sum(axis=1)
    he4_per_cell_std  = np.sqrt((he4_std**2).sum(axis=1))

    # Add He-3 and He-4 (means add; std add in quadrature if independent)
    he_per_cell_mean = he3_per_cell_mean + he4_per_cell_mean
    he_per_cell_std  = np.sqrt(he3_per_cell_std**2 + he4_per_cell_std**2)

    # OB slice
    he_ob_src_mean = he_per_cell_mean[:ncells_ob]
    he_ob_src_std  = he_per_cell_std[:ncells_ob]

    # Convert to appm/y
    he_ob = np.zeros(ncells_ob)
    he_ob_std = np.zeros(ncells_ob)

    for i, cid in enumerate(cell_ids):
        factor = neutron_source_rate * s_in_y / bm.cell_total_atoms_solid[cid] * 1e6
        he_ob[i]     = he_ob_src_mean[i] * factor
        he_ob_std[i] = he_ob_src_std[i]  * factor

    print('Maximum helium appm/y (OB): ', np.max(he_ob))

    lower_he = np.maximum(he_ob - he_ob_std, 1e-30)
    upper_he = he_ob + he_ob_std

    fig, ax = plt.subplots()
    ax.set_yscale('log')

    ax.step(
        x_edges_ob,
        np.r_[he_ob, he_ob[-1]],
        where='post',
        color='k',
        label='Helium [appm/y]'
    )
    ax.fill_between(
        x_edges_ob,
        np.r_[lower_he, lower_he[-1]],
        np.r_[upper_he, upper_he[-1]],
        step='post',
        alpha=0.3,
        edgecolor='none',
        facecolor='gray',
        label='Std. Dev.'
    )

    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.set_ylabel('Helium [appm/y]')
    ax.set_xlabel('Radial Position [cm]')
    ax.legend()
    fig.savefig('he_OB.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    # ------------------------------------------------------
    # DPA (OB) 
    # ------------------------------------------------------
    dpa_per_y_ob_mean = np.zeros(ncells_ob)
    dpa_per_y_ob_std  = np.zeros(ncells_ob)

    # Loop over all dpa tallies
    for tally in bm.dpa_tallies:
        # get tally
        dpa_tally = sp.get_tally(id=tally.id)

        # find cell_id
        cell_filter_dpa = next(
            f for f in dpa_tally.filters if isinstance(f, openmc.CellFilter)
        )
        cid = cell_filter_dpa.bins[0]
        idx = cell_ids.index(cid)

        # Sum damage-energy value of all the nuclides of this dpa tally (1 element)
        dpa_mean = dpa_tally.summation(nuclides=dpa_tally.nuclides).mean.squeeze()
        dpa_std  = dpa_tally.summation(nuclides=dpa_tally.nuclides).std_dev.squeeze()

        # Convert damage-energy to dpa for that element
        Ed = materials.Ed(materials.element(dpa_tally.nuclides))
        scale_factor = 0.8 / (2 * Ed) * neutron_source_rate * s_in_y

        # Compute total displacements
        dpa_y_mean_cell = scale_factor * dpa_mean
        dpa_y_std_cell  = scale_factor * dpa_std

        # Scale displacements per total atoms (in the solid materials list)
        dpa_per_y_ob_mean[idx] += dpa_y_mean_cell / bm.cell_total_atoms_solid[cid]
        dpa_per_y_ob_std[idx]  += dpa_y_std_cell  / bm.cell_total_atoms_solid[cid]

    dpa_per_y_ob_mean = np.asarray(dpa_per_y_ob_mean)
    dpa_per_y_ob_std  = np.asarray(dpa_per_y_ob_std)

    print("Maximum DPA/y (OB): ", np.max(dpa_per_y_ob_mean))

    lower_dpa = np.maximum(dpa_per_y_ob_mean - dpa_per_y_ob_std, 1e-30)
    upper_dpa = dpa_per_y_ob_mean + dpa_per_y_ob_std

    fig, ax = plt.subplots()
    ax.set_yscale('log')

    ax.step(
        x_edges_ob,
        np.r_[dpa_per_y_ob_mean, dpa_per_y_ob_mean[-1]],
        where='post',
        color='k',
        label='DPA/y (OB)'
    )
    ax.fill_between(
        x_edges_ob,
        np.r_[lower_dpa, lower_dpa[-1]],
        np.r_[upper_dpa, upper_dpa[-1]],
        step='post',
        alpha=0.3,
        color='gray',
        label='Std. Dev.'
    )
    gap_start = 100.0

    gap_end   = 102.2

    ax.axvline(gap_start, color="k", linestyle="--", linewidth=1)
    ax.axvline(gap_end,   color="k", linestyle="--", linewidth=1)

    ax.text(
        0.5 * (gap_start + gap_end),
        ax.get_ylim()[1] * 0.7,
        "Gap",
        ha="right",
        va="top",
        rotation=90,
        fontsize=9
    )

    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.set_ylabel('DPA/y')
    ax.set_xlabel('Radial Position [cm]')
    ax.legend()
    fig.savefig('dpa_OB.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    # ------------------------------------------------------
    # Albedo post-processing 
    # ------------------------------------------------------
    info = bm.info      
    t_currents = bm.t_current_tally
    p_currents = bm.p_current_tallies
    
    # surface -> "owner" cell (first one wins) 
    # This should not be necessary if my previous surface list is correct
    # (should remove surfaces that are shared between cells in cell_ids)
    surf_to_cell = {}
    for cid in cell_ids:
        for s in info[cid]:
            sid = s["surface_id"]
            surf_to_cell.setdefault(sid, cid)

    # Total current tally (Jnet) 
    t_tot = sp.get_tally(id=t_currents.id)
    surf_filter = next(f for f in t_tot.filters if isinstance(f, openmc.SurfaceFilter))
    surf_ids_tot = list(surf_filter.bins)

    Jnet_mean = np.atleast_1d(t_tot.mean.squeeze())
    Jnet_std  = np.atleast_1d(t_tot.std_dev.squeeze())

    partial_mean = {}
    partial_std  = {}

    # Build tables for outgoing partial current Jout per (cell_id, surface_id).
    # Result: one mean/std value per surface for that cell.
    for cid, t in p_currents.items():
        tally = sp.get_tally(id=t.id)

        sfilter = next(
            f for f in tally.filters if isinstance(f, openmc.SurfaceFilter)
        )
        sids = list(sfilter.bins)

        mean = np.atleast_1d(tally.mean.squeeze())
        std  = np.atleast_1d(tally.std_dev.squeeze())

        # One surface per bin → 1D arrays
        for sid, m, s in zip(sids, mean, std):
            partial_mean[(cid, sid)] = float(m)
            partial_std[(cid, sid)]  = float(s)

    # Albedo per surface
    rows = []
    eps = 1e-15

    # Calculate the albedo mean and std for each surface in each cell.
    for idx, sid in enumerate(surf_ids_tot):
        if sid not in surf_to_cell:
            continue

        cid = surf_to_cell[sid]

        Jnet_m = float(Jnet_mean[idx])
        Jnet_s = float(Jnet_std[idx])

        Jout_m = partial_mean.get((cid, sid), 0.0)
        Jout_s = partial_std.get((cid, sid), 0.0)

        Jnet = abs(Jnet_m)
        Jout = abs(Jout_m)
        Jin  = abs(Jout - Jnet)

        if (Jout < eps):  # covers both special cases cleanly
            A_mean, A_std = 0.0, 0.0
        else:
            A_mean = Jin / Jout
            dA_dJout = Jnet / (Jout ** 2)
            dA_dJnet = -1.0 / Jout
            A_var = (dA_dJout ** 2) * (Jout_s ** 2) + (dA_dJnet ** 2) * (Jnet_s ** 2)
            A_std = math.sqrt(max(A_var, 0.0))

        rows.append(
            dict(
                surface_id=sid,
                cell_id=cid,
                J_total_mean=Jnet_m,
                J_total_std=Jnet_s,
                J_out_mean=Jout_m,
                J_out_std=Jout_s,
                albedo_mean=A_mean,
                albedo_std=A_std,
            )
        )

    df_alb = pd.DataFrame(rows).sort_values("surface_id")

    # provide csv file with all albedo values per surface
    #df_alb.to_csv("surface_albedo_all.csv", index=False)

    # Filter excluded surfaces
    # CSV of albedo surfaces (without filtered surfaces)
    df_alb_filt = df_alb[~df_alb["surface_id"].isin(EXCLUDED_SURFACES)].copy()
    df_alb_filt.to_csv("surface_albedo_filtered.csv", index=False)

    # Cell summaries
    # provide a csv file with the total albedo per cell (lateral surfaces)
    cell_rows = []
    for cid, grp in df_alb_filt.groupby("cell_id"):
        A_mean = grp["albedo_mean"].mean()
        A_std  = np.sqrt((grp["albedo_std"] ** 2).sum()) / max(len(grp), 1)
        cell_rows.append({
            "cell_id": cid,
            "cell_name": cell_id_to_name.get(cid, f"Cell {cid}"),
            "albedo_mean": A_mean,
            "albedo_std": A_std,
        })

    df_cell = pd.DataFrame(cell_rows).sort_values("cell_id")
    df_cell.to_csv("cell_albedo_summary.csv", index=False)

    # Special surfaces 
    special_surface_map = {
        507: "Armor",
        896: "OB_8",
        1075: "VV",
    }

    special_points = []
    for surf_id, cell_name in special_surface_map.items():
        row = df_alb[df_alb["surface_id"] == surf_id]
        if row.empty:
            print(f"Warning: surface {surf_id} not found in df_alb, skipping.")
            continue
        r = row.iloc[0]
        special_points.append(
            dict(
                surface_id=surf_id,
                cell_name=cell_name,
                albedo_mean=float(r["albedo_mean"]),
                albedo_std=float(r["albedo_std"]),
            )
        )

    # Plot: cell-average ± uncertainty + special surfaces (aligned x)
    ordered_layers = [cell_id_to_name[cid] for cid in cell_ids]
    summary = df_cell.set_index("cell_name")
    summary = summary.loc[[name for name in ordered_layers if name in summary.index]].reset_index()

    # numeric x positions for each layer label 
    xlabels = list(summary["cell_name"])
    xpos = np.arange(len(xlabels))
    xmap = {name: i for i, name in enumerate(xlabels)}

    plt.figure(figsize=(11, 6))

    # -------------------------------
    # 1) Cell average 
    # -------------------------------
    plt.errorbar(
        xpos,
        summary["albedo_mean"].values,
        yerr=summary["albedo_std"].values,
        fmt="none",
        elinewidth=1,
        capsize=4,
    )
    plt.scatter(
        xpos,
        summary["albedo_mean"].values,
        marker="x",
        s=90,
        linewidths=2,
        label="Cell average (lateral surfaces)",
    )

    # -------------------------------
    # 2) Special surfaces 
    # -------------------------------
    special_surface_labels = {
        507:  "Frontal surface (Armor, x−)",
        896:  "Backward surface (OB_8, x+)",
        1075: "Backward surface (VV, x+)",
    }

    special_colors = {
    507:  "orange",
    896:  "green",
    1075: "red",
    }

    for sp_pt in special_points:
        sid  = sp_pt["surface_id"]
        name = sp_pt["cell_name"]

        if name not in xmap:
            print(f"Warning: special surface {sid} has cell_name={name} not in plot labels, skipping.")
            continue

        x = xmap[name]
        y = sp_pt["albedo_mean"]
        e = sp_pt["albedo_std"]

        label = special_surface_labels.get(sid, f"Surface {sid}")
        color = special_colors.get(sid, "black")

        # error bars 
        plt.errorbar(
            [x],
            [y],
            yerr=[[e], [e]],
            fmt="none",
            ecolor=color,       
            elinewidth=1,
            capsize=4,
        )

        # X marker 
        plt.scatter(
            [x],
            [y],
            marker="x",
            color=color,
            s=140,
            linewidths=2.5,
            label=label,
        )

    # ----------------------
    # Axes formatting
    # ----------------------
    plt.xticks(xpos, xlabels, rotation=45, fontsize=11)
    plt.xlabel("Layer", fontsize=12)
    plt.ylabel("Albedo", fontsize=12)
    plt.title("Albedo per Layer", fontsize=15)
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    plt.tight_layout()

    # de-duplicate legend entries
    handles, labels = plt.gca().get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    plt.legend(unique.values(), unique.keys(), fontsize=11)

    plt.savefig("albedo_layers_filtered.png", dpi=300)
    plt.close()
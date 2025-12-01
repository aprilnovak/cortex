# built-in modules
from collections import defaultdict
from enum import Enum
import math
import re

import openmc
import numpy as np
import matplotlib.pyplot as plt
#import materials
from openmc_plasma_source import tokamak_source
import pandas as pd
import pydagmc
import os
import sys

module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'materials'))
sys.path.append(module_path)
import materials

_DAGMC_MODEL_FILE = 'EUDEMO_11_10d.h5m'

model = openmc.Model()

# --------------------------------
# GEOMETRY
# -------------------------------- 

# Import DAGMC geometry
dagmc_universe = openmc.DAGMCUniverse(filename='EUDEMO_11_10d.h5m')

pydagmc_model = pydagmc.Model(str(dagmc_universe.filename))

# reserve volume and surface IDs in the DAGMC model to avoid overlaps
openmc.reserve_ids([v.id for v in pydagmc_model.volumes], cls=openmc.Cell)
openmc.reserve_ids([s.id for s in pydagmc_model.surfaces], cls=openmc.Surface)

# Create planes for reflective BC
def azimuthal_plane(theta_deg, boundary_type=None, name=None, surface_id=None):
    theta = math.radians(theta_deg)
    if abs(theta_deg - 90.0) < 1e-6:
        a = 1.0   
        b = 0.0
        c = 0.0
        d = 0.0
        return openmc.Plane(a=a, b=b, c=c, d=d,
                            boundary_type=boundary_type,
                            name=name, surface_id=surface_id)
    
    a = -math.tan(theta)
    b = 1.0
    c = 0.0
    d = 0.0

    return openmc.Plane(a=a, b=b, c=c, d=d,
                        boundary_type=boundary_type,
                        name=name, surface_id=surface_id)

# Angles for planes
number_sectors = 16
theta0 = 0.0            # starting angle in degrees
theta1 = theta0 + 360/number_sectors  # 1/16 sector

# Reflective azimuthal cuts
cut_lo = azimuthal_plane(theta0, boundary_type='reflective', surface_id=100_001)
cut_hi = azimuthal_plane(theta1, boundary_type='reflective', surface_id=100_002)

# Outer BC graveyard (top/bottom)
z_min = openmc.ZPlane(z0=-2500, boundary_type='vacuum', surface_id=100_003) 
z_max = openmc.ZPlane(z0=+1500, boundary_type='vacuum', surface_id=100_004)

# Outer BC graveyard (radial)
r_out = openmc.ZCylinder(r=2500, boundary_type='vacuum', surface_id=100_005) 
sector_region = (+cut_lo & -cut_hi) & (+z_min & -z_max) & (-r_out)

# Put the DAGMC universe inside the CSG cell
sector_cell = openmc.Cell(cell_id=10_000, region=sector_region, fill=dagmc_universe, name='sector_container')
model.geometry = openmc.Geometry(root=[sector_cell])

# --------------------------------
#  Materials 
# --------------------------------
# This list will change depending on the BREEDER MODEL
ss316 = materials.ss316(7.93) # SS316-LN DOI 10.1088/1741-4326/ac2a6b
ccz = materials.CuCrZr(8.9) #CuCrZr DOI 10.1088/1741-4326/ac2a6b
eurofer = materials.eurofer97(7.87) # EUROFER97 DOI 10.1088/1741-4326/ac2a6b
t = materials.W(19.3) #Tungsten DOI 10.1088/1741-4326/ac2a6b
water = materials.Water(0.866) #H20 https://doi.org/10.1016/j.fusengdes.2020.111833
h = materials.Helium(0.0001785) # Helium PNNL COMPENDIUM
# Still missing definition
nb3sn = materials.Nb3Sn(5.7) #Nb3Sn
epoxy = materials.Epoxy(1.207) #Epoxy resin
bronze = materials.Bronze(8.8775) #Bronze
be = materials.Be(1.85) # Beryllium PNNL COMPENDIUM
nbti = materials.NbTi(6.538) #NbTi
c = materials.Cu(8.96) #Cu PNNL COMPENDIUM
# check function
PbLi = materials.PbLi(0.60, 9.8)

# ATRIBUTE MATERIALS TO DAGMC geometry
# Several references use Void as a substitute for He cooling elements
# Armor
Armor = materials.W(19.3) #Tungsten DOI 10.1088/1741-4326/ac2a6b
Armor.name = 'Armor'
# Cryostat (From Bluemira Model)
CS = materials.ss316(7.93) # SS316-LN DOI 10.1088/1741-4326/ac2a6b
CS.name = 'CS'
# Radiation shielding # where did they use concrete? (find reference)
RS = materials.concrete(2.3) # PNNL COMPENDIUM
RS.name = 'RS'
# Divertor (From Bluemira Model)
divertor = openmc.Material.mix_materials([ccz, c, eurofer, t, water], [0.00552, 0.00438, 0.5238,  0.01026, 0.45604], 'vo',name="Divertor")

# --------------------------------------------------------------------
# homogenization from https://doi.org/10.1016/j.fusengdes.2020.111833
# Breeder Design WCLL
# First Wall
pf = openmc.Material.mix_materials([t, water, eurofer],[0.0027, 0.14268, 0.85468],'vo',name='FW')
# Inner Breeder Blanket
matIB1 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.833, 0.025, 0.139], 'vo',name="IB1")
matIB2 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.858, 0.018, 0.124], 'vo',name="IB2")
matIB3 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.8132, 0.0158, 0.171], 'vo',name="IB3")
matIB4 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.427, 0.016, 0.558], 'vo',name="IB4")
matIB5 = openmc.Material.mix_materials([water, eurofer], [0.486, 0.514], 'vo',name="IB5")

# Outer Breeder Blanket
matOB1 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.833, 0.025, 0.139], 'vo',name="OB1")
matOB2 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.858, 0.018, 0.124], 'vo',name="OB2")
matOB3 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.8132, 0.0158, 0.171], 'vo',name="OB3")
matOB4 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.8132, 0.0158, 0.171], 'vo',name="OB4")
matOB5 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.8132, 0.0158, 0.171], 'vo',name="OB5")
matOB6 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.8132, 0.0158, 0.171], 'vo',name="OB6")
matOB7 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.427, 0.016, 0.558], 'vo',name="OB7")
matOB8 = openmc.Material.mix_materials([water, eurofer], [0.486, 0.514], 'vo',name="OB8")
# --------------------------------------------------------------------

# Vacuum Vessel
VV = openmc.Material.mix_materials([ss316, water],[0.6, 0.4],'vo',name='VV')

# Port fillings 
portf = openmc.Material.mix_materials([ss316, water], [0.6, 0.4], 'vo', name='PF')

# Poloidal coils (From Bluemira Model)
PC = openmc.Material.mix_materials([nbti, c, epoxy, bronze, h, ss316], [0.02895, 0.1169, 0.18, 0.0735, 0.1682, 0.4319], 'vo',name="PC")

# Toroidal Coils (From Bluemira Model)
TFcoil = openmc.Material.mix_materials([nb3sn, c, epoxy, bronze, h, ss316], [0.02895, 0.1169, 0.18,  0.0735, 0.1682, 0.4319], 'vo',name="TFC")

# Build materials model
model.materials = openmc.Materials([Armor, pf, portf, VV, divertor, TFcoil,
                                    matIB1, matIB2, matIB3, matIB4, matIB5,
                                    matOB1, matOB2, matOB3, matOB4, matOB5,
                                    matOB6, matOB7, matOB8 , PC, CS, RS])

# --------------------------------
#  SOURCE
# --------------------------------
# Openmc_plasma_source
my_source = tokamak_source(
    angles=(0.0, math.pi/8),
    elongation=1.739,
    ion_density_centre=6.8e19,
    ion_density_pedestal=5.78e19,
    ion_density_peaking_factor=1,
    ion_density_separatrix=1.02e19,
    ion_temperature_centre=23.7e3,
    ion_temperature_pedestal=5.5e3,
    ion_temperature_separatrix=0.1e3,
    ion_temperature_peaking_factor=8.06,
    ion_temperature_beta=6,
    major_radius=840.67,
    minor_radius=300.2,
    pedestal_radius=0.94 * 300.2,
    mode="H",
    shafranov_factor=0.44789,
    triangularity=0.333,
    fuel={"D": 0.5, "T": 0.5},
)

# --------------------------------
# Settings
# --------------------------------

# Simulation Settings
model.settings = openmc.Settings()
model.settings.dagmc = True
model.settings.photon_transport = True
model.settings.batches = 10
model.settings.particles = 10_000
model.settings.run_mode = "fixed source"
model.settings.source = my_source

# --------------------------------
# VOLUME CALCULATION
# --------------------------------
# reset the OpenMC ID space for cells and surfaces
openmc.Cell.reset_ids()
openmc.Surface.reset_ids()

# apply volumes from PyDAGMC to the OpenMC model cells
model.init_lib(output=False)
model.sync_dagmc_universes()
model.finalize_lib()

# IS THIS NECESSARY?
# reclaim the ID space for all cells and surfaces currently in the model
# to avoid clashes later in the model
openmc.reserve_ids([c_id for c_id in model.geometry.get_all_cells()], cls=openmc.Cell)
openmc.reserve_ids([s_id for s_id in model.geometry.get_all_surfaces()], cls=openmc.Surface)

dagmc_universe_cells = dagmc_universe.get_all_cells()
for volume in pydagmc_model.volumes:
    dagmc_universe_cells[volume.id].volume = volume.volume

# list of equatorial cells in Bluemira_11_10d.h5m model. (will change with other model)
# OB - outer blankets  (IB version can provided, currently removed the IB1 for depletion calculations)
cell_list_ob_1 = [49, 45, 79, 66, 82, 94, 103, 112, 124, 132, 2]
cell_ids = cell_list_ob_1

# centroids of layers
# if outer blanket (Armor, FW, OB1, OB2, OB3, OB4, OB5, OB6, OB7, OB8, VV)
xcentroids_ob = (0.1, 1.1, 6.0, 15.0, 25.0, 35.0, 45.0, 55.0, 70.0, 90.0, 157.2)

# INSERT ROUTINE TO DETERMINE THE CENTROIDS AUTOMATICALLY
#
#


# -------------
# Updates on Patrick's version of dagmc_surface functions
# SURFACES
class Orientation(Enum):
    FORWARD = 1 # indicates normal vector points outwards from the volume
    REVERSE = -1 # indicates normal vector points inwards to the volume

def dagmc_volume_surface_info(pydagmc_model, volume_ids):
    """
    For each volume in volume_ids, returns info for all its surfaces:
    - surface_id
    - orientation w.r.t. that volume
    - list of adjacent volume IDs for that surface, EXCLUDING any
      adjacent volumes that are in volume_ids.

    Returns
    -------
    dict : vol_id with surface_ids, orientation, and adjacent volumes
 
    """
    volume_ids_set = set(volume_ids)
    result = {}
    surface_ids = set()

    for vol_id in volume_ids:
        volume = pydagmc_model.volumes_by_id[vol_id]
        surf_list = []
        
        for surface in volume.surfaces:
            # --- Orientation relative to this volume ---
            parent_volumes = surface.senses
            if parent_volumes[0].id == vol_id:
                orientation = Orientation.FORWARD
            else:
                orientation = Orientation.REVERSE

            # --- Adjacent volumes via this surface (other than vol_id) ---
            all_adj_ids = [v.id for v in surface.volumes if v.id != vol_id]

            # Adjacent volumes we care about: NOT in the input volume_ids
            ext_adj_ids = [aid for aid in all_adj_ids if aid not in volume_ids_set]

            # Store surfaces that are facing external surfaces to the radial tallies block. 
            # - (Keep BC or gaps) -> If there are no adjacents at all
            if not all_adj_ids:
                # Pure boundary surface (vacuum, graveyard, etc.)
                surf_list.append({
                    "surface_id": surface.id,
                    "adjacent_volumes": [],
                    "orientation": orientation
                })
                surface_ids.add(surface.id)
            # - (Keep surfaces with other volumes that are not in the cell_ids list) -> If there are adjacents and at least one is outside volume_ids
            elif ext_adj_ids:
                # Has at least one neighbor outside the region of interest
                surf_list.append({
                    "surface_id": surface.id,
                    "adjacent_volumes": ext_adj_ids,
                    "orientation": orientation
                })
                surface_ids.add(surface.id)
            # - (Skip surfaces that have both adjacent volumes are inside the cell_ids list) -> If all adjacents are inside volume_ids
            else:
                # All neighbors are inside volume_ids -> skip
                continue

        result[vol_id] = surf_list

    return result, sorted(surface_ids)

# Dict and surface list
info, all_surface_ids = dagmc_volume_surface_info(pydagmc_model, cell_ids)

#Print
#for cid, surfaces in info.items():
#    print(f"\nCell {cid}:")
#    for s in surfaces:
#        adj_str = ", ".join(str(a) for a in s["adjacent_volumes"]) or "None"
#        print(f"  Surface {s['surface_id']}: "
#              f"{adj_str},"
#              f"{s['orientation'].name} ")

# --------------------------------
# TALLY
# --------------------------------
# Source intensity -> scaling
total_power = 2e9 # 2000 MW
section_power = total_power / number_sectors
ev_to_joule = 1.60218e-19
ev_fusion = 17.6e6 
convert_e = ev_to_joule * ev_fusion
neutron_source_rate = section_power / convert_e
s_in_y = (365 * 24 * 60 * 60) 

# FILTERS
cell_filter = openmc.CellFilter(cell_ids)
particle_filter = openmc.ParticleFilter(bins=['neutron', 'photon']) # for neutron and photon flux
t_surf_filter = openmc.SurfaceFilter(all_surface_ids) # current tallies
n_particle_filter = openmc.ParticleFilter(bins=['neutron']) # current tallies
energies = openmc.mgxs.GROUP_STRUCTURES['CCFE-709']
energy_filter = openmc.EnergyFilter(energies)
unit_lethargy = [np.log(energies[i+1]/energies[i]) for i in range(len(energies)-1)]

model.tallies = openmc.Tallies()

# Tallies
# flux score, with an energy filter
# Neutron and gamma flux spectrum
flux_tally = openmc.Tally()
flux_tally.filters = [cell_filter, particle_filter, energy_filter]
flux_tally.scores = ['flux']
model.tallies.append(flux_tally)

# Current
# Total current
t_current_tally = openmc.Tally()
t_current_tally.filters = [t_surf_filter, n_particle_filter] # just do for neutrons
t_current_tally.scores  = ['current']
model.tallies.append(t_current_tally)

# Partial current
# The idea is to create a tally for each surface from a cell_id
# Only creates the tally if surface and cellfromfilter match
p_current_tallies = {}
for cid in cell_ids:
    # get cells from cell_ids
    ocell = dagmc_universe_cells[cid]
    # find surfaces for that cell
    surf_ids_for_cell = [s["surface_id"] for s in info[cid]]

    if not surf_ids_for_cell:
        continue

    # create filter
    cell_from_filter = openmc.CellFromFilter([ocell])
    surf_filter = openmc.SurfaceFilter(surf_ids_for_cell)

    current_tally = openmc.Tally()
    current_tally.filters = [cell_from_filter, surf_filter, n_particle_filter]
    current_tally.scores  = ["current"]

    model.tallies.append(current_tally)
    p_current_tallies[cid] = current_tally

# heating
heating_tally = openmc.Tally()
heating_tally.filters = [cell_filter]
heating_tally.scores = ['heating']
model.tallies.append(heating_tally)

# H1, He3, and He4 production
# Create a list of nuclide for the EUROFER, tungstein, SS316 for He3,He4 and H1 production
# Avoid using most of the elements in the breeder material and water for gas generation
# However, there are still some Oxygen and other elements that are common with other materials 
# (Current estimations are better than before)
material_list = [eurofer, t, ss316]
solid_nuclides = set()
for mat in material_list:
    nuclides = mat.get_nuclide_densities()
    for name in nuclides.keys():
        solid_nuclides.add(name)

# Sorted list of nuclide names
solid_nuclides = sorted(solid_nuclides) 

# proton(hydrogen) production tallies
h1_tally = openmc.Tally()
h1_tally.filters = [cell_filter] 
h1_tally.scores = ['H1-production']
h1_tally.nuclides = solid_nuclides
model.tallies.append(h1_tally)

# helium production tallies
he3_tally = openmc.Tally()
he3_tally.filters = [cell_filter] 
he3_tally.scores = ['He3-production']
he3_tally.nuclides = solid_nuclides
model.tallies.append(he3_tally)

he4_tally = openmc.Tally()
he4_tally.filters = [cell_filter] 
he4_tally.scores = ['He4-production']
he4_tally.nuclides = solid_nuclides
model.tallies.append(he4_tally)

# -------------------------------
# DPA tallies
# -------------------------------

# Collect all nuclides over selected cells
all_cells = model.geometry.get_all_cells()
all_nuclides = set()
# Atoms calculations
cell_nuclide_atoms = {}   # {cell_id: {nuclide: mean_atoms}}
cell_total_atoms   = {}   # {cell_id: mean_total_atoms}
cell_total_atoms_solid = {} # {cell_id: mean_total_atoms_solid}

# Obtain cell_total_atoms and cell_total_atoms_solid
for cid in cell_ids:
    cell = all_cells[cid]
    all_nuclides.update(cell.fill.get_nuclides())
    
    cell.nuclide_atoms = cell.fill.get_nuclide_atoms(volume=cell.volume)
    cell.total_atoms = sum(cell.nuclide_atoms.values())

    cell_nuclide_atoms[cid] = cell.nuclide_atoms
    cell_total_atoms[cid]   = cell.total_atoms

    # total sum of atoms in solid nuclides in cell [cid] 
    solid_total_atoms = sum(
        atoms for nuc, atoms in cell.nuclide_atoms.items()
        if nuc in solid_nuclides
    )
    cell_total_atoms_solid[cid] = solid_total_atoms
    ratio_solid = cell_total_atoms_solid[cid] / cell_total_atoms[cid]
    #print('solid ratio:', ratio_solid)


# CREATE DPA TALLIES
# Current version is generating tallies for each element existent in their respective cell.
# Also filters to the elements present in the solid_nuclides list.
# (Previous version had all elements for all cells)
cell_filters_by_cid = {cid: openmc.CellFilter([cid]) for cid in cell_ids}

dpa_tallies = []
for cid in cell_ids:
    nuclide_atoms = cell_nuclide_atoms[cid]
    el_to_isos = defaultdict(set)
    for nuc in nuclide_atoms:
        el = re.sub(r'\d+', '', nuc)
        el_to_isos[el].add(nuc)

    for el, iso_set in el_to_isos.items():
        dpa_tally = openmc.Tally()
        dpa_tally.scores   = ["damage-energy"]
        dpa_tally.filters  = [cell_filters_by_cid[cid]]
        dpa_tally.nuclides = sorted(iso_set)
        model.tallies.append(dpa_tally)
        dpa_tallies.append(dpa_tally)

# Export tallies.xml
model.tallies.export_to_xml()    

# --------------------------------
# EXECUTE (and plot)
# --------------------------------
statepoint = model.run()
with openmc.StatePoint(statepoint) as sp:

    # OB cells 
    # Current model doesn't have IB (it was simplified to run with depletion)
    idx_ob_start = 0
    idx_ob_end   = len(cell_ids)        # slice [0 : N_ob]

    # Define the scaling for each cell.volume in the tally
    scaling_by_cell = {}
    for cid in cell_ids:
        cell_cid = all_cells[cid]
        vol = cell_cid.volume
        scaling = 1.0 / vol * neutron_source_rate
        scaling_by_cell[cid] = scaling    

    # Neutron and Gamma Spectrum
    particle_tally = sp.get_tally(id=flux_tally.id)
    flux_mean = particle_tally.get_reshaped_data()                 
    flux_std  = particle_tally.get_reshaped_data(value='std_dev') 

    # Neutrons lists
    neutron_flux        = flux_mean[:, 0, :]   
    neutron_flux_std   = flux_std[:, 0, :]    

    # OB
    neutron_flux_ob = neutron_flux[idx_ob_start:idx_ob_end, :]
    neutron_flux_ob_std = neutron_flux_std[idx_ob_start:idx_ob_end, :]

    # prepare csv flux data
    E_mid = 0.5 * (energies[:-1] + energies[1:])
    flux_data = {'Energy [eV]': E_mid}

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
    plt.savefig('n_flux_spectrum_ob.png')
    plt.close()

    # Neutron flux csv file for students
    df_flux = pd.DataFrame(flux_data)
    # Just the flux on the first cell (Armor)
    df_first = df_flux.iloc[:, [0]] 
    df_first.to_csv('neutron_flux_spectrum.csv', index=False)

    # Photon Spectrum
    photon_flux     = flux_mean[:, 1, :]   
    photon_flux_std = flux_std[:, 1, :]

    # OB
    photon_flux_ob  = photon_flux[idx_ob_start:idx_ob_end, :]
    photon_flux_ob_std  = photon_flux_std[idx_ob_start:idx_ob_end, :]
    
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
    plt.savefig('p_flux_spectrum_ob.png')
    plt.close()

    # ------------------------------------------------------
    # Total fluxes by integrating over energy (OB)
    # ------------------------------------------------------
    
    # some left over from previous model with IB (will work on monday for the clean-up)
    ncells_ob = len(cell_ids)
    n_cells = len(cell_ids)

    # I would like to talk briefly about the total flux in our meeting
    # I did some testing with the total flux tally (no energy tally) and obtained the same results
    def integrate_over_energy(flux_E, std_E, cell_ids_group):
        """Integrate over energy for each cell group (OB)."""
        n_cells = flux_E.shape[0]
        total_flux = np.zeros(n_cells)
        total_std  = np.zeros(n_cells)

        for i, cid in enumerate(cell_ids_group):
            scaling = scaling_by_cell[cid]

            f = flux_E[i].flatten()
            s = std_E[i].flatten()

            total_flux[i] = np.sum(f) * scaling
            total_std[i]  = np.sqrt(np.sum(s**2)) * scaling

        return total_flux, total_std

    # ---- Neutron / Photon integrated flux: OB ----
    total_neutron_flux_ob, total_neutron_std_ob = integrate_over_energy(
        neutron_flux_ob, neutron_flux_ob_std, cell_ids
    )
    total_photon_flux_ob, total_photon_std_ob = integrate_over_energy(
        photon_flux_ob, photon_flux_ob_std, cell_ids
    )

    print('Maximum neutron flux (OB): ', np.max(total_neutron_flux_ob))
    print('Maximum photon flux  (OB): ', np.max(total_photon_flux_ob))

    x_ob = np.asarray(xcentroids_ob)

    fig, ax = plt.subplots()

    # Neutrons OB
    ax.step(x_ob, total_neutron_flux_ob, where='mid', label='Neutron flux (OB)')
    neut_lower_ob = np.maximum(total_neutron_flux_ob - total_neutron_std_ob, 1e-30)
    neut_upper_ob = total_neutron_flux_ob + total_neutron_std_ob
    ax.fill_between(
        x_ob, neut_lower_ob, neut_upper_ob,
        step='mid', alpha=0.3, label='_nolegend_'
    )

    # Photons OB
    ax.step(x_ob, total_photon_flux_ob, where='mid', label='Photon flux (OB)')
    phot_lower_ob = np.maximum(total_photon_flux_ob - total_photon_std_ob, 1e-30)
    phot_upper_ob = total_photon_flux_ob + total_photon_std_ob
    ax.fill_between(
        x_ob, phot_lower_ob, phot_upper_ob,
        step='mid', alpha=0.3, label='_nolegend_'
    )

    ax.set_yscale('log')
    ax.set_ylabel('Total Flux [1/cm$^2$/s]')
    ax.set_xlabel('Radial Position [cm]')
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.legend()
    fig.savefig('flux_int_OB.png', bbox_inches='tight')
    plt.close(fig)

    # ------------------------------------------------------
    # Current / Albedo Calculation
    # ------------------------------------------------------
    # Total Current per surface
    t_tot = sp.get_tally(id=t_current_tally.id)

    # Surface IDs for total current
    t_surf_filter = next(f for f in t_tot.filters if isinstance(f, openmc.SurfaceFilter))
    surf_ids = list(t_surf_filter.bins)

    J_total = t_tot.mean.squeeze()  

    # Obtain Partial Currents
    partial_currents = {}

    for cid in cell_ids:
        t = sp.get_tally(id=p_current_tallies[cid].id)


        sfilter = next(f for f in t.filters if isinstance(f, openmc.SurfaceFilter))
        sids    = list(sfilter.bins)
        mean    = t.mean.squeeze()

        for sid, J in zip(sids, mean):
            partial_currents[(cid, sid)] = float(J)

    surf_to_cell = {}
    for cid in cell_ids:
        for s in info[cid]:
            sid = s["surface_id"]
            surf_to_cell.setdefault(sid, cid)
    
    # ALBEDO Calculation
    albedo_by_surface = {}   # sid -> dict of info
    eps = 1e-12  # Small number to avoid division by zero

    for sid, Jnet in zip(surf_ids, J_total):
        cid = surf_to_cell[sid]
        J_out = partial_currents.get((cid, sid), 0.0)

        # Partial current inward the volume 
        J_in = J_out - Jnet

        # Absolute value 
        J_in_mag = max(abs(J_in), eps) # (avoid zeros)
        J_out_mag = abs(J_out)

        # what happens if J_in > J_out?
        # Would like feedback on this one
        A = J_out_mag / J_in_mag

        albedo_by_surface[sid] = {
            "cell_id": cid,
            "J_total": Jnet,
            "J_out": J_out,
            "J_in": J_in,
            "albedo": A,
        }

    rows = []

    for sid, data in albedo_by_surface.items():
        rows.append({
            "surface_id": sid,
            "cell_id": data["cell_id"],
            "albedo": data["albedo"],
        })

    df = pd.DataFrame(rows).sort_values("surface_id")
    df.to_csv("surface_albedo.csv", index=False)

    # ------------------------------------------------------
    # Heating (OB)
    # ------------------------------------------------------
    h_tally = sp.get_tally(id=heating_tally.id)
    heating     = h_tally.get_values().flatten()                 
    heating_std = h_tally.get_values(value='std_dev').flatten()

    # OB (in case introduced IB later)
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

    # OB
    heat_ob, heat_std_ob = convert_heating(heating_ob, heating_std_ob, cell_ids)
    print('Maximum heating (OB): ', np.max(heat_ob))

    plt.figure()
    plt.yscale('log')
    plt.step(x_ob, heat_ob, where='mid', color='k', label='Heating (OB)')
    plt.fill_between(
        x_ob,
        np.maximum(heat_ob - heat_std_ob, 1e-16),
        heat_ob + heat_std_ob,
        step='mid',
        alpha=0.3,
        edgecolor='none',
        facecolor='gray',
        label='Std. Dev.'
    )
    plt.grid()
    plt.ylabel('Heating [W/cm³]')
    plt.xlabel('Radial Position [cm]')
    plt.legend()
    plt.savefig("heating_uncertainty_OB.png", dpi=300)
    plt.close()

    # ------------------------------------------------------
    # Hydrogen production (OB)
    # ------------------------------------------------------

    hydrogen1_tally = sp.get_tally(id=h1_tally.id)
    n_nuclides = len(solid_nuclides)
    
    # reshape to (cells, nuclides)
    h1_mean = hydrogen1_tally.mean.squeeze().reshape(n_cells, n_nuclides)

    # sum over parent nuclides -> per-cell production per source neutron
    h1_per_cell = h1_mean.sum(axis=1)

    # OB 
    h1_ob = h1_per_cell[:ncells_ob]

    h_ob = []
    for i, cid in enumerate(cell_ids):
        h_per_s = h1_ob[i] * neutron_source_rate
        h_per_y = h_per_s * s_in_y
        h_appm_per_y = h_per_y / cell_total_atoms_solid[cid] * 1e6
        h_ob.append(h_appm_per_y)

    h_ob = np.array(h_ob, dtype=float)
    print('Maximum proton appm/y (OB): ', np.max(h_ob))

    plt.figure()
    plt.yscale('log')
    plt.step(x_ob, h_ob, marker='o', color='k', markersize=1.5)
    plt.grid()
    plt.ylabel('Proton [appm/y]')
    plt.xlabel('Radial Position [cm]')
    plt.savefig('h1_OB.png')
    plt.close()

    # ------------------------------------------------------
    # Helium production (OB)
    # ------------------------------------------------------
    helium3_tally = sp.get_tally(id=he3_tally.id)
    helium4_tally = sp.get_tally(id=he4_tally.id)

    # reshape to (cells, nuclides)
    he3_mean = helium3_tally.mean.squeeze().reshape(n_cells, n_nuclides)
    he4_mean = helium4_tally.mean.squeeze().reshape(n_cells, n_nuclides)

    # sum over parent nuclides -> per-cell He3/He4 production per source neutron
    he3_per_cell = he3_mean.sum(axis=1)   
    he4_per_cell = he4_mean.sum(axis=1)   

    # total helium production (He3 + He4) per source neutron, per cell
    he_per_cell = he3_per_cell + he4_per_cell

    # OB 
    he_ob_src = he_per_cell[:ncells_ob]   # per source neutron, OB cells

    he_ob = []
    for i, cid in enumerate(cell_ids):
        he_per_s = he_ob_src[i] * neutron_source_rate
        he_per_y = he_per_s * s_in_y
        he_appm_per_y = he_per_y / cell_total_atoms_solid[cid] * 1e6
        he_ob.append(he_appm_per_y)

    he_ob = np.array(he_ob, dtype=float)
    print('Maximum helium appm/y (OB): ', np.max(he_ob))

    plt.figure()
    plt.yscale('log')
    plt.step(x_ob, he_ob, marker='o', color='k', markersize=1.5)
    plt.grid()
    plt.ylabel('Helium [appm/y]')
    plt.xlabel('Radial Position [cm]')
    plt.savefig('he_OB.png')
    plt.close()

    # ------------------------------------------------------
    # DPA (OB) 
    # ------------------------------------------------------
    dpa_per_y_ob_mean = np.zeros(ncells_ob)
    dpa_per_y_ob_std  = np.zeros(ncells_ob)

    for tally in dpa_tallies:
        dpa_tally = sp.get_tally(id=tally.id)

        # Find cell for that Tally
        cell_filter_dpa = next(
            f for f in dpa_tally.filters if isinstance(f, openmc.CellFilter)
        )
        cid = cell_filter_dpa.bins[0]  
        idx = cell_ids.index(cid)  # index in OB cell_ids list

        # Mean / Std per source neutron
        dpa_mean = dpa_tally.summation(
            nuclides=dpa_tally.nuclides
        ).mean.squeeze()

        dpa_std = dpa_tally.summation(
            nuclides=dpa_tally.nuclides
        ).std_dev.squeeze()

        # ----------------------------
        # Compute DPA/y mean and std
        # Formula: DPA = (0.8 / (2 Ed)) * dpa 
        # scale by neutron_source_rate and seconds_in_year #
        # uncertainty propagation: std_y = scale_factor * std
        #  ----------------------------
        Ed = materials.Ed(materials.element(dpa_tally.nuclides))

        # Scale to DPA per year 
        scale_factor = 0.8 / (2 * Ed) * neutron_source_rate * s_in_y

        dpa_y_mean_cell = scale_factor * dpa_mean
        dpa_y_std_cell  = scale_factor * dpa_std

        # Normalize by atoms in this cell (or solid_materials_atoms )
        atoms = cell_total_atoms_solid[cid]  # or cell_total_atoms[cid]

        dpa_per_y_ob_mean[idx] += dpa_y_mean_cell / atoms
        dpa_per_y_ob_std[idx]  += dpa_y_std_cell  / atoms

    dpa_per_y_ob_mean = np.asarray(dpa_per_y_ob_mean)
    dpa_per_y_ob_std  = np.asarray(dpa_per_y_ob_std)

    print("Maximum DPA/y (OB): ", np.max(dpa_per_y_ob_mean))

    plt.figure()
    plt.yscale('log')

    plt.step(
        x_ob,
        dpa_per_y_ob_mean,
        where='mid',           
        color='k',
        marker='o',
        markersize=1.5,
        label='DPA/y (OB)'
    )

    plt.fill_between(
        x_ob,
        np.maximum(dpa_per_y_ob_mean - dpa_per_y_ob_std, 1e-30),
        dpa_per_y_ob_mean + dpa_per_y_ob_std,
        step='mid',            
        alpha=0.3,
        color='gray',
        label='Std. Dev.'
    )

    plt.grid()
    plt.ylabel('DPA/y')
    plt.xlabel('Radial Position [cm]')
    plt.legend()
    plt.savefig('dpa_OB.png')
    plt.close()


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
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'materials'))
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
RS = materials.concrete(2.35) # PNNL COMPENDIUM
RS.name = 'RS'
# Divertor (From Bluemira Model)
divertor = openmc.Material.mix_materials([ccz, c, eurofer, t, water], [0.00552, 0.00438, 0.5238,  0.01026, 0.45604], 'vo',name="Divertor")

# --------------------------------------------------------------------
# homogenization from https://doi.org/10.1016/j.fusengdes.2020.111833
# Breeder Design WCLL
# First Wall
pf = openmc.Material.mix_materials([t, water, eurofer],[0.0027, 0.14268, 0.85462],'vo',name='FW')

# Inner Breeder Blanket
matIB1 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.835509, 0.025075, 0.139416], 'vo',name="IB1")
matIB2 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.858, 0.018, 0.124], 'vo',name="IB2")
matIB3 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.8132, 0.0158, 0.171], 'vo',name="IB3")
matIB4 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.426574, 0.015984, 0.557442], 'vo',name="IB4")
matIB5 = openmc.Material.mix_materials([water, eurofer], [0.486, 0.514], 'vo',name="IB5")
# Outer Breeder Blanket
matOB1 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.835509, 0.025075, 0.139416], 'vo',name="OB1")
matOB2 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.858, 0.018, 0.124], 'vo',name="OB2")
matOB3 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.8132, 0.0158, 0.171], 'vo',name="OB3")
matOB4 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.8132, 0.0158, 0.171], 'vo',name="OB4")
matOB5 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.8132, 0.0158, 0.171], 'vo',name="OB5")
matOB6 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.8132, 0.0158, 0.171], 'vo',name="OB6")
matOB7 = openmc.Material.mix_materials([PbLi, water, eurofer], [0.426574, 0.015984, 0.557442], 'vo',name="OB7")
matOB8 = openmc.Material.mix_materials([water, eurofer], [0.486, 0.514], 'vo',name="OB8")
# --------------------------------------------------------------------

# Vacuum Vessel
VV = openmc.Material.mix_materials([ss316, water],[0.6, 0.4],'vo',name='VV')

# Port fillings 
portf = openmc.Material.mix_materials([ss316, water], [0.6, 0.4], 'vo', name='PF')

# Poloidal coils (From Bluemira Model)
PC = openmc.Material.mix_materials([nbti, c, epoxy, bronze, h, ss316], [0.02895, 0.1169, 0.18, 0.0735, 0.1682, 0.43245], 'vo',name="PC")

# Toroidal Coils (From Bluemira Model)
TFcoil = openmc.Material.mix_materials([nb3sn, c, epoxy, bronze, h, ss316], [0.02895, 0.1169, 0.18,  0.0735, 0.1682, 0.43245], 'vo',name="TFC")

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
model.settings.particles = 1_000
model.settings.run_mode = "fixed source"
model.settings.source = my_source
#model.settings.statepoint_batch = [20, 40, 60, 80, 100]

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

# reclaim the ID space for all cells and surfaces currently in the model
# to avoid clashes later in the model
openmc.reserve_ids([c_id for c_id in model.geometry.get_all_cells()], cls=openmc.Cell)
openmc.reserve_ids([s_id for s_id in model.geometry.get_all_surfaces()], cls=openmc.Surface)

# -------------

# SURFACES
class Orientation(Enum):
    FORWARD = 1 # indicates normal vector points outwards from the volume
    REVERSE = -1 # indicates normal vector points inwards to the volume

def dagmc_bounding_box(pydagmc_model, volume_id):
    """Returns the bounding box of a given volume in a PyDAGMC model.

    Args:
        pydagmc_model: PyDAGMC model object
        volume_id: ID of the volume for which to find the bounding box

    Returns:
        openmc.BoundingBox: Bounding box of the volume
    """
    volume = pydagmc_model.volumes_by_id[volume_id]
    triangle_coords = volume.triangle_coords
    min_coords = np.min(triangle_coords, axis=0)
    max_coords = np.max(triangle_coords, axis=0)
    return openmc.BoundingBox(min_coords, max_coords)

# Updates on Patrick's version of dagmc_surface functions
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


dagmc_universe_cells = dagmc_universe.get_all_cells()
for volume in pydagmc_model.volumes:
    dagmc_universe_cells[volume.id].volume = volume.volume
    dagmc_universe_cells[volume.id].bounding_box = dagmc_bounding_box(pydagmc_model, volume.id)

# list of equatorial cells in Bluemira_11_10d.h5m model. (will change with other model)
# OB - outer blankets  (IB version can provided, currently removed the IB1 for depletion calculations)
cell_list_ob_1 = [49, 45, 79, 66, 82, 94, 103, 112, 124, 132, 2]
cell_ids = cell_list_ob_1

# centroids of layers
# if outer blanket (Armor, FW, OB1, OB2, OB3, OB4, OB5, OB6, OB7, OB8, VV)
xcentroids_ob = (0.1, 1.1, 6.0, 15.0, 25.0, 35.0, 45.0, 55.0, 70.0, 90.0, 157.2)

# Obtain a dict and surface list (for further albedo calculations)
info, all_surface_ids = dagmc_volume_surface_info(pydagmc_model, cell_ids)

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

# Total flux
flux_tally_total = openmc.Tally()
flux_tally_total.filters = [cell_filter, particle_filter]
flux_tally_total.scores = ['flux']
model.tallies.append(flux_tally_total)

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

    p_current_tally = openmc.Tally()
    p_current_tally.filters = [cell_from_filter, surf_filter, n_particle_filter]
    p_current_tally.scores  = ["current"]

    model.tallies.append(p_current_tally)
    p_current_tallies[cid] = p_current_tally

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
cell_atomic_ratio = {}

# Obtain cell_total_atoms and cell_total_atoms_solid
for cid in cell_ids:
    cell = all_cells[cid]
    all_nuclides.update(cell.fill.get_nuclides())
    
    cell.nuclide_atoms = cell.fill.get_nuclide_atoms(volume=cell.volume)
    cell.total_atoms = sum(cell.nuclide_atoms.values())

    cell_nuclide_atoms[cid] = cell.nuclide_atoms
    cell_total_atoms[cid]   = cell.total_atoms

    # Per-nuclide atomic fractions
    # Added so I can compare nuclide with depletion results
    cell_atomic_ratio[cid] = {
        nuc: atoms / cell_total_atoms[cid]
        for nuc, atoms in cell_nuclide_atoms[cid].items()
    }
    
    # total sum of atoms in solid nuclides in cell [cid] 
    solid_total_atoms = sum(
        atoms for nuc, atoms in cell.nuclide_atoms.items()
        if nuc in solid_nuclides
    )
    cell_total_atoms_solid[cid] = solid_total_atoms
    ratio_solid = cell_total_atoms_solid[cid] / cell_total_atoms[cid]


# CREATE DPA TALLIES
# Current version is generating tallies for each element existent in their respective cell.
# Also filters (maybe not the right word) to the elements present in the solid_nuclides list.
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
# built-in modules
from collections import defaultdict
from enum import Enum
import math
import re

import openmc
import numpy as np
import matplotlib.pyplot as plt
import materials
from openmc_plasma_source import tokamak_source
import pandas as pd
import pydagmc

_DAGMC_MODEL_FILE = 'bluemira_10_22.h5m'

model = openmc.Model()

# --------------------------------
# GEOMETRY
# --------------------------------

# Import DAGMC geometry
dagmc_universe = openmc.DAGMCUniverse(filename='bluemira_10_22.h5m')

pydagmc_model = pydagmc.Model(str(dagmc_universe.filename))

# reserve volume and surface IDs in the DAGMC model to avoid overlaps
openmc.reserve_ids([v.id for v in pydagmc_model.volumes], cls=openmc.Cell)
openmc.reserve_ids([s.id for s in pydagmc_model.surfaces], cls=openmc.Surface)

# Create planes for reflective BC
def azimuthal_plane(theta_deg, boundary_type=None, name=None, surface_id=None):
    theta = math.radians(theta_deg)
    a = -math.tan(theta)
    b = 1.0
    c = 0.0
    d = 0.0
    return openmc.Plane(a=a, b=b, c=c, d=d, boundary_type=boundary_type,
                        name=name, surface_id=surface_id)

# Angles for planes
theta0 = 0.0            # starting angle in degrees
theta1 = theta0 + 22.5  # 1/16 sector

# Reflective azimuthal cuts
cut_lo = azimuthal_plane(theta0, boundary_type='reflective', name='phi_lo')
cut_hi = azimuthal_plane(theta1, boundary_type='reflective', name='phi_hi')

# Outer BC graveyard (top/bottom)
z_min = openmc.ZPlane(z0=-2500, boundary_type='vacuum', name='z_min')
z_max = openmc.ZPlane(z0=+1500, boundary_type='vacuum', name='z_max')

# Outer BC graveyard (radial)
r_out = openmc.ZCylinder(r=2500, boundary_type='vacuum', name='r_out')
sector_region = (+cut_lo & -cut_hi) & (+z_min & -z_max) & (-r_out)

# Put the DAGMC universe inside the CSG cell
sector_cell = openmc.Cell(cell_id=10_000, region=sector_region, fill=dagmc_universe, name='sector_container')
model.geometry = openmc.Geometry(root=[sector_cell])

# --------------------------------
#  Materials
# --------------------------------
# Defined these materials in materials.py (still missing DE)
# This list will change depending on the materials/mixtures
# Rough estimate of density
ss316 = materials.ss316(8.0) # ss316 (for mixing materials)
ccz = materials.CuCrZr(8.9) #CuCrZr
c = materials.Cu(8.92) #Cu
t = materials.W(19.3) #Tungsten
water = materials.Water(1.0) #H20
nb3sn = materials.Nb3Sn(5.7) #Nb3Sn
epoxy = materials.Epoxy(1.207) #Epoxy resin
bronze = materials.Bronze(8.8775) #Bronze
h = materials.Helium(0.1785) #Helium
be = materials.Be(1.85) #Beryllium
li4si04 = materials.Li4SiO4(2.4) #Li4SiO4
nbti = materials.NbTi(6.538) #NbTi
be12ti = materials.be12ti(2.4) #Be12Ti
kalos = materials.kalos_cb(2.4) #kalos_cb
eurofer = materials.eurofer97(7.87)

# ATRIBUTE MATERIALS TO DAGMC geometry
# Mateirals with simple structure(repetitive - simplify this ?)
# Vacuum Vessel (VV)
VV = materials.ss316(8.0)
VV.name = 'VV'
# Radiation shielding
RS = materials.ss316(8.0)
RS.name = 'RS'
# Cryostat
CR = materials.ss316(8.0)
CR.name = 'CR'

# DEFINE MIXTURES
# There is a certain amount of void in each of them (account for He)
# Several references use Void as a substitute for He cooling elements
divertor = openmc.Material.mix_materials([ccz, c, ss316, t, water], [0.00552, 0.00438, 0.5238,  0.01026, 0.45604], 'vo',name="Divertor")
TFcoil = openmc.Material.mix_materials([nb3sn, c, epoxy, bronze, h, ss316], [0.02895, 0.1169, 0.18,  0.0735, 0.1682, 0.4319], 'vo',name="TFC")
breeder1 = openmc.Material.mix_materials([be12ti, kalos, eurofer], [0.555, 0.088, 0.1230], 'vo',name="IB")
breeder2 = openmc.Material.mix_materials([be12ti, kalos, eurofer], [0.493, 0.1030, 0.1430], 'vo',name="OB")
PC = openmc.Material.mix_materials([nbti, c, epoxy, bronze, h, ss316], [0.02895, 0.1169, 0.18, 0.0735, 0.1682, 0.4319], 'vo',name="PC")
portf = openmc.Material.mix_materials([ss316, water], [0.6, 0.4], 'vo', name='PF')
pf = openmc.Material.mix_materials([eurofer],[0.5730],'vo',name='FW')

# Build materials model
model.materials = openmc.Materials([pf, portf, VV,divertor, TFcoil, breeder1, breeder2, PC, CR, RS])

# --------------------------------
#  SOURCE
# --------------------------------
# Openmc_plasma_source
my_source = tokamak_source(
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
model.settings.particles = 10000
model.settings.run_mode = "fixed source"
model.settings.source = my_source

# --------------------------------
# VOLUME CALCULATION
# --------------------------------

# apply volumes from PyDAGMC to the OpenMC model cells
model.init_lib()
model.sync_dagmc_universes()
model.finalize_lib()

dagmc_universe_cells = dagmc_universe.get_all_cells()
for volume in pydagmc_model.volumes:
    dagmc_universe_cells[volume.id].volume = volume.volume

# list of equatorial cells in Bluemira_10_22.h5m model. (will change with other model)
# OB - outer blankets # IB - inner blankets
# Mostly focusing now on Cell_list_ob_1
cell_list_ob_1 = [99, 47, 51, 49, 53, 3]
cell_list_ob_2 = [105, 62, 66, 64, 68]
cell_list_ob_3 = [102, 55, 59, 57, 61, 5]
cell_list_ib_4 = [101, 21, 25, 19, 23, 4]
cell_list_ib_5 = [107, 29, 33, 27, 31, 6]


def dagmc_adjacent_cells(pydagmc_model, volume_id):
    """Returns all volumes adjacent to a given volume in a PyDAGMC model.

    Args:
        pydagmc_model: PyDAGMC model object
        volume_id: ID of the volume for which to find adjacent volumes
    """
    volume = pydagmc_model.volumes_by_id[volume_id]
    adj_volumes = set()
    for surface in volume.surfaces:
        adj_volumes.update(surface.volumes)

    return list(adj_volumes - {volume})


class Orientation(Enum):
    FORWARD = 1 # indicates normal vector points outwards from the volume
    REVERSE = -1 # indicates normal vector points inwards to the volume


def dagmc_surface_orientations(pydagmc_model, volume_id):
    """Returns a dictionary mapping surface IDs to their orientations
    with respect to a given volume in a PyDAGMC model.

    Args:
        pydagmc_model: PyDAGMC model object
        volume_id: ID of the volume for which to find surface orientations

    Returns:
        dict: Mapping of surface IDs to Orientation enum values
    """
    volume = pydagmc_model.volumes_by_id[volume_id]
    orientations = {}
    for surface in volume.surfaces:
        parent_volumes = surface.senses
        if parent_volumes[0].id == volume_id:
            orientations[surface.id] = Orientation.FORWARD
        else:
            orientations[surface.id] = Orientation.REVERSE

    return orientations


# --------------------------------
# TALLY
# --------------------------------
# Source intensity -> normalization
total_power = 2e9 # 2000 MW
section = 1/16
section_power = total_power * section
ev_to_joule = 1.60218e-19
ev_fission = 17.6e6 # eV
convert_e = ev_to_joule * ev_fission
neutron_source_rate = section_power / convert_e

# centroids of layers
# (FW, breeder1, breeder2, breeder3, breeder4, VV)
xcentroids = (1.0, 12.0, 32.0, 50.375, 78.5, 150.0)

# FILTERS
cell_filter = openmc.CellFilter(cell_list_ob_1)
n_particle_filter = openmc.ParticleFilter(bins=['neutron'])
p_particle_filter = openmc.ParticleFilter(bins=['photon'])
energies = openmc.mgxs.GROUP_STRUCTURES['CCFE-709']
energy_filter = openmc.EnergyFilter(energies)
unit_lethargy = [np.log(energies[i+1]/energies[i]) for i in range(len(energies)-1)]


# Filter solid/fluid materials to account correctly gas production in solids.
# Ideas:
# Material filters  -> Works on MATERIALS (full mixture)
#                   -> Don't work with specific elements inside mixtures

# Nuclide filters   -> Works with any nuclide inside materials/mixture
#                   -> Be careful with overlapping nuclides from solid and fluid.

# Also - greater detail in model to consider cooling channels separately

# add tallies for solution quantities
model.tallies = openmc.Tallies()

# Tallies
# flux score, with an energy filter
# neutron flux
n_flux_tally = openmc.Tally(name='n_flux_tally')
n_flux_tally.filters = [cell_filter, n_particle_filter, energy_filter]
n_flux_tally.scores = ['flux']
model.tallies.append(n_flux_tally)

# gamma flux
p_flux_tally = openmc.Tally(name='p_flux_tally')
p_flux_tally.filters = [cell_filter, p_particle_filter, energy_filter]
p_flux_tally.scores = ['flux']
model.tallies.append(p_flux_tally)

# heating
heating_tally = openmc.Tally(name='heating_tally')
heating_tally.filters = [cell_filter]
heating_tally.scores = ['heating']
model.tallies.append(heating_tally)

# proton(hydrogen) production tallies
h1_tally = openmc.Tally(name='H1_tally')
h1_tally.filters = [cell_filter] #,material_solid_filter]
h1_tally.scores = ['H1-production']
model.tallies.append(h1_tally)

# helium production tallies
he3_tally = openmc.Tally(name='He3_tally')
he3_tally.filters = [cell_filter] #,material_solid_filter]
he3_tally.scores = ['He3-production']
model.tallies.append(he3_tally)

he4_tally = openmc.Tally(name='He4_tally')
he4_tally.filters = [cell_filter] #,material_solid_filter]
he4_tally.scores = ['He4-production']
model.tallies.append(he4_tally)

# DPA tallies
# Moved DPA tallies for later where I load the DAGMC file

# -------------------
# Plot
# -------------------
plot = openmc.Plot()
plot.basis = 'xz'
plot.origin = (1000.0, 5.0, 0.0)              # x=0 slice
plot.width = (2000.0, 2000.0)            # 20 m x 20 m (in cm)
plot.pixels = (1000, 1000)
plot.color_by = 'cell'
plot.filename = 'xz_geometry'              # will create yz_geometry.png (or .png on some builds)

plots = openmc.Plots([plot])
model.plots = plots

# --------------------------------
# Load Information (save volumes/ num of atoms per cell inside DAGMC cells)
# --------------------------------
# Allocate volume calculations and atoms per cell
all_cells = model.geometry.get_all_cells()

# Collect all nuclides over selected cells
all_nuclides = set()
for c in cell_list_ob_1:
    cell = all_cells[c]
    all_nuclides.update(cell.fill.get_nuclides())

# group nuclides per element
nuclides_by_element = {}
for n in sorted(all_nuclides):
    el = re.sub(r'\d+', '', n)
    nuclides_by_element.setdefault(el, []).append(n)

# DPA tallies
dpa_tallies = []
for nuc_list in nuclides_by_element.values():
    if not nuc_list:
        continue
    dpa_tally = openmc.Tally()
    dpa_tally.scores = ["damage-energy"]
    dpa_tally.filters = [cell_filter]
    dpa_tally.nuclides = nuc_list
    model.tallies.append(dpa_tally)
    dpa_tallies.append(dpa_tally)

# Export tallies.xml
model.tallies.export_to_xml()

# Atoms calculatiosn
cell_nuclide_atoms = {}   # {cell_id: {nuclide: mean_atoms}}
cell_total_atoms   = {}   # {cell_id: mean_total_atoms}

# VOLUME CELL UPDATE
for cid in cell_list_ob_1:
    cell = all_cells[cid]
    print(cell)
    # Attach to the Cell object (for runtime only)
    cell.nuclide_atoms = cell.fill.get_nuclide_atoms(volume=cell.volume)
    cell.total_atoms = sum(cell.nuclide_atoms.values())
    # save for scaling in post-processing
    cell_nuclide_atoms[cid] = cell.nuclide_atoms
    cell_total_atoms[cid]   = cell.total_atoms

model.finalize_lib()

# --------------------------------
# EXECUTE (and plot)
# --------------------------------
statepoint = model.run(threads=4)
with openmc.StatePoint(statepoint) as sp:

    # Define the scaling for each cell.volume in the tally
    scaling_by_cell = {}
    for cid in cell_list_ob_1:
        cell_cid = all_cells[cid]
        vol = cell_cid.volume
        scaling = 1.0 / vol * neutron_source_rate
        scaling_by_cell[cid] = scaling

    # Neutron Spectrum
    n_tally = sp.get_tally(name="n_flux_tally")
    neutron_flux = n_tally.get_reshaped_data()
    neutron_flux_std_dev = n_tally.get_reshaped_data(value='std_dev')

    # csv file for students
    E_mid = 0.5 * (energies[:-1] + energies[1:])
    flux_data = {'Energy [eV]': E_mid}

    for i, cid in enumerate(cell_list_ob_1):
        scaling = scaling_by_cell[cid]
        flux_scaled = neutron_flux[i].flatten() * scaling / unit_lethargy
        flux_data[f'Flux_cell_{cid} [n/cm^2/s/lethargy]'] = flux_scaled
        plt.loglog(
            energies[:-1],
            flux_scaled,
            label=f"Depth = {xcentroids[i]:.2f} cm"
        )

    plt.legend()
    plt.grid()
    plt.ylabel('Neutron Flux Per Unit Lethargy [1/cm$^2$/s]')
    plt.xlabel('Energy [eV]')
    plt.xlim([1, 100e6])
    plt.savefig('n_flux_spectrum.png')
    plt.close()

    # save the CSV file
    df_flux = pd.DataFrame(flux_data)
    df_flux.to_csv('neutron_flux_spectrum.csv', index=False)

    # Photon Spectrum
    p_tally = sp.get_tally(name='p_flux_tally')
    photon_flux = p_tally.get_reshaped_data()

    for i, cid in enumerate(cell_list_ob_1):
        scaling = scaling_by_cell[cid]
        plt.loglog(
            energies[:-1],
            photon_flux[i].flatten() * scaling / unit_lethargy,
            label=f"Depth = {xcentroids[i]:.2f} cm"
        )

    plt.legend()
    plt.grid()
    plt.ylabel('Photon Flux [1/cm$^2$/s/eV]')
    plt.xlabel('Energy [eV]')
    plt.xlim([1, 100e6])
    plt.savefig('p_flux_spectrum.png')
    plt.close()

    # Total fluxes by integrating over energy
    # Neutron flux
    n_cells = len(cell_list_ob_1)
    total_neutron_flux = np.zeros(n_cells)
    for i, cid in enumerate(cell_list_ob_1):
        scaling = scaling_by_cell[cid]
        total_neutron_flux[i] = np.sum(neutron_flux[i].flatten()) * scaling

    plt.semilogy(xcentroids, total_neutron_flux, label='Neutron flux')
    # Gamma flux
    total_photon_flux = np.zeros(n_cells)
    for i, cid in enumerate(cell_list_ob_1):
        scaling = scaling_by_cell[cid]
        total_photon_flux[i] = np.sum(photon_flux[i].flatten()) * scaling

    plt.semilogy(xcentroids, total_photon_flux, label='Photon flux')
    plt.legend()
    plt.grid()
    plt.ylabel('Flux [1/cm$^2$/s]')
    plt.xlabel('Radial Position [cm]')
    plt.savefig('flux.png')
    plt.close()

    # heating
    h_tally = sp.get_tally(name='heating_tally')
    heating = h_tally.get_values().flatten()

    heat = []
    for i, cid in enumerate(cell_list_ob_1):
        scaling = scaling_by_cell[cid]
        heating_per_cc = heating[i] * scaling * ev_to_joule
        heat.append(heating_per_cc)

    plt.semilogy(xcentroids, heat, marker='o', color='k', markersize=1.5)
    plt.grid()
    plt.ylabel('Heating [W/cm^3]')
    plt.xlabel('Radial Position [cm]')
    plt.savefig('heating.png')
    plt.close()

    # Hydrogen production - No material filter yet
    hydrogen1_tally = sp.get_tally(name='H1_tally')
    h1 = hydrogen1_tally.get_values().flatten()

    h = []
    for i, cid in enumerate(cell_list_ob_1):
        h_per_s = h1[i] * neutron_source_rate
        h_per_y = h_per_s * (365 * 24 * 60 * 60)
        h_appm_per_y = h_per_y / cell_total_atoms[cid] * 1e6
        h.append(h_appm_per_y)

    print('Maximum proton appm/y: ', np.max(h))

    plt.semilogy(xcentroids, h, marker='o', color='k', markersize=1.5)
    plt.grid()
    plt.ylabel('Proton [appm/y]')
    plt.xlabel('Radial Position [cm]')
    plt.savefig('h1.png')
    plt.close()

    # Helium production - No material filter yet
    helium3_tally = sp.get_tally(name='He3_tally')
    he3 = helium3_tally.get_values().flatten()
    helium4_tally = sp.get_tally(name='He4_tally')
    he4 = helium4_tally.get_values().flatten()

    he = []
    for i, cid in enumerate(cell_list_ob_1):
        he_per_s = (he3[i] + he4[i]) * neutron_source_rate
        he_per_y = he_per_s * (365 * 24 * 60 * 60)
        he_appm_per_y = he_per_y / cell_total_atoms[cid] * 1e6
        he.append(he_appm_per_y)

    print('Maximum helium appm/y: ', np.max(he))

    plt.semilogy(xcentroids, he, marker='o', color='k', markersize=1.5)
    plt.grid()
    plt.ylabel('Helium [appm/y]')
    plt.xlabel('Radial Position [cm]')
    plt.savefig('he.png')
    plt.close()

    # create radial plots of the dpa; each of the tallies is for a particular element
    dpa_per_y = np.zeros(n_cells)
    atoms_vec = np.array([cell_total_atoms[cid] for cid in cell_list_ob_1], dtype=float)

    for tally in dpa_tallies:
        dpa_tally = sp.get_tally(id=tally.id)

        dpa = dpa_tally.summation(nuclides=dpa_tally.nuclides).mean.flatten()

        # get the Ed for this element
        Ed = materials.Ed(materials.element(dpa_tally.nuclides))

        displacements_per_source = 0.8 * dpa / (2 * Ed)
        displacements_per_s = displacements_per_source * neutron_source_rate
        displacements_per_y = displacements_per_s * (365 * 24 * 60 * 60)

        for i, cid in enumerate(cell_list_ob_1):
            dpa_per_y[i] += displacements_per_y[i] / float(cell_total_atoms[cid])

    print('Maximum dpa: ', np.max(dpa_per_y))

    xcentroids = np.asarray(xcentroids).squeeze()
    dpa_per_y  = np.asarray(dpa_per_y, float).squeeze()

    plt.semilogy(xcentroids, dpa_per_y, marker='o', color='k', markersize=1.5)
    plt.grid()
    plt.ylabel('DPA/y')
    plt.xlabel('Radial Position [cm]')
    plt.savefig('dpa.png')
    plt.close()


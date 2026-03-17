import matplotlib.pyplot as plt
import math
import numpy as np

from openmc_plasma_source import (
    tokamak_convert_a_alpha_to_R_Z,
    tokamak_ion_density,
    tokamak_ion_temperature,
    tokamak_neutron_source_density,
)

sample_size = 20000
mode = "H"

# [1] https://iopscience.iop.org/article/10.1088/1741-4326/ada6d9/pdf
# [2] Fausser

# from OUT.DAT
minor_radius = 300.2

# from OUT.DAT
major_radius = 840.67

elongation = 1.739

# [1], table 2, \alpha_n
# ion_density_peaking_factor: Ion density peaking factor, (referred in [2] as ion density exponent, \alpha_n where it is shown a value of 1.0 for an A-mode plasma)
ion_density_peaking_factor = 1

ion_density_separatrix = 1.02e19
ion_density_centre = 6.8e19 # big impact
ion_density_pedestal = 5.78e19

# OUT.DAT lists Electron temperature on axis (keV) as 23.775 keV; based on [1] the ion and electron temperatures are the same
ion_temperature_centre = 23.775e3

# OUT.DAT lists electron temperature at the pedestal height as 5.5 keV (teped)
ion_temperature_pedestal = 5.5e3

# [1], page 10: "Ion temperature is assumed to be equal to electron temperature... The electron temperature at the separatrix is taken as 0.1 keV"
# OUT.DAT lists electron temperature at the separatrix as 0.1 keV (tesep)
ion_temperature_separatrix = 0.1e3

# [1], T-profile exponent in Table 2 has 1.45
# OUT.DAT has 1.45 for a quantity alphat, the "Temperature profile index"
# ion_temperature_peaking_factor: Ion temperature peaking factor (referred in [2] as ion temperature exponent alpha_T). Value of 8.06 is given for example of A-mode plasma.
#ion_temperature_peaking_factor = 8.06
ion_temperature_peaking_factor = 1.45

# [1] has 2.0 in Table 1
# value of 6 comes from the github example
#ion_temperature_beta = 6
ion_temperature_beta = 2.0

# TODO: do we need to take at plasma edge, or at the 95th percentile flux?
triangularity = 0.333

shafranov_factor = 0.44789

pedestal_radius = 0.94 * minor_radius

# create a sample of (a, alpha) coordinates
a = np.random.random(sample_size) * minor_radius
alpha = np.random.random(sample_size) * 2 * np.pi

temperatures = tokamak_ion_temperature(
    r=a,
    ion_temperature_centre=ion_temperature_centre/1e3,
    ion_temperature_pedestal=ion_temperature_pedestal/1e3,
    ion_temperature_separatrix=ion_temperature_separatrix/1e3,
    ion_temperature_peaking_factor=ion_temperature_peaking_factor,
    ion_temperature_beta=ion_temperature_beta,
    major_radius=major_radius,
    pedestal_radius=pedestal_radius,
    mode=mode
)

densities = tokamak_ion_density(
    r=a,
    ion_density_centre=ion_density_centre,
    ion_density_pedestal=ion_density_pedestal,
    ion_density_peaking_factor=ion_density_peaking_factor,
    ion_density_separatrix=ion_density_separatrix,
    major_radius=major_radius,
    pedestal_radius=pedestal_radius,
    mode=mode
)

neutron_source_density = tokamak_neutron_source_density(densities, temperatures, "DD")

RZ = tokamak_convert_a_alpha_to_R_Z(
    a=a,
    alpha=alpha,
    elongation=elongation,
    major_radius=major_radius,
    minor_radius=minor_radius,
    shafranov_factor=shafranov_factor,
    triangularity=triangularity
)

plt.scatter(RZ[0], RZ[1], c=neutron_source_density, s=1)
plt.ylim([-1000, 1000])
plt.xlim([0, 2000])
plt.gca().set_aspect("equal")
plt.xlabel("R [cm]")
plt.ylabel("Z [cm]")
plt.colorbar(label="neutron source density")

plt.savefig("tokamak_source_neutron_source_density.png")
print("written tokamak_source_neutron_source_density.png")

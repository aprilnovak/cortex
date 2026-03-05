import openmc

# Script created by April to check the fractions used to screen out non-solid contributions
# to dpa, H, and He production

# check for the first wall, looks reasonable
eurofer_density = 7.87
water_density = 0.866
eurofer_vol_fraction = 0.85462
water_vol_fraction = 0.14268
volume = 1.0
Na = 6.022e23

# CHANGE THIS LINE TO MATCH EUROFER97 definition in materials.py
eurofer_o_wt_percent = 0.01

oxygen_molar_mass = openmc.data.atomic_weight('O')
hydrogen_molar_mass = openmc.data.atomic_weight('H')
water_molar_mass = oxygen_molar_mass + 2 * hydrogen_molar_mass

g_eurofer = volume * eurofer_vol_fraction * eurofer_density
g_water = volume * water_vol_fraction * water_density

# certain percent of the eurofer is oxygen
g_O_in_eurofer = g_eurofer * eurofer_o_wt_percent * 1e-2
mol_O_in_eurofer = g_O_in_eurofer / oxygen_molar_mass
mol_O16_in_eurofer = mol_O_in_eurofer * openmc.data.isotopes('O')[0][1]
atoms_O16_in_eurofer = mol_O16_in_eurofer * Na

mol_O_in_water = g_water / water_molar_mass
mol_O16_in_water = mol_O_in_water * openmc.data.isotopes('O')[0][1]
atoms_O16_in_water = mol_O16_in_water * Na

print('Fraction for O16: ', atoms_O16_in_eurofer / (atoms_O16_in_water + atoms_O16_in_eurofer))

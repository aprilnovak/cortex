import matplotlib.pyplot as plt
import csv
import openmc
import sys

# read the armor neutron energy spectrum from the neutronics_results. This script
# can only be run after neutronics_post.py has been run because it reads in the CSV files

E = []
flux = []
with open('../neutronics_results/OB_1_b6/neutron_spectrum_OB_1_b6.csv') as f:
  reader = csv.reader(f)
  next(reader)
  for row in reader:
    E.append(float(row[0]))
    flux.append(float(row[1]))

# now, read in tungsten cross section data to see which reaction(s) are visible here
path_to_my_cross_sections = '/gpfs/fs1/home/anovak/cross_sections/endfb-vii.1-hdf5'

nuclides = ['W182', 'W183', 'W184', 'W186']

T = '294K'

# [5] Reaction flag - these follow the ENDF-6 formatting convention
# (https://www.oecd-nea.org/dbdata/data/manual-endf/endf102.pdf). Some common
# reactions of interest:
#    1: total
#    2: elastic scatter
#    4: inelastic scatter
#   18: fission
#   27: total absorption
#  102: radiative capture
reactions = [2, 102]
labels = ['(n, n)', '(n, $\\gamma$)']
colors = ['b', 'mediumslateblue', 'r', 'salmon', 'green', 'limegreen', 'darkviolet', 'plum']
j = 0

for nn, nuclide in enumerate(nuclides):
  try:
    data = openmc.data.IncidentNeutron.from_hdf5(path_to_my_cross_sections + '/neutron/' + nuclide + '.h5')
  except Exception as e:
    print("The path to your cross sections, '", \
      path_to_my_cross_sections + "/neutron/" + nuclide + \
      ", was not found. Double check the location of your file download " + \
      "or whether you have a typo in the nuclide.")
    sys.exit()

  for i, r in enumerate(reactions):
    # get the cross section data
    xsr = data[r].xs[T]
    energies = data.energy[T]
    plt.loglog(energies, xsr(energies) * 1e10, label=nuclide + ' ' + labels[i], color=colors[j])
    j += 1

plt.loglog(E, flux, color='k', label='Flux', lw=2.5)
plt.grid()
plt.xlim([1e0,1e3])
plt.ylim([1e9, 1e16])
plt.ylabel('Neutron flux per unit lethargy [1/cm$^2$/s]')
plt.xlabel('Energy [ev]')
plt.legend(ncol=5, loc='upper center', columnspacing=0.7, fontsize=9, handlelength=0.7)
plt.savefig('../neutronics_results/OB_1_b6/n_spectrum_ss.png', bbox_inches="tight")
plt.close()

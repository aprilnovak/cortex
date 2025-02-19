import openmc
import openmc.data

# This file defines default materials to populate into existing models

def PbLi(li6_enrichment, density):
  """ Return an OpenMC material for PbLi. We first add the impurity concentrations
      from ITER recommendations [10.1016/j.nme.2022.101146, table 2] (at their
      maximum values). Then, of the remaining material, we take Li at 0.62 weight %
      and Pb at 99.38 weight % [10.1016/j.nme.2022.101146, page 2].

      Keyword arguments:
      li6_enrichment -- the Li6 enrichment in Li [atomic percent]
      density -- the overall density of the material [g/cc]
  """

  print('Adding PbLi material...')

  PbLi = openmc.Material()

  # Impurities are given in wppm (weight part per million)
  PbLi.add_element('Ag', 10e-6, 'wo')
  PbLi.add_element('Cu', 10e-6, 'wo')
  PbLi.add_element('Nb', 10e-6, 'wo')
  PbLi.add_element('Pd', 10e-6, 'wo')
  PbLi.add_element('Zn', 10e-6, 'wo')
  PbLi.add_element('Cr', 50e-6, 'wo')
  PbLi.add_element('Fe', 50e-6, 'wo')
  PbLi.add_element('Mn', 50e-6, 'wo')
  PbLi.add_element('Ni', 50e-6, 'wo')
  PbLi.add_element('Pd', 50e-6, 'wo')
  PbLi.add_element('V' , 50e-6, 'wo')
  PbLi.add_element('Al', 100e-6, 'wo')
  PbLi.add_element('Si', 100e-6, 'wo')
  PbLi.add_element('Bi', 200e-6, 'wo')
  PbLi.add_element('Sn', 200e-6, 'wo')
  PbLi.add_element('W', 200e-6, 'wo')

  weight_sum = 0
  for nuclide in PbLi.nuclides:
    weight_sum += nuclide.percent

  print('\tImpurities (weight %): ', weight_sum * 100)

  # TODO: I don't know if these papers are using weight or atomic percents.
  # Try to find which one. Here, I am assuming atomic.
  molar_mass_lithium = li6_enrichment * openmc.data.atomic_mass('Li6') + \
    (1 - li6_enrichment) * openmc.data.atomic_mass('Li7')

  # total Li concentration taken as 0.62% by weight of whatever remains after impurities
  li_weight_percent = 0.62e-2 * (1 - weight_sum)
  li6_grams_in_1_g_li = li6_enrichment * openmc.data.atomic_mass('Li6') / molar_mass_lithium
  li7_grams_in_1_g_li = 1.0 - li6_grams_in_1_g_li
  print('\tLithium (weight %):    ', li_weight_percent * 100)

  PbLi.add_nuclide('Li6', li6_grams_in_1_g_li * li_weight_percent, 'wo')
  weight_sum += li6_grams_in_1_g_li * li_weight_percent
  PbLi.add_nuclide('Li7', li7_grams_in_1_g_li * li_weight_percent, 'wo')
  weight_sum += li7_grams_in_1_g_li * li_weight_percent

  # total Pb concentration taken as 99.38% by weight of whatever remains after impurities
  PbLi.add_element('Pb', 0.9938 * (1 - weight_sum), 'wo')
  print('\tLithium (weight %):    ', (1 - weight_sum) * 100)

  PbLi.set_density('g/cc', density)


import openmc
import openmc.data
import re

# This file defines default materials to populate into existing models

def PbLi(li6_enrichment, density):
  """ Return an OpenMC material for PbLi. We first add the impurity concentrations
      from ITER recommendations [10.1016/j.nme.2022.101146, table 2] (at their
      maximum values). Then, of the remaining material, we take Li at 0.62 weight %
      and Pb at 99.38 weight % [10.1016/j.nme.2022.101146, page 2].

      We assume that this material is not depletable, since it will be flowing on
      a timescale much faster than depletion.

      TODO: non-metallic impurities (O, C, N, H) are not included, even though
      some will remain when ingots are melted. Do we need to include these?

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
  print('\tLead (weight %):    ', (1 - weight_sum) * 100)

  PbLi.set_density('g/cc', density)
  PbLi.depletable = False
  return PbLi

def tungsten(density):
  """ Return an OpenMC material for Tungsten. TODO: add impurities
  """
  tungsten = openmc.Material()
  tungsten.add_element('W', 1.0, 'wo')
  tungsten.set_density('g/cc', density)
  return tungsten

def ods(density):
  """ Return an OpenMC material for oxide dispersion strengthened steel
      (10.1088/1741-4326/ac2523). This paper gives actual composition from a demonstration. Composition is given in weight %.
  """
  ods = openmc.Material()
  ods.add_element('Cr', 13.5, 'wo')
  ods.add_element('C', 0.0164, 'wo')
  ods.add_element('Mn', 0.0917, 'wo')
  ods.add_element('W', 0.99, 'wo')
  ods.add_element('N', 0.0078, 'wo')
  ods.add_element('O', 0.117, 'wo')
  ods.add_element('Y', 0.146, 'wo')
  ods.add_element('Ti', 0.16, 'wo')
  ods.add_element('Ni', 0.0599, 'wo')
  ods.add_element('H', 0.0057, 'wo')

  weight_sum = 0
  for nuclide in ods.nuclides:
    weight_sum += nuclide.percent

  print('\tAlloying elements (weight %): ', weight_sum)

  ods.add_element('Fe', 100 - weight_sum, 'wo')

  print('\tIron (weight %):             ', 100 - weight_sum)

  ods.set_density('g/cc', density)
  return ods

def eurofer97(density):
  """ Return an OpenMC material for Eurofer97 RAFM steel
      (10.1016/j.fusengdes.2018.06.027). This paper gives min, max, and target values
      for the main alloying elements, and maximum values for impurities (radiologically
      undesirable quantities, which may or may not be the same as other impurities
      which are not undesirable, per se). Composition is given in weight %.

      The table in the reference includes min weight %, max weight %, and target
      weight %, but not all elements have all three numbers provided. We make the
      following assumption:
        - if a target % is given, we use it
        - if only a range, between min and max is given, we take the average
        - if only one numeric value is given, we use it
        - if the minimum is as-low-as-possible, we take the maximum given
        - for the As, Sn, Sb, and Zr (these are lumped together), we assume the 0.05 weight % is evenly split among these
  """
  eurofer97 = openmc.Material()
  eurofer97.add_element('C', 0.11, 'wo')
  eurofer97.add_element('Cr', 9.0, 'wo')
  eurofer97.add_element('W', 1.1, 'wo')
  eurofer97.add_element('Mn', 0.4, 'wo')
  eurofer97.add_element('V', (0.25 + 0.15)/2, 'wo')
  eurofer97.add_element('Ta', 0.12, 'wo')
  eurofer97.add_element('N', 0.03, 'wo')
  eurofer97.add_element('P', 0.005, 'wo')
  eurofer97.add_element('S', 0.005, 'wo')
  eurofer97.add_element('B', 0.002, 'wo')
  eurofer97.add_element('O', 0.01, 'wo')

  alloy_weight_sum = 0
  for nuclide in eurofer97.nuclides:
    alloy_weight_sum += nuclide.percent

  print('\tAlloying elements (weight %): ', alloy_weight_sum)

  # impurities
  eurofer97.add_element('Nb', 0.005, 'wo')
  eurofer97.add_element('Mo', 0.005, 'wo')
  eurofer97.add_element('Ni', 0.01, 'wo')
  eurofer97.add_element('Cu', 0.01, 'wo')
  eurofer97.add_element('Al', 0.01, 'wo')
  eurofer97.add_element('Ti', 0.02, 'wo')
  eurofer97.add_element('Si', 0.05, 'wo')
  eurofer97.add_element('Co', 0.01, 'wo')
  eurofer97.add_element('As', 0.05/4, 'wo')
  eurofer97.add_element('Sn', 0.05/4, 'wo')
  eurofer97.add_element('Sb', 0.05/4, 'wo')
  eurofer97.add_element('Zr', 0.05/4, 'wo')

  weight_sum = 0
  for nuclide in eurofer97.nuclides:
    weight_sum += nuclide.percent

  print('\tImpurity elements (weight %): ', weight_sum - alloy_weight_sum)

  eurofer97.add_element('Fe', 100 - weight_sum, 'wo')

  print('\tIron (weight %):             ', 100 - weight_sum)

  eurofer97.set_density('g/cc', density)
  return eurofer97

def water(density):
  """ Return an OpenMC material for water. TODO: add impurities
  """

  water = openmc.Material()
  water.add_element('H', 2.0)
  water.add_element('O', 1.0)
  water.set_density('g/cc', density)
  return water

def helium(density):
  """ Return an OpenMC material for helium. TODO: add impurities
  """

  helium = openmc.Material()
  helium.add_element('He', 1.0)
  helium.set_density('g/cc', density)
  return helium

def W(density):
  """ Return an OpenMC material for pure tungsten.
  """

  tungsten = openmc.Material()
  tungsten.add_element('W', 100, 'wo')
  tungsten.set_density('g/cc', density)
  return tungsten

def atoms(material):
  Na = 6.023e23

  atoms = 0
  for n in material.nuclides:
    nuclide_mass = material.get_mass(nuclide=n.name)
    atoms += nuclide_mass / openmc.data.atomic_mass(n.name) * Na

  return atoms

def atoms_of_each_element(material):
  Na = 6.023e23

  atoms_of_element = {}
  for n in material.nuclides:
    nuclide_mass = material.get_mass(nuclide=n.name)
    element = re.sub(r'[0-9]+', '', n.name)
    n_atoms = nuclide_mass / openmc.data.atomic_mass(n.name) * Na
    if (element in atoms_of_element):
      atoms_of_element[element] += n_atoms
    else:
      atoms_of_element[element] = n_atoms

  return atoms_of_element

def nuclides_for_each_element(material):
  nuclides = {}
  # get the nuclide names corresponding to each element in a material
  for n in material.nuclides:
    element = re.sub(r'[0-9]+', '', n.name)
    if (element in nuclides):
      nuclides[element].append(n.name)
    else:
      nuclides[element] = []
      nuclides[element].append(n.name)

  return nuclides

def Ed(element):
  ed = {}
  ed['Be'] = 31
  ed['C'] = 31
  ed['Mg'] = 25
  ed['Al'] = 27
  ed['Si'] = 25
  ed['Ca'] = 40
  ed['Ti'] = 40
  ed['V'] = 40
  ed['Cr'] = 40
  ed['Mn'] = 40
  ed['Fe'] = 40
  ed['Co'] = 40
  ed['Ni'] = 40
  ed['Cu'] = 40
  ed['Zr'] = 40
  ed['Nb'] = 40
  ed['Mo'] = 60
  ed['Ag'] = 60
  ed['Ta'] = 90
  ed['W'] = 90
  ed['Au'] = 30
  ed['Pb'] = 25

  default = 25
  if element in ed:
    return ed[element]
  else:
    return default

def element(nuclides):
  elem = 'UNKNOWN'
  for n in nuclides:
    result = re.sub(r'\d+', '', n)
    if (elem == 'UNKNOWN'):
      elem = result
    elif (result != elem):
      raise ValueError("Nuclide list not a single nuclide!")
  return elem


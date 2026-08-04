"""
Stainless Steel Uranus 65 (X1CrNi25-21)

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
4.2571. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Stainless Steel Uranus 65 (X1CrNi25-21)'

density = 7.9

composition = {
    'C': {'wt': 0.02, 'type': 'alloy'},
    'Cr': {'wt': 25.0, 'type': 'alloy'},
    'Mn': {'wt': 2.0, 'type': 'alloy'},
    'Mo': {'wt': 0.3, 'type': 'alloy'},
    'Ni': {'wt': 20.5, 'type': 'alloy'},
}

remainder_element = 'Fe'

Cmp_map = {
    'CNC': 4.2571,
}

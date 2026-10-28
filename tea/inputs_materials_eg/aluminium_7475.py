"""
Aluminium 7475

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
1.4900. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Aluminium 7475'

density = 2.7

composition = {
    'Cu': {'wt': 1.5499999999999998, 'type': 'alloy'},
    'Fe': {'wt': 0.12, 'type': 'alloy'},
    'Mg': {'wt': 2.25, 'type': 'alloy'},
    'Mn': {'wt': 0.06, 'type': 'alloy'},
    'Si': {'wt': 0.1, 'type': 'alloy'},
    'Ti': {'wt': 0.06, 'type': 'alloy'},
    'Zn': {'wt': 5.7, 'type': 'alloy'},
}

remainder_element = 'Al'

Cmp_map = {
    'CNC': 1.4900,
}

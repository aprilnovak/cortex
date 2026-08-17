"""
Material C93200 (G-CuSn75Pb)

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
2.2406. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Material C93200 (G-CuSn75Pb)'

density = 8.6

composition = {
    'Al': {'wt': 0.0025, 'type': 'alloy'},
    'Fe': {'wt': 0.1, 'type': 'alloy'},
    'Ni': {'wt': 0.5, 'type': 'alloy'},
    'P': {'wt': 0.075, 'type': 'alloy'},
    'Pb': {'wt': 7.0, 'type': 'alloy'},
    'S': {'wt': 0.04, 'type': 'alloy'},
    'Sb': {'wt': 0.175, 'type': 'alloy'},
    'Si': {'wt': 0.0025, 'type': 'alloy'},
    'Sn': {'wt': 6.9, 'type': 'alloy'},
    'Zn': {'wt': 2.5, 'type': 'alloy'},
}

remainder_element = 'Cu'

Cmp_map = {
    'CNC': 2.2406,
}

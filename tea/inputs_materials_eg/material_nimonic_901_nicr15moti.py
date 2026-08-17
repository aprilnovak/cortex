"""
Material Nimonic 901 (NiCr15MoTi)

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
14.1905. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Material Nimonic 901 (NiCr15MoTi)'

density = 8.3

composition = {
    'C': {'wt': 0.05, 'type': 'alloy'},
    'Cr': {'wt': 12.5, 'type': 'alloy'},
    'Fe': {'wt': 35.6, 'type': 'alloy'},
    'Mo': {'wt': 6.0, 'type': 'alloy'},
    'Ti': {'wt': 2.8, 'type': 'alloy'},
}

remainder_element = 'Ni'

Cmp_map = {
    'CNC': 14.1905,
}

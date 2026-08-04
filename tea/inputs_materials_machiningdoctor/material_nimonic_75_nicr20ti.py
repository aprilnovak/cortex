"""
Material Nimonic 75 (NiCr20Ti)

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
10.6429. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Material Nimonic 75 (NiCr20Ti)'

density = 8.3

composition = {
    'Al': {'wt': 0.1, 'type': 'alloy'},
    'C': {'wt': 0.45, 'type': 'alloy'},
    'Cr': {'wt': 20.0, 'type': 'alloy'},
    'Fe': {'wt': 4.0, 'type': 'alloy'},
    'Mn': {'wt': 0.45, 'type': 'alloy'},
    'Si': {'wt': 0.45, 'type': 'alloy'},
    'Ti': {'wt': 0.35, 'type': 'alloy'},
}

remainder_element = 'Ni'

Cmp_map = {
    'CNC': 10.6429,
}

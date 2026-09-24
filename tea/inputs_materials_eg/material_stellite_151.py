"""
Material Stellite 151

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
29.8000. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Material Stellite 151'

density = 8.3

composition = {
    'C': {'wt': 0.5, 'type': 'alloy'},
    'Cr': {'wt': 20.0, 'type': 'alloy'},
    'Fe': {'wt': 2.0, 'type': 'alloy'},
    'Ni': {'wt': 1.0, 'type': 'alloy'},
    'W': {'wt': 13.0, 'type': 'alloy'},
}

remainder_element = 'Co'

Cmp_map = {
    'CNC': 29.8000,
}

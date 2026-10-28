"""
Material RA 330

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
4.9667. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Material RA 330'

density = 8.3

composition = {
    'C': {'wt': 0.04, 'type': 'alloy'},
    'Cr': {'wt': 18.5, 'type': 'alloy'},
    'Fe': {'wt': 43.8, 'type': 'alloy'},
    'Mn': {'wt': 1.0, 'type': 'alloy'},
    'P': {'wt': 0.01, 'type': 'alloy'},
    'S': {'wt': 0.01, 'type': 'alloy'},
    'Si': {'wt': 1.13, 'type': 'alloy'},
}

remainder_element = 'Ni'

Cmp_map = {
    'CNC': 4.9667,
}

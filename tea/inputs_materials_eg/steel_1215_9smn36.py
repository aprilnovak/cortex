"""
Steel 1215 (9SMn36)

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
1.2521. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Steel 1215 (9SMn36)'

density = 7.85

composition = {
    'C': {'wt': 0.09, 'type': 'alloy'},
    'Mn': {'wt': 0.9, 'type': 'alloy'},
    'P': {'wt': 0.065, 'type': 'alloy'},
    'S': {'wt': 0.305, 'type': 'alloy'},
}

remainder_element = 'Fe'

Cmp_map = {
    'CNC': 1.2521,
}

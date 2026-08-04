"""
Steel 4340 (35CrNiMo6)

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
2.9505. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Steel 4340 (35CrNiMo6)'

density = 7.85

composition = {
    'C': {'wt': 0.405, 'type': 'alloy'},
    'Cr': {'wt': 0.8, 'type': 'alloy'},
    'Mn': {'wt': 0.7, 'type': 'alloy'},
    'Mo': {'wt': 0.25, 'type': 'alloy'},
    'Ni': {'wt': 1.825, 'type': 'alloy'},
}

remainder_element = 'Fe'

Cmp_map = {
    'CNC': 2.9505,
}

"""
Cast Iron GGG-80 (GGG-80)

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
5.1379. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Cast Iron GGG-80 (GGG-80)'

density = 7.2

composition = {
    'C': {'wt': 3.85, 'type': 'alloy'},
    'Cu': {'wt': 0.82, 'type': 'alloy'},
    'Mn': {'wt': 0.5, 'type': 'alloy'},
    'Mo': {'wt': 0.039, 'type': 'alloy'},
    'P': {'wt': 0.035, 'type': 'alloy'},
    'S': {'wt': 0.015, 'type': 'alloy'},
    'Si': {'wt': 2.5, 'type': 'alloy'},
}

remainder_element = 'Fe'

Cmp_map = {
    'CNC': 5.1379,
}

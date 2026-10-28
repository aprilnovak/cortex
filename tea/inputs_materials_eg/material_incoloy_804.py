"""
Material Incoloy 804

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
6.2083. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Material Incoloy 804'

density = 8.3

composition = {
    'Al': {'wt': 0.25, 'type': 'alloy'},
    'C': {'wt': 0.06, 'type': 'alloy'},
    'Cr': {'wt': 29.5, 'type': 'alloy'},
    'Fe': {'wt': 25.4, 'type': 'alloy'},
    'Mn': {'wt': 0.75, 'type': 'alloy'},
    'Si': {'wt': 0.5, 'type': 'alloy'},
    'Ti': {'wt': 0.6, 'type': 'alloy'},
}

remainder_element = 'Ni'

Cmp_map = {
    'CNC': 6.2083,
}

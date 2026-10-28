"""
Material Inconel 800H

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

name = 'Material Inconel 800H'

density = 8.3

composition = {
    'Al': {'wt': 0.4, 'type': 'alloy'},
    'C': {'wt': 0.08, 'type': 'alloy'},
    'Cr': {'wt': 21.0, 'type': 'alloy'},
    'Cu': {'wt': 0.3, 'type': 'alloy'},
    'Fe': {'wt': 44.0, 'type': 'alloy'},
    'Mn': {'wt': 1.0, 'type': 'alloy'},
    'P': {'wt': 0.02, 'type': 'alloy'},
    'S': {'wt': 0.01, 'type': 'alloy'},
    'Si': {'wt': 0.35, 'type': 'alloy'},
    'Ti': {'wt': 0.4, 'type': 'alloy'},
}

remainder_element = 'Ni'

Cmp_map = {
    'CNC': 4.9667,
}

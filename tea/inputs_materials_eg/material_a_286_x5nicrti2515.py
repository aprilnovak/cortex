"""
Material A-286 (X5NiCrTi2515)

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
4.6562. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Material A-286 (X5NiCrTi2515)'

density = 8.3

composition = {
    'Al': {'wt': 0.2, 'type': 'alloy'},
    'C': {'wt': 0.05, 'type': 'alloy'},
    'Cr': {'wt': 14.5, 'type': 'alloy'},
    'Mn': {'wt': 0.2, 'type': 'alloy'},
    'Mo': {'wt': 1.3, 'type': 'alloy'},
    'Ni': {'wt': 25.0, 'type': 'alloy'},
    'Si': {'wt': 0.2, 'type': 'alloy'},
    'Ti': {'wt': 2.1, 'type': 'alloy'},
}

remainder_element = 'Fe'

Cmp_map = {
    'CNC': 4.6562,
}

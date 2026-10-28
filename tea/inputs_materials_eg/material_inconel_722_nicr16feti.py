"""
Material Inconel 722 (NiCr16FeTi)

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
12.9565. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Material Inconel 722 (NiCr16FeTi)'

density = 8.3

composition = {
    'Al': {'wt': 0.7, 'type': 'alloy'},
    'C': {'wt': 0.04, 'type': 'alloy'},
    'Cr': {'wt': 16.0, 'type': 'alloy'},
    'Fe': {'wt': 7.0, 'type': 'alloy'},
    'Ti': {'wt': 2.4, 'type': 'alloy'},
}

remainder_element = 'Ni'

Cmp_map = {
    'CNC': 12.9565,
}

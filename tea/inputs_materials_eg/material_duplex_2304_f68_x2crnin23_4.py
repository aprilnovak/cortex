"""
Material Duplex 2304 (F68) (X2CrNiN23-4)

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
8.5143. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Material Duplex 2304 (F68) (X2CrNiN23-4)'

density = 7.9

composition = {
    'C': {'wt': 0.03, 'type': 'alloy'},
    'Cr': {'wt': 22.5, 'type': 'alloy'},
    'Mn': {'wt': 2.0, 'type': 'alloy'},
    'Ni': {'wt': 4.5, 'type': 'alloy'},
    'P': {'wt': 0.04, 'type': 'alloy'},
    'S': {'wt': 0.02, 'type': 'alloy'},
    'Si': {'wt': 1.0, 'type': 'alloy'},
}

remainder_element = 'Fe'

Cmp_map = {
    'CNC': 8.5143,
}

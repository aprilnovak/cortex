"""
Aluminium 2014 (AlCuSiMn)

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
1.6742. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Aluminium 2014 (AlCuSiMn)'

density = 2.7

composition = {
    'Cr': {'wt': 0.1, 'type': 'alloy'},
    'Cu': {'wt': 4.45, 'type': 'alloy'},
    'Fe': {'wt': 0.7, 'type': 'alloy'},
    'Mg': {'wt': 0.5, 'type': 'alloy'},
    'Mn': {'wt': 0.8, 'type': 'alloy'},
    'Si': {'wt': 0.85, 'type': 'alloy'},
    'Ti': {'wt': 0.15, 'type': 'alloy'},
    'Zn': {'wt': 0.25, 'type': 'alloy'},
}

remainder_element = 'Al'

Cmp_map = {
    'CNC': 1.6742,
}

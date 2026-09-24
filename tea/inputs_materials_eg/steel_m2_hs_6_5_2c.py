"""
Steel M2 (HS-6-5-2C)

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
4.3824. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Steel M2 (HS-6-5-2C)'

density = 7.85

composition = {
    'C': {'wt': 1.0, 'type': 'alloy'},
    'Cr': {'wt': 4.125, 'type': 'alloy'},
    'Mn': {'wt': 0.275, 'type': 'alloy'},
    'Mo': {'wt': 5.0, 'type': 'alloy'},
    'Ni': {'wt': 0.15, 'type': 'alloy'},
    'Si': {'wt': 0.325, 'type': 'alloy'},
    'V': {'wt': 1.975, 'type': 'alloy'},
    'W': {'wt': 6.125, 'type': 'alloy'},
}

remainder_element = 'Fe'

Cmp_map = {
    'CNC': 4.3824,
}

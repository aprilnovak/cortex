"""
Material Ti-99.7 (Grade 3) (Ti-99.7)

Auto-generated from the machiningdoctor.com scrape by
generate_material_inputs.py, for exercising the TEA cost-model toolkit
against a large material set. NOT a validated input file - see
inputs_materials/example_material_input_file.py for the full field spec.

density: coarse per-family default (g/cm3), not validated for this alloy.
Cmp_map['CNC']: max(turning_vc_m_min in db) / turning_vc_m_min(material) =
5.3214. A machinability proxy scaled so the most machinable material
in the database is 1.0 and harder-to-machine materials are >1, matching
the convention used in other Cmp_map files (>1 means costlier than ideal);
treat as a placeholder for testing, not a validated cost coefficient.
"""

name = 'Material Ti-99.7 (Grade 3) (Ti-99.7)'

density = 4.5

composition = {
    'C': {'wt': 0.1, 'type': 'alloy'},
    'Fe': {'wt': 0.3, 'type': 'alloy'},
}

remainder_element = 'Ti'

Cmp_map = {
    'CNC': 5.3214,
}

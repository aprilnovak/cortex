# ===========================================
# Example usage
# Sould be run from within the tea directory
# ===========================================

import tea_level1

import sys
import os
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'materials'))
sys.path.append(module_path)
print(module_path)
import materials

#--------------------------------
# Example usage of TEA Level 1
#--------------------------------   
eurofer97 = materials.eurofer97(8.0)  # density in g/cc
eurofer97_cost = tea_level1.TEA_Lv1(eurofer97)
print(f'Eurofer97 Cost: ${eurofer97_cost:.2f}/kg')
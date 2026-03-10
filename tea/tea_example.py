# ===========================================
# Example usage
# Sould be run from within the tea directory
# ===========================================

print("Importing tea module... ")
import tea

#--------------------------------
# Example usage of TEA Level 1
#--------------------------------
componet1_material = tea.ManufacturingMaterials_Database["Eurofer97"]   
component1_Lv1_cost = componet1_material.price_mass
print(f'Level 1 Cost: ${component1_Lv1_cost:.2f}/kg')

#--------------------------------
# Example usage of TEA Level 2
#--------------------------------   
# TODO: Convert cost to $/kg
component1_process1 = tea.ManufacturingProcesses_Database["HCEM"]
component1_process2 = tea.ManufacturingProcesses_Database["AM"] 
component1 = tea.ManufacturingComponent(name = "Component1", 
                                        volume_mm3 = 100000,
                                        material = componet1_material,
                                        processes = [component1_process1,component1_process2], 
                                        production_qty = 1000,
                                        cc=1.0, 
                                        cs=1.0, 
                                        ct=1.0, 
                                        cf=1.0, 
                                        waste_coeff=1.0
                                        )
component1_Lv2_cost = component1.manufacturing_cost()
print(f'Level 2 Cost: ${component1_Lv2_cost:.2f}/$/kg')
# ===========================================
# Example usage
# Sould be run from within the tea directory
# ===========================================

print("Importing tea module... ")
import tea

#--------------------------------
# Example usage of TEA Level 1
#--------------------------------
component1_material = tea.ManufacturingMaterials_Database["Eurofer97"]   
component1_Lv1_cost = component1_material.price_mass
print(f'Level 1 Cost: ${component1_Lv1_cost:.2f}/kg')

#--------------------------------
# Example usage of TEA Level 2
#--------------------------------   
# TODO: Convert cost to $/kg
component1_process1 = tea.ManufacturingProcesses_Database["HCEM"]
component1_process2 = tea.ManufacturingProcesses_Database["AM"] 
component1 = tea.ManufacturingComponent(name = "Component1", 
                                        volume_mm3 = 100000,
                                        material = component1_material,
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
#--------------------------------
# Example usage of TEA Level 2 (First wall stuctural material))
#--------------------------------   
component3_material = tea.ManufacturingMaterials_Database["Eurofer97"]     
component3_process1 = tea.ManufacturingProcesses_Database["Electron Beam"]
component3_process2 = tea.ManufacturingProcesses_Database["Cold Rolling"]
component3_process3 = tea.ManufacturingProcesses_Database["HIP"]
component3 = tea.ManufacturingComponent(name = "First Wall", 
                                        volume_mm3 = 3000000, #TODO: update to actual volume of first wall
                                        material = component3_material,
                                        processes = [component3_process1, component3_process2, component3_process3], 
                                        production_qty = 1000,
                                        cc=1.0, 
                                        cs=1.0, 
                                        ct=1.0, 
                                        cf=1.0, 
                                        waste_coeff=1.0
                                        )
component3_cost = component3.manufacturing_cost()
print(f'Component 3 Cost: ${component3_cost:.2f}/kg')

#--------------------------------
# Example usage of TEA Level 2 (Tungsten Spray Deposition vs ODS Hot Rolling + Cold Rolling + HIP)
#--------------------------------   
component2a_material1 = tea.ManufacturingMaterials_Database["Tungsten"]   
component2a_material2 = tea.ManufacturingMaterials_Database["Eurofer97"]   
component2a_process1 = tea.ManufacturingProcesses_Database["Spray Deposition"]
component2a = tea.ManufacturingComponent(name = "Tungsten", 
                                        volume_mm3 = 3000000*.75,  # 75% tungsten, 25% eurofer97    
                                        material = component2a_material1,
                                        processes = [component2a_process1], 
                                        production_qty = 1000,
                                        cc=1.0, 
                                        cs=1.0, 
                                        ct=1.0, 
                                        cf=1.0, 
                                        waste_coeff=1.0
                                        )
component2a_material2_process1_cost = component2a.manufacturing_cost()
component2a_material1_cost = component2a_material1.price_mass*(3000000*.25/1e9*8000)
component2a_cost = component2a_material1_cost + component2a_material2_process1_cost
print(f'Component 2a Cost: ${component2a_cost:.2f}/kg')

component2b_material = tea.ManufacturingMaterials_Database["ODS"]     
component2b_process1 = tea.ManufacturingProcesses_Database["Hot Rolling"]
component2b_process2 = tea.ManufacturingProcesses_Database["Cold Rolling"]
component2b_process3 = tea.ManufacturingProcesses_Database["HIP"]
component2b = tea.ManufacturingComponent(name = "ODS", 
                                        volume_mm3 = 3000000,
                                        material = component2b_material,
                                        processes = [component2b_process1, component2b_process2, component2b_process3], 
                                        production_qty = 1000,
                                        cc=1.0, 
                                        cs=1.0, 
                                        ct=1.0, 
                                        cf=1.0, 
                                        waste_coeff=1.0
                                        )
component2b_cost = component2b.manufacturing_cost()
print(f'Component 2b Cost: ${component2b_cost:.2f}/kg')

summary_a = component2a.summary()
summary_b = component2b.summary()

import pandas as pd

df = pd.DataFrame([
    {"Component": summary_a["Component"], **summary_a["Cost Breakdown"]},
    {"Component": summary_b["Component"], **summary_b["Cost Breakdown"]}
]).fillna(0)


import matplotlib.pyplot as plt

# --- CSU color palette (your version) ---
csu_colors = {
    "Material": "#1E4D2B",
    "Spray Deposition": "#E56954",
    "Hot Rolling": "#E77103",
    "Cold Rolling": "#7E5373",
    "HIP": "#008FB2"
}

# --- Prepare data ---
df_plot = df.set_index("Component")

# Ensure Material is first (bottom of stack)
cols = ["Material"] + [c for c in df_plot.columns if c != "Material"]
df_plot = df_plot[cols]

# Apply colors in correct order
colors = [csu_colors.get(col, "#999999") for col in df_plot.columns]

# --- Plot ---
plt.figure(figsize=(5, 6))

ax = df_plot.plot(
    kind="bar",
    stacked=True,
    color=colors,
    edgecolor="none"
)

# --- Labels (UPDATED TO $/component) ---
plt.ylabel("Cost ($/component)", fontsize=18, color="#59595B")
plt.xlabel("")

# --- Legend styling ---
legend = plt.legend(
    fontsize=14,
    frameon=False,
    loc="upper right"
)

for text in legend.get_texts():
    text.set_color("#59595B")

# --- Axis styling ---
for spine in ax.spines.values():
    spine.set_color("#59595B")

ax.tick_params(axis='both', colors="#59595B", labelsize=18)

plt.xticks(rotation=0)
plt.yticks(color="#59595B")

# --- Add total labels on top ---
totals = df_plot.sum(axis=1)
for i, total in enumerate(totals):
    ax.text(
        i, total,
        f"{total:.1f}",
        ha='center',
        va='bottom',
        fontsize=14,
        color="#59595B"
    )

plt.tight_layout()
plt.show()

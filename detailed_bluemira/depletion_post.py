# sdr_bluemira_postprocess.py
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import openmc
import openmc.deplete
from neutronics_model import cell_atomic_ratio, xcentroids_ob

#from full_setup import all_cells, xcentroids_ob, cell_ids  # geometry info

# Some reference papers:
# https://doi.org/10.1016/j.fusengdes.2021.112428
# https://doi.org/10.1016/j.fusengdes.2021.112338
# https://iopscience.iop.org/article/10.1088/1741-4326/aca61f

# ------------------------------------------------------------------
# Load depletion results
# ------------------------------------------------------------------
results_path = Path("r2s/activation/depletion_results.h5").resolve()
if not results_path.exists():
    raise FileNotFoundError(f"Results file not found: {results_path}")

chain_file = Path("bluemira_chain.xml").resolve()

# MUST CHANGE: set chain file BEFORE constructing Results (more consistent)
openmc.config['chain_file'] = str(chain_file)

results = openmc.deplete.Results(str(results_path))

source_rates = results.get_source_rates()

# Time conversions
s_in_d = 86400.0
d_in_y = 365.25
s_in_y = d_in_y * s_in_d

# Region names:
cell_ids = [49, 45, 79, 66, 82, 94, 103, 112, 124, 132, 2]
cell_id_to_name = {
    49:  "Armor",
    45:  "First Wall",
    79:  "OB_1",
    66:  "OB_2",
    82:  "OB_3",
    94:  "OB_4",
    103: "OB_5",
    112: "OB_6",
    124: "OB_7",
    132: "OB_8",
    2:   "VV",
}
 
# load cell_id - material mapping
map_df = pd.read_csv("r2s/activation/cell_material_map.csv")

# MUST CHANGE: normalize mat_id once
map_df["mat_id"] = map_df["mat_id"].astype(str)

# cell_id -> mat_id (1-to-1, safe)
cell_to_mat = dict(zip(map_df["cell_id"], map_df["mat_id"]))

# MUST CHANGE: mat_id -> list[cell_id] (SAFE replacement for mat_to_cell)
cells_by_mat = (
    map_df.groupby("mat_id")["cell_id"]
          .apply(list)
          .to_dict()
)

# mat_id -> material_name (use groupby-first to avoid silent overwrite)
mat_id_to_name = (
    map_df.groupby("mat_id")["material_name"]
          .first()
          .to_dict()
)

# Get volumes per material from first step
step0 = results[0]
vol_by_mat = step0.volume 

# Create depletion folder for output plots and results
out_dir = Path("./depletion_results")
out_dir.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------
# Activity
# ------------------------------------------------------------------

# Total Activity 
total_activity_all = {}   # cid → (t_rel_years, total_act_plot)

# -------------------------------------
# Total (and per nuclide) activity vs time
# For each cell:
#   1 - Find the 'n' top activity nuclides in each time step
#   2 - Save all top activity nuclides of each time step
#   3 - Plot the activity of all nuclides saved for all cooling time steps
# -------------------------------------
for cid in cell_ids:
    # get deplet_material_id and volume
    mat_id = str(cell_to_mat[cid])
    vol    = vol_by_mat[mat_id]

    # total activity
    _, total_act = results.get_activity(
        mat=mat_id,
        by_nuclide=False,
        units='Bq/kg',
        volume=vol
    )
    
    # Relevant top nuclides (arbitrary number)
    top_n = 5

    # Per-nuclide activity for this cell
    time_nuc, act_dict_list = results.get_activity(
        mat=mat_id,
        by_nuclide=True,
        units='Bq/kg',
        volume=vol
    )

    # Build top_n sets for each timestep (t > 0) 
    selected_nucs = set()

    # Obtain the top_n nuclides in decreasing order (activity)
    for istep in range(1, len(time_nuc)):   # skip t=0 (initial)
        act_dict = act_dict_list[istep]
        if not act_dict:
            continue

        step_sorted = sorted(
            act_dict.items(),
            key=lambda x: x[1],
            reverse=True
        )
        step_topN = [nuc for nuc, val in step_sorted[:top_n]]
        selected_nucs.update(step_topN)

    selected_nucs = sorted(selected_nucs)

    # Find end-of-irradiation index (last step with non-zero source)
    irr_intervals = np.where(source_rates > 0.0)[0]

    if irr_intervals.size > 0:
        idx_shutdown = irr_intervals[-1] + 1
        if idx_shutdown >= len(time_nuc):
            idx_shutdown = len(time_nuc) - 1
    else:
        # fallback: no irradiation, treat step 0 as shutdown
        idx_shutdown = 0

    # Adjust time vector to remove initial irradiation time (plotting purposes)
    t_rel = time_nuc - time_nuc[idx_shutdown]   # t=0 at shutdown | time_nuc[1] end of irradiation
    mask = t_rel > 0.0                          # boolean filter to exclude t<=0
    t_rel_plot = t_rel[mask] / s_in_y           # convert into years
    total_act_plot = total_act[mask]

    # store for combined plot later
    total_activity_all[cid] = (t_rel_plot, total_act_plot)

    # Plotting
    # create unique colors for each nuclide (too many nuclides, needed more colors)
    num_nucs = len(selected_nucs)
    colors = plt.cm.tab20(np.linspace(0, 1, num_nucs))

    plt.figure()

    # plot each individual nuclide activity 
    for i, nuc in enumerate(selected_nucs):
        vals = np.array([d.get(nuc, 0.0) for d in act_dict_list])
        vals_plot = vals[mask]

        if np.all(vals_plot == 0.0):
            continue

        plt.loglog(t_rel_plot, vals_plot, marker='o', label=nuc, color=colors[i])

    # Plot total activity on same axes
    plt.loglog(t_rel_plot, total_act_plot, color='black', linewidth=2.5,
               label="Total Activity")

    plt.xlabel("Time after irradiation [years]")
    plt.ylabel("Activity [Bq/kg]")
    region = cell_id_to_name[cid]  # get correct cell for that depleted_material
    plt.title(f"Top nuclides — {region} (mat {mat_id})")
    plt.grid(True, which='both')

    ymax_local = np.max(total_act_plot) # max for addapting plot y_lim
    plt.ylim(1e0, 10*ymax_local) 
    #plt.xlim(1e-10,2e3)

    # vertical lines
    plt.axvline(1.0/365, color='gray', linestyle='--', linewidth=1)
    plt.axvline(1,       color='gray', linestyle='--', linewidth=1)
    plt.axvline(100,     color='gray', linestyle='--', linewidth=1)

    plt.legend(
        fontsize=8,
        title="Nuclides",
        loc='center left',
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0.
    )
    plt.subplots_adjust(right=0.75)
    region = cell_id_to_name[cid]  
    plt.savefig(f"./depletion_results/activity_nuclides_cell_{region}.png", dpi=200)
    plt.close()

# -------------------------------------
# Combined total activity plot for all cells
# -------------------------------------
plt.figure(figsize=(10, 6))

num_cells = len(total_activity_all)
colors = plt.cm.tab20(np.linspace(0, 1, num_cells))

for i, (cid, (t_rel_y, total_act_plot)) in enumerate(total_activity_all.items()):
    region = cell_id_to_name[cid] 
    plt.loglog(t_rel_y, total_act_plot, marker='o',
               label=region, color=colors[i])

plt.xlabel("Time after irradiation [years]")
plt.ylabel("Activity [Bq/kg]")
plt.title("Total Activity over cooling time")
plt.grid(True, which='both')

plt.axvline(1.0/365, color='gray', linestyle='--', linewidth=1)  # 1 day
plt.axvline(1,       color='gray', linestyle='--', linewidth=1)  # 1 year
plt.axvline(100,     color='gray', linestyle='--', linewidth=1)  # 100 years

ymin, ymax = plt.ylim()
y_top = ymax * 0.85  
plt.xlim(1e-10,2e3)

plt.text(1.0/365, y_top, "1 day",   rotation=90,
         va='top', ha='center', color='gray')
plt.text(1,       y_top, "1 year",  rotation=90,
         va='top', ha='center', color='gray')
plt.text(100,     y_top, "100 year", rotation=90,
         va='top', ha='center', color='gray')

plt.legend(title="Regions", fontsize=8,
           loc='center left', bbox_to_anchor=(1.02, 0.5))
plt.subplots_adjust(right=0.75)

plt.savefig("./depletion_results/activity_all_cells.png", dpi=200)
plt.close()

# --------------------------
# Total activity over radial depth (for different time_steps)
# For each time_step:
#       1 - Plot total activity for each cell_id
# --------------------------

# Build a matrix: radial_activity[i_cell, j_time] ---
# xcentroids_ob = (0.1, 1.1, 6.0, 15.0, 25.0, 35.0, 45.0, 55.0, 70.0, 90.0, 157.2) -> imported from bluemira_model
x = np.array(xcentroids_ob) # radial positions
n_cells = len(cell_ids)

# use the time axis from the first cell
first_cid = cell_ids[0]
t_rel_y, _ = total_activity_all[first_cid]   # years after shutdown
n_times = len(t_rel_y)

radial_activity = np.zeros((n_cells, n_times))

for i, cid in enumerate(cell_ids):
    t_y, act = total_activity_all[cid]
    radial_activity[i, :] = act

plt.figure(figsize=(8, 6))

# Choose which time indices to plot
idx_to_plot = [0, 3, 6, 9, 12, 15]  # arbritary time_steps
idx_to_plot = [i for i in idx_to_plot if i < n_times] # filter out invalid indices

colors = plt.cm.tab20(np.linspace(0, 1, len(idx_to_plot)))

for k, j in enumerate(idx_to_plot):
    y = radial_activity[:, j]
    plt.semilogy(x, y, marker='o', color=colors[k],
                 label=f"t = {t_rel_y[j]:.2e} y")


plt.xlabel("Radial position [cm]")
plt.ylabel("Total activity [Bq/kg]")
plt.title("Radial profile of Total activity per cooling time")
plt.grid(True, which='both')

# legend outside
plt.legend(
    fontsize=8,
    title="Timesteps",
    loc='center left',
    bbox_to_anchor=(1.02, 0.5),
    borderaxespad=0.
)
plt.subplots_adjust(right=0.75)

plt.savefig("./depletion_results/radial_activity_profiles.png", dpi=200)
plt.close()


# ------------------------------------------------------------------
# Decay heat
# ------------------------------------------------------------------

total_heat_all = {}   # cid → (t_rel_years, total_heat_plot)
max_decay_comb = 0.0
top_n = 5 

# -------------------------------------
# Total (and per nuclide) Decay Heat vs time
# For each cell:
#   1 - Find the 'n' top decay heat nuclides in each time step
#   2 - Save all top decay heat nuclides of each time step
#   3 - Plot the decay heat of all nuclides saved for all cooling time steps
# -------------------------------------
for cid in cell_ids:
    mat_id = str(cell_to_mat[cid])
    vol    = vol_by_mat[mat_id]

    # total decay heat (no per-nuclide breakdown) 
    time_total_h, total_h = results.get_decay_heat(
        mat=mat_id,
        by_nuclide=False,
        units='W/cm3',
        volume=vol
    )

    # Per-nuclide decay heat (to get top nuclides)
    time_nuc_h, heat_dict_list = results.get_decay_heat(
        mat=mat_id,
        by_nuclide=True,
        units='W/cm3',
        volume=vol
    )

    # build set of nuclides that ever appear in top_n at any timestep (t>0)
    selected_nucs = set()

    # Obtain the top_n nuclides in decreasing order (decay heat)
    for istep in range(1, len(time_nuc_h)):  # skip t=0
        heat_dict = heat_dict_list[istep]
        if not heat_dict:
            continue

        step_sorted = sorted(
            heat_dict.items(),
            key=lambda x: x[1],
            reverse=True
        )
        step_topN = [nuc for nuc, val in step_sorted[:top_n]]
        selected_nucs.update(step_topN)

    selected_nucs = sorted(selected_nucs)
    
    # remove irradiation time
    irr_intervals = np.where(source_rates > 0.0)[0]
    if irr_intervals.size > 0:
        idx_shutdown = irr_intervals[-1] + 1
        if idx_shutdown >= len(source_rates) + 1:
            idx_shutdown = len(source_rates)
    else:
        # fallback: no irradiation
        idx_shutdown = 0

    # time axis: seconds after end of irradiation
    t_rel = time_total_h - time_total_h[idx_shutdown]    # t=0 at shutdown | time_total_h[1] end of irradiation
    mask = t_rel > 0.0                        # boolean filter to exclude t<=0
    t_rel_plot = t_rel[mask] / s_in_y        # convert into years
    total_h_plot = total_h[mask]

    # store for later combined / radial plots
    total_heat_all[cid] = (t_rel_plot, total_h_plot)

    # per-cell decay heat vs time (top nuclides + total) 
    num_nucs = len(selected_nucs)
    colors = plt.cm.tab20(np.linspace(0, 1, num_nucs))

    plt.figure()
    # plot per-nuclide decay heat
    for i, nuc in enumerate(selected_nucs):
        # heat for this nuclide at all timesteps
        vals = np.array([d.get(nuc, 0.0) for d in heat_dict_list])
        vals_plot = vals[mask]

        if np.all(vals_plot == 0.0):
            continue

        plt.loglog(t_rel_plot, vals_plot, marker='o',
                   label=nuc, color=colors[i])

    # plot total decay heat
    plt.loglog(t_rel_plot, total_h_plot, marker='o',
               color='black', linewidth=2.5, label='Total')

    plt.xlabel("Time after irradiation [years]")
    plt.ylabel("Decay heat [W/cm3]")
    region = cell_id_to_name[cid] 
    plt.title(f"Decay heat — {region} (mat {mat_id})")
    plt.grid(True, which='both')

    # Find maximum decay_heat (over all times and materials)
    # to change plot range
    max_decay = np.max(total_h_plot)
    if (max_decay >= max_decay_comb):
        max_decay_comb = max_decay
    plt.xlim(1e-10,2e3)
    plt.ylim(1e-16, 10 * max_decay)

    # legend outside
    plt.legend(
        fontsize=8,
        title="Nuclides",
        loc='center left',
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0.
    )
    plt.subplots_adjust(right=0.75)
    plt.savefig(f"./depletion_results/decayheat_nuclides_{region}.png", dpi=200)
    plt.close()

# -------------------------------------
# Total decay heat vs time for all cells 
# -------------------------------------
plt.figure()
num_cells = len(total_heat_all)
colors = plt.cm.tab20(np.linspace(0, 1, num_cells))

for i, (cid, (t_rel_y, total_h_plot)) in enumerate(total_heat_all.items()):
    plt.loglog(t_rel_y, total_h_plot, marker='o',
               label=f"Cell {cid}", color=colors[i])

plt.xlabel("Time after irradiation [years]")
plt.ylabel("Decay heat [W/cm3]")
plt.title("Total decay heat")
plt.grid(True, which='both')

plt.axvline(1.0/365, color='gray', linestyle='--', linewidth=1)
plt.axvline(1,       color='gray', linestyle='--', linewidth=1)
plt.axvline(100,     color='gray', linestyle='--', linewidth=1)

plt.xlim(1e-10,2e3)
plt.ylim(1e-12, 10*max_decay_comb)

plt.legend(title="Cells", fontsize=8,
           loc='center left', bbox_to_anchor=(1.02, 0.5))
plt.subplots_adjust(right=0.75)

plt.savefig("./depletion_results/decayheat_all_cells.png", dpi=200)
plt.close()

# -------------------------------------
# Radial profile of decay heat at different times
# -------------------------------------

x = np.array(xcentroids_ob)    # radial positions (cm), aligned with cell_ids
n_cells = len(cell_ids)

# time axis from the first cell
first_cid = cell_ids[0]
t_rel_y, _ = total_heat_all[first_cid]   
n_times = len(t_rel_y)

radial_heat = np.zeros((n_cells, n_times))

for i, cid in enumerate(cell_ids):
    t_y, h = total_heat_all[cid]
    if not np.allclose(t_y, t_rel_y):
        print(f"WARNING: time grid mismatch for cell {cid}")
    radial_heat[i, :] = h

plt.figure(figsize=(8, 6))

# which time indices to plot
idx_to_plot = [0, 3, 6, 9, 12, 15] 
idx_to_plot = [i for i in idx_to_plot if i < n_times]

colors = plt.cm.tab20(np.linspace(0, 1, len(idx_to_plot)))

for k, j in enumerate(idx_to_plot):
    y = radial_heat[:, j]
    plt.semilogy(x, y, marker='o', color=colors[k],
                 label=f"t = {t_rel_y[j]:.2e} y")

plt.xlabel("Radial position [cm]")
plt.ylabel("Decay heat [W/cm3]")
plt.title("Radial profile of decay heat per cooling time")
plt.grid(True, which='both')

plt.legend(
    fontsize=8,
    title="Timesteps",
    loc='center left',
    bbox_to_anchor=(1.02, 0.5),
    borderaxespad=0.
)
plt.subplots_adjust(right=0.75)

plt.savefig("./depletion_results/radial_decayheat_profiles.png", dpi=200)
plt.close()

# ------------------------------------------------------------------
# Material transmutation
# ------------------------------------------------------------------
# For each activation material mapped from cell_ids, we output:
#   1 - original neutronics composition -> (from materials.xml, no depletion)
#   2 - composition after irradiation   -> (end of first non-zero source step)
#   3 - final composition               -> (last depletion step, after cooling)

# COMMENTS:
# results.get_atoms get all atoms or just the requested nuclide?
# (mat: Material | str, nuc: str, nuc_units: str = 'atoms', time_units: str = 's')
# → tuple[ndarray, ndarray]
# maybe change logic for the material transmutation/storage/plotting
# TBD: discuss with Patrick and April what would be the best way to display this information

def atom_fractions_from_dep_step(results_obj, step_index, mat_id):
    """
    Get atom fractions for material `mat_id` at depletion step `step_index`.
    Handles both dict and pandas Series outputs from get_nuclide_atom_densities().
    """
    step = results_obj[step_index]
    dep_mat = step.get_material(mat_id)

    data = dep_mat.get_nuclide_atom_densities()

    # Case 1: dict {nuc: density}
    if isinstance(data, dict):
        total = sum(data.values())
        if total == 0.0:
            return {}
        return {nuc: val / total for nuc, val in data.items()}

    # Case 2: pandas Series
    total = data.sum()
    if total == 0.0:
        return {}
    return (data / total).to_dict()


# Determine depletion step indices
times_y = results.get_times(time_units='a')     # times for each stored state
source_rates = results.get_source_rates()       # source rates per interval
n_steps = len(times_y)

# initial (pre-depletion / t=0)
idx_initial = 0

# find last interval with non-zero source;
# composition "after irradiation" is at the end of that interval → +1
irr_intervals = np.where(source_rates > 0.0)[0]
if irr_intervals.size > 0:
    idx_after_irr = irr_intervals[-1] + 1
    # safety cap
    if idx_after_irr >= n_steps:
        idx_after_irr = n_steps - 1
else:
    # fallback: if all source_rates are zero, treat step 0 as "after irr"
    idx_after_irr = 0

# final step: last stored time
idx_final = n_steps - 1

print("Transmutation step indices:")
print(f"  initial (neutronics) ")
print(f"  after irradiation step = {idx_after_irr}, t = {times_y[idx_after_irr]:.6e} y")
print(f"  final step             = {idx_final}, t = {times_y[idx_final]:.6e} y")

# Materials of interest: those in cell_ids 
mat_ids_of_interest = sorted({str(cell_to_mat[cid]) for cid in cell_ids})

# Get for each material in the cell_ids find:
#   1 - Initial composition
#   2 - Composition after irradiation
#   3 - Composition after cooling
for mat_id in mat_ids_of_interest:
    # Get name (from the cell_material_map.csv)
    mat_name = mat_id_to_name.get(mat_id, f"mat_{mat_id}")

    # Original neutronics composition (no depletion)
    # Use cells_by_mat to map material -> cell(s)
    if mat_id not in cells_by_mat:
        print(f"[Transmutation] WARNING: mat_id {mat_id} not in cells_by_mat; skipping.")
        continue

    # representative cell for this material
    cids = sorted(cells_by_mat[mat_id])
    cid = cids[0]

    if cid not in cell_atomic_ratio:
        print(f"[Transmutation] WARNING: cell_atomic_ratio missing for cell {cid}; skipping mat_id {mat_id}.")
        continue

    af_initial_neutronics = cell_atomic_ratio[cid]

    # Composition after irradiation (from depletion)
    af_after_irr = atom_fractions_from_dep_step(results, idx_after_irr, mat_id)

    # Final composition (last cooling step)
    af_final = atom_fractions_from_dep_step(results, idx_final, mat_id)

    # Union of all nuclides seen in any of the three states
    all_nucs = sorted(set(af_initial_neutronics) | set(af_after_irr) | set(af_final))

    rows = []
    # Provide composition at different time_steps and
    # also differnce between then for each nuclide.
    for nuc in all_nucs:
        f0 = af_initial_neutronics.get(nuc, 0.0)
        f1 = af_after_irr.get(nuc, 0.0)
        fF = af_final.get(nuc, 0.0)

        rows.append({
            "nuclide": nuc,
            "atom_frac_initial": f0,
            "atom_frac_after_irradiation": f1,
            "atom_frac_final": fF,
            "delta_after_irr": (f1 - f0) * 1e6,
            "delta_final":     (fF - f0) * 1e6,
        })

    df = pd.DataFrame(rows)
    df = df.sort_values("atom_frac_initial", ascending=False)

    safe_name = mat_name.replace(" ", "_").replace("/", "_")
    out_csv = out_dir / f"transmutation_{safe_name}_mat{mat_id}_cells{len(cids)}.csv"
    df.to_csv(out_csv, index=False)

    print(f"[Transmutation] Wrote {out_csv}")
    # CREATE SECOND CSV SORTED BY FINAL FRACTION
    df_sorted_final = df.sort_values("atom_frac_final", ascending=False)

    out_csv_sorted = out_dir / f"transmutation_{safe_name}_mat{mat_id}_cells{len(cids)}_sorted_final.csv"
    df_sorted_final.to_csv(out_csv_sorted, index=False)

    print(f"[Transmutation] Wrote sorted table: {out_csv_sorted}")


# ADDITIONAL PLANNING 
# Add better visionlization for material composition/composition changes
# (which plots would be better) -> cleaner look into major materials and impurities.
# ADD INDIVIDUAL EVOLUTION OF H1, He3, He4

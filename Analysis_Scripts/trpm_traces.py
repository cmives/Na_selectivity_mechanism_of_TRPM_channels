#!/usr/bin/env python3

import sys
import random
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from itertools import islice, cycle
from plottify import autosize

sns.set_context("paper", font_scale=1.4)
sns.set_style("whitegrid")
sns.palplot(sns.color_palette("colorblind"))

#################################################################################
# SPECIFYING THE PROTEIN OF INTEREST AND THE PATH"S TO THE RELEVANT DIRECTORIES #
#################################################################################

protein = ["trpm5"]
interval_range = 50
start_frame = 0
last_frame = 250
V_tm = 0.35
all_or_permeating = "permeating"

################################################################################
#  #
################################################################################

for protein in protein:
    if protein == "trpm5":
        construct = "7MBS"
    ion_solution = ["kcl","nacl","cacl2"]

    # Reading in the position of the gates...
    df_gate = pd.read_csv(f"../binding_site_cog/binding_site_z_coordinates_{protein}_{construct}_nacl.csv")
    df_gate = df_gate.drop(df_gate[df_gate.BS == "BS_B"].index)
    gates_values = df_gate["Z_mean"].values.tolist()

    for ion_solution in ion_solution:
        if ion_solution == "cacl2":
            ionic_species = ["CAM"]
        if ion_solution == "charmm_cacl2":
            ionic_species = ["Cal"]
        elif ion_solution == "nacl":
            ionic_species = ["NA"]
        elif ion_solution == "kcl":
            ionic_species = ["K"]
        elif ion_solution == "mixed": 
            ionic_species = ["CAM","NA"]

        replicate = ["repl_a","repl_b","repl_c"]
        for replicate in replicate:
            df_data_list = []
            df_resid_list = []
            resid_colors_list = []
            for s in range(len(ionic_species)):
                # Importing the coordinate data...

                if all_or_permeating == "permeating":
                    try:
                        df_data = pd.read_csv(f"../permeation_identification/permeation_identification_pbc_{protein}_{construct}_{ion_solution}_{ionic_species[s]}_{replicate}.csv")
                        df_data = df_data.set_index("Timestep")
                        df_data = df_data.loc[:, df_data.columns.str.endswith("z")]
                        df_data.where(df_data > -45, np.nan, inplace=True)
                        df_data.where(df_data < 45, np.nan, inplace=True)
                        df_data = df_data.reset_index()
                        df_data_list.append(df_data)
                    except FileNotFoundError:
                        sys.exit(f"FILE NOT FOUND: {protein} {ion_solution} {replicate} {ionic_species[s]}")
                    
                    # Importing the permeating ion resids...
                    permeation_df = pd.read_csv(f"../permeation_identification/permeating_ion_resids_{protein}_{construct}_{ion_solution}_{ionic_species[s]}_{replicate}.csv")
                    permeating_ions = permeation_df.iloc[:,0].tolist()
                    df_resid_list.append(permeating_ions)

                    # Making a list of colors for plotting
                    #color_list = ["tab:blue","tab:orange","tab:green","tab:red","tab:purple","tab:brown","tab:pink","tab:cyan","tab:olive"]
                    color_list = ["gray"]
                    color_list_mod = list(islice(cycle(color_list), len(permeating_ions)))
                    resid_colors_list.append(color_list_mod)

                if all_or_permeating == "all":
                    try:
                        df_data = pd.read_csv(f"../all_cation_cartesians/cartesian_coordinates_{protein}_{construct}_{ion_solution}_{ionic_species[s]}_{replicate}.csv")
                        # Removing PBC effects...
                        df_data = df_data.set_index("Timestep")
                        df_data = df_data.loc[:, df_data.columns.str.endswith("z")]
                        df_data.where(df_data > -45, np.nan, inplace=True)
                        df_data.where(df_data < 45, np.nan, inplace=True)
                        df_data = df_data.reset_index()
                        df_data_list.append(df_data)
                    except FileNotFoundError:
                        sys.exit(f"FILE NOT FOUND: {protein} {ion_solution} {replicate} {ionic_species[s]}")

                    # Making a list of colors for plotting
                    #color_list = ["tab:blue","tab:green","tab:red","tab:purple","tab:brown","tab:cyan","tab:olive"]
                    color_list = ["tab:cyan"]
                    color_list_mod = list(islice(cycle(color_list), (len(df_data.columns)-1)))
                    resid_colors_list.append(color_list_mod)


            # Plotting the permeation trace...
            fig, ax = plt.subplots(figsize=(15,10))
            #plt.rcParams.update({'font.size':30})
            sns.set(font_scale = 2)

            for s in range(len(ionic_species)):
                if all_or_permeating == "permeating":
                    if ion_solution != "mixed":
                        for n in range(len(df_resid_list[s])):
                            ax.plot("Timestep", f"{df_resid_list[s][n]}_z", data=df_data_list[s], color=resid_colors_list[s][n], alpha=0.7)
                    elif (ion_solution == "mixed") and (ionic_species[s] == "CAM"):
                        for n in range(len(df_resid_list[s])):
                            ax.plot("Timestep", f"{df_resid_list[s][n]}_z", data=df_data_list[s], color="tab:orange", alpha=0.7)
                    elif (ion_solution == "mixed") and (ionic_species[s] == "NA"):
                        for n in range(len(df_resid_list[s])):
                            ax.plot("Timestep", f"{df_resid_list[s][n]}_z", data=df_data_list[s], color="tab:blue", alpha=0.7)

                if all_or_permeating == "all":
                    if ion_solution != "mixed":
                        for n in np.arange(1,len(df_data_list[s].columns)-1):
                            ax.plot("Timestep", f"{df_data_list[s].columns[n]}", data=df_data_list[s], color=resid_colors_list[s][n], alpha=0.7)
                    elif (ion_solution == "mixed") and (ionic_species[s] == "CAM"):
                        for n in np.arange(1,len(df_data_list[s].columns)-1):
                            ax.plot("Timestep", f"{df_data_list[s].columns[n]}", data=df_data_list[s], color="tab:orange", alpha=0.7)
                    elif (ion_solution == "mixed") and (ionic_species[s] == "NA"):
                        for n in np.arange(1,len(df_data_list[s].columns)-1):
                            ax.plot("Timestep", f"{df_data_list[s].columns[n]}", data=df_data_list[s], color="tab:blue", alpha=0.7)
                
            for gate_z, gate_colour in zip(gates_values,["tab:gray","tab:gray"]):
                ax.axhspan(gate_z-1, gate_z+1, alpha=0.35, facecolor=gate_colour)

            ax.grid()
            ax.set(ylim=(-30,25))
            ax.set(xlim=(0,250))
            ax.set_xlabel("Simulated time (ns)", fontsize=25)
            ax.set_ylabel("Distance along pore axis\n(Å)", fontsize=25)

            plt.title(f"Permeation traces of cations through {protein.upper()} in {ion_solution} ({replicate})\n",weight='bold')
            plt.xticks(fontsize=25)
            plt.yticks(fontsize=25)
            #autosize()
            ax.set_facecolor("w")

            if all_or_permeating == "permeating":
                fig.savefig(f"../permeation_plots/permeation_traces_{protein}_{construct}_{ion_solution}_{replicate}_permeating_ions_POSTER.png")
            if all_or_permeating == "all":
                fig.savefig(f"../permeation_plots/permeation_traces_{protein}_{construct}_{ion_solution}_{replicate}_all_ions_POSTER.png")

            plt.cla()


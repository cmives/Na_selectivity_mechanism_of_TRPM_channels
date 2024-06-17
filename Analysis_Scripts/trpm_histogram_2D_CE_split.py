#!/usr/bin/env python3
#$ -cwd
#$ -V
#$ -N m5_his_2c
#$ -e /cluster/uz_lab/cmives/TRP_channels/MDS/permeation_analysis/z_histogram/m5_his_ce.err
#$ -o /cluster/uz_lab/cmives/TRP_channels/MDS/permeation_analysis/z_histogram/m5_his_ce.out
#$ -pe smp 2
#$ -jc long

import os
import sys
import time
import numpy as np
import pandas as pd
import seaborn as sns
import MDAnalysis as mda
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#import matplotlib.ticker as ticker

#################################################################################
# SPECIFYING THE PROTEIN OF INTEREST AND THE PATH'S TO THE RELEVANT DIRECTORIES #
#################################################################################

protein = "trpm5"
forcefield = "charmm"
trpmds = "/cluster/uz_lab/cmives/TRP_channels/MDS/"
result_dir = "/cluster/uz_lab/cmives/TRP_channels/MDS/permeation_analysis/z_histogram"
voltage = "800mv_conc"
desired_time_step = 0.2     # in ns
bin_value = 1           # in A

#################################################################################
# Binning round #
#################################################################################

def bin_round(x, base=bin_value):
    return base * round(x/base)

################################################################################
# DEFINING FUNCTIONS #
###############################################################################
sns.set_context("paper", font_scale=1.4)
sns.set_style("whitegrid")

################################################################################
# STARTING THE ANALYSIS OF THE TRAJECTORIES #
################################################################################

if protein == "trpm5":
    construct = "7MBS"

ion_solution = ["mixed"]
for ion_solution in ion_solution:
    # Specifiying the names of the ionic species...
    if ion_solution == "cacl2":
        ionic_species = ["CAM"]

    elif ion_solution == "nacl":
        ionic_species = ["NA"]

    elif ion_solution == "kcl":
        ionic_species = ["K"]

    elif ion_solution == "mixed": 
        ionic_species = ["CAM"]#,"NA"]

    for ionic_species in ionic_species:

        df_list = []

        print("############################################################################################\nSCRIPT STARTED ON %s FOR %s %s %s\n############################################################################################\n" %(time.asctime().upper(),protein.upper(), construct.upper(), ion_solution.upper()))

        ############################################################
        #  DEFINING VARIABLES BASED UPON THE INPUT VARIABLES #
        ############################################################

        dir_stem = os.path.join(trpmds,protein,forcefield,construct,ion_solution,f"prodmd/compel/{voltage}/traj_splitting")

        traj_ts = int(desired_time_step/0.2)

        for system in ["ch0","ch1"]:

            for replicate in ["repl_a","repl_b","repl_c"]:

                ############################################################
                # IMPORTING THE COORDINATE AND TRAJECTORY FILES #
                ############################################################
                print("READING DATA FROM %s" %dir_stem)
                GRO = os.path.join(dir_stem,f"{system}_system.gro")
                XTC = GRO.replace(f"{system}_system.gro", f"{system}_{replicate}.xtc")
                print("STARTING TO LOAD TRAJECTORY")
                u = mda.Universe(GRO,XTC)
                print("TRAJECTORY LOADED\n")                        

                ############################################################
                # IDENTIFYING IONS THAT ENTER THE PORE #
                ############################################################
                #Specifying the pore region, so we can identify ions that are within the pore. This selection is set to update, so will update for every frame in the trajectory. I'm going to define the pore as the cylinder between the upper and lower gate.
                if protein == 'trpm5':
                    upper_gate_resid = 906  # Gln
                    lower_gate_resid = 966  # Ile


                ############################################################
                # READING IN THE PERMEATING ION RESIDS # 
                ###########################################################

                # Selecting cations within the pore...
                reference_pore = u.select_atoms(f"cyzone 10 15 -15 (backbone and (resid {upper_gate_resid} or resid {lower_gate_resid})) and not (resname TIP3 or resname POPC)", updating=True)
                pore = u.select_atoms(f"(cyzone 5 40 -40 (backbone and (resid {upper_gate_resid} or resid {lower_gate_resid}))) and resname {ionic_species}", updating=True) # Changed this to 5 A to catch the effect of the sponge sites.

                z_coordinate = []
                r_coordinate = []
                for frame in u.trajectory[::traj_ts]:
                    print(f"Analysing frame for {frame.time / 1000} ns")
                    for cation in pore.residues:
                        # Getting the z-coordinate of the cation...
                        if system == "ch0":
                            z_coordinate.append(cation.atoms.center_of_geometry()[2]-reference_pore.center_of_geometry()[2])
                        elif system == "ch1":
                            z_coordinate.append(-1*(cation.atoms.center_of_geometry()[2]-reference_pore.center_of_geometry()[2]))
                            
                        r_coordinate.append(np.sqrt((cation.atoms.center_of_geometry()[0]-reference_pore.center_of_geometry()[0])**2 + (cation.atoms.center_of_geometry()[1]-reference_pore.center_of_geometry()[1])**2))

                d = {'Z':z_coordinate,'R':r_coordinate}
                df = pd.DataFrame(data=d)
                df_list.append(df)

            print("############################################################################################\nSCRIPT FINISHED ON %s FOR %s %s %s\n############################################################################################\n" %(time.asctime().upper(),protein.upper(), construct.upper(), ion_solution.upper()))

            df_cat = pd.concat(df_list)
            df_cat['Z'] = bin_round(df_cat['Z'],bin_value)
            df_cat['R'] = bin_round(df_cat['R'],bin_value)
            df_cat["Count"] = df_cat['Z']
            df_histogram = df_cat.groupby(['Z','R']).count().reset_index()
            df_histogram = df_histogram.sort_values(by=['Z','R'],ascending=False)

            # Making an empty df...
            z_list = []
            r_list = []
            c_list = []
            for z in np.arange(-40,40.000001,bin_value):
                for r in np.arange(0,5.000001,bin_value):
                    z_list.append(z)
                    r_list.append(r)
                    c_list.append(0)
            d = {"Z":z_list, "R":r_list, "Count":c_list}
            df_empty = pd.DataFrame(data=d)

            df_merge = df_empty.merge(df_histogram, on=["Z","R"], how="outer", suffixes=("_l",""))
            cols = ["Z","R","Count"]
            df_histogram = df_merge[cols]
            df_histogram = df_histogram.fillna(0)
            df_histogram = df_histogram.sort_values(by=["Z","R"],ascending=False)



            df_histogram_file = f"count_analysis_2D_CE_{voltage}_{protein}_{construct}_{ion_solution}_{ionic_species}.csv"
            df_histogram_path = os.path.join(result_dir,df_histogram_file)
            df_histogram.to_csv(os.path.join(result_dir,df_histogram_path), index=False)

    print("############################################################################################\n############################################################################################\nSCRIPT COMPLETELY FINISHED ON %s\n############################################################################################\n############################################################################################\n" %time.asctime().upper())

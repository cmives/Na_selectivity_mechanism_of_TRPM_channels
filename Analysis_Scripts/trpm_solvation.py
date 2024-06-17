#!/usr/bin/env python3
#$ -cwd
#$ -V
#$ -N m5m_sol_na
#$ -e /cluster/uz_lab/cmives/TRP_channels/MDS/permeation_analysis/solvation_analysis/m5_solv_2024.err
#$ -o /cluster/uz_lab/cmives/TRP_channels/MDS/permeation_analysis/solvation_analysis/m5_solv_2024_na.out
#$ -pe smp 1
#$ -q c6100.q

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
result_dir = "/cluster/uz_lab/cmives/TRP_channels/MDS/permeation_analysis/solvation_analysis"
v_tm = 0.35
desired_time_step = 0.2      # in ns
bin_value = 0.25              # in A

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

ion_solution = ["kcl"]#,"cacl2"]#,"mixed"]
for ion_solution in ion_solution:
    # Specifiying the names of the ionic species...
    # Distances = International Journal of Thermophysics, Vol. 28, No. 2, April 2007, DOI: 10.1007/s10765-007-0154-6, Solvation Structure of Ions in Water, Raymond D. Mountain 
    if ion_solution == "cacl2":
        ionic_species = ["CAM"]
        ion_radius = 3.25

    elif ion_solution == "nacl":
        ionic_species = ["NA"]
        ion_radius = 3.29

    elif ion_solution == "kcl":
        ionic_species = ["K"]
        ion_radius = 3.57

    elif ion_solution == "mixed": 
        ionic_species = ["CAM","NA"]

    df_list = []
    replicate = ["repl_a","repl_b","repl_c"]
    for replicate in replicate:

        print("############################################################################################\nSCRIPT STARTED ON %s FOR %s %s %s %s\n############################################################################################\n" %(time.asctime().upper(),protein.upper(), construct.upper(), ion_solution.upper(), replicate.upper()))

        ############################################################
        #  DEFINING VARIABLES BASED UPON THE INPUT VARIABLES #
        ############################################################
        voltage = str(int(v_tm*1000)) + "mv"

        dir_stem = os.path.join(trpmds,protein,forcefield,construct,ion_solution,"prodmd/applied_field",voltage,replicate)

        traj_ts = int(desired_time_step/0.2)

        ############################################################
        # IMPORTING THE COORDINATE AND TRAJECTORY FILES #
        ############################################################
        print("READING DATA FROM %s" %dir_stem)
        GRO = os.path.join(dir_stem,"prerun_confout.gro")
        XTC = GRO.replace("prerun_confout.gro", "traj_analysis.xtc")
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
        if ion_solution != "mixed":
            pore = u.select_atoms(f"(cyzone 4 40 -40 (backbone and (resid {upper_gate_resid} or resid {lower_gate_resid}))) and resname {ionic_species[0]}", updating=True)
        else:
            pore = u.select_atoms(f"(cyzone 4 40 -40 (backbone and (resid {upper_gate_resid} or resid {lower_gate_resid}))) and (resname {ionic_species[0]} or resname {ionic_species[1]})", updating=True)

        cation_resname = []
        z_coordinate = []
        solvation_number = []
        protein_solvation_number = []
        total_solvation_number = []
        for frame in u.trajectory[::traj_ts]:
            print(f"Analysing frame for {frame.time / 1000} ns")
            for cation in pore.residues:
                # Getting the name of the cation so I can split it if necessary...
                cation_resname.append(cation.resname)
                # Getting the z-coordinate of the cation...
                z_coordinate.append(cation.atoms.center_of_geometry()[2]-reference_pore.center_of_geometry()[2])
                # Getting the number of water molecules in the first solvation shell...
                solvation_shell = u.select_atoms(f"(sphzone {ion_radius} resid {cation.resid} and resname {cation.resname} and name {cation.atoms[0].name}) and name OH2")
                solvation_number.append(len(solvation_shell.residues))
                # Getting the number of protein contacts in the first solvation shell...
                protein_shell = u.select_atoms(f"(sphzone {ion_radius} resid {cation.resid} and resname {cation.resname} and name {cation.atoms[0].name}) and name O* and protein")
                protein_solvation_number.append(len(protein_shell.residues))
                # Getting the total number of contacts in the first solvation shell...
                total_shell = u.select_atoms(f"(sphzone {ion_radius} resid {cation.resid} and resname {cation.resname} and name {cation.atoms[0].name}) and name O*")
                total_solvation_number.append(len(total_shell.residues))


        d = {'Z':z_coordinate,'Nc_TIP3':solvation_number,'Nc_Protein':protein_solvation_number,'Nc_Total':total_solvation_number,'cation_resname':cation_resname}
        df = pd.DataFrame(data=d)
        cols = ['Z','Nc_TIP3','Nc_Protein','Nc_Total','cation_resname']
        df = df[cols]
        df_list.append(df)

        print("############################################################################################\nSCRIPT FINISHED ON %s FOR %s %s %s %s\n############################################################################################\n" %(time.asctime().upper(),protein.upper(), construct.upper(), ion_solution.upper(), replicate.upper()))

    df_cat = pd.concat(df_list)
    df_cat['Z'] = bin_round(df_cat['Z'],bin_value)
    if ion_solution != "mixed":
        df_mean = df_cat.groupby('Z').mean().reset_index()
        df_mean = df_mean.sort_values(by='Z',ascending=False)
        df_sem = df_cat.groupby('Z').sem().reset_index()
        df_sem = df_sem.sort_values(by='Z',ascending=False)
        df_sem_tip3 = df_sem['Nc_TIP3'].tolist()
        df_sem_prot = df_sem['Nc_Protein'].tolist()
        df_sem_tot = df_sem['Nc_Total'].tolist()
        df_solvation = df_mean
        df_solvation['TIP3_sem'] = df_sem_tip3
        df_solvation['Protein_sem'] = df_sem_prot
        df_solvation['Total_sem'] = df_sem_tot

        df_solvation_file = f"solvation_analysis_2024_{protein}_{construct}_{ion_solution}_{ionic_species[0]}.csv"
        df_solvation_path = os.path.join(result_dir,df_solvation_file)
        df_solvation.to_csv(os.path.join(result_dir,df_solvation_path), index=False)

    else:
        for s in range(len(ionic_species)):
            df_cat_ion = df_cat[df_cat['cation_resname'] == ionic_species[s]]
            df_mean = df_cat_ion.groupby('Z').mean().reset_index()
            df_mean = df_mean.sort_values(by='Z',ascending=False)
            df_sem = df_cat_ion.groupby('Z').sem().reset_index()
            df_sem = df_sem.sort_values(by='Z',ascending=False)
            df_sem_tip3 = df_sem['Nc_TIP3'].tolist()
            df_sem_prot = df_sem['Nc_Protein'].tolist()
            df_sem_tot = df_sem['Nc_Total'].tolist()
            df_solvation = df_mean
            df_solvation['TIP3_sem'] = df_sem_tip3
            df_solvation['Protein_sem'] = df_sem_prot
            df_solvation['Total_sem'] = df_sem_tot

            df_solvation_file = f"solvation_analysis_2024_{protein}_{construct}_{ion_solution}_{ionic_species[s]}.csv"
            df_solvation_path = os.path.join(result_dir,df_solvation_file)
            df_solvation.to_csv(os.path.join(result_dir,df_solvation_path), index=False)

print("############################################################################################\n############################################################################################\nSCRIPT COMPLETELY FINISHED ON %s\n############################################################################################\n############################################################################################\n" %time.asctime().upper())

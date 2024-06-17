#!/usr/bin/env python3
#$ -cwd
#$ -V
#$ -N SFlex
#$ -e /cluster/uz_lab/cmives/TRP_channels/MDS/permeation_analysis/sf_flexibility/SFlex.err
#$ -o /cluster/uz_lab/cmives/TRP_channels/MDS/permeation_analysis/sf_flexibility/SFlex.out
#$ -pe smp 12
#$ -jc short

import os
import sys
import itertools
import numpy as np
import pandas as pd
from scipy import stats
import MDAnalysis as mda

#################################################################################
# SPECIFYING THE PROTEIN OF INTEREST AND THE PATH'S TO THE RELEVANT DIRECTORIES #
#################################################################################

protein = ["trpm5"]
forcefield = "charmm"
trpmds = "/cluster/uz_lab/cmives/TRP_channels/MDS/"
result_dir = "/cluster/uz_lab/cmives/TRP_channels/MDS/permeation_analysis/sf_flexibility"
v_tm = 0.35
calculate_sf_area = True
analysis_time_step = 0.2   # in ns
calculate_bb_flexibility = False

################################################################################
# DEFINING METHODS #
###############################################################################

voltage = str(int(v_tm*1000)) + "mv"

analysis_ts = int(analysis_time_step/0.2)

# Taken from Stack Overflow answer by Jamie Bull
# https://stackoverflow.com/questions/12642256/find-area-of-polygon-from-xyz-coordinates

#unit normal vector of plane defined by points a, b, and c
def unit_normal(a, b, c):
    x = np.linalg.det([[1,a[1],a[2]],
         [1,b[1],b[2]],
         [1,c[1],c[2]]])
    y = np.linalg.det([[a[0],1,a[2]],
         [b[0],1,b[2]],
         [c[0],1,c[2]]])
    z = np.linalg.det([[a[0],a[1],1],
         [b[0],b[1],1],
         [c[0],c[1],1]])
    magnitude = (x**2 + y**2 + z**2)**.5
    return (x/magnitude, y/magnitude, z/magnitude)

#area of polygon poly
def poly_area(poly):
    if len(poly) < 3: # not a plane - no area
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i+1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return abs(result/2)

def xvg_to_list(file_name):
    x, y = [], []
    with open(file_name) as f:
        next(f)
        for line in f:
            cols = line.split()
            if len(cols) == 2:
                x.append(float(cols[0]))
                y.append(float(cols[1]))
        return x, y

sf_res_names = ["alpha","beta","gamma"]

################################################################################
# STARTING THE ANALYSIS OF THE TRAJECTORIES #
###############################################################################

summary_file = os.path.join(result_dir,f"sf_area_summaries_trpm5.txt")
try:
    os.remove(summary_file)
except FileNotFoundError:
    pass
f = open(summary_file, 'w')

for protein in protein:

    if protein == 'trpm5':
        construct = "7MBS"
        sf_res_list = ["Q906","G905","F904"]

    #####################################################################################################

    replicate = ["repl_a","repl_b","repl_c"]
    df_area_list = [ [] for i in sf_res_list ]
    df_bb_list = []

    for ion_solution in ["nacl","cacl2"]:

        for replicate in replicate:
            ############################################################
            # IMPORTING THE COORDINATE AND TRAJECTORY FILES #
            ############################################################

            dir_stem = os.path.join(trpmds,protein,forcefield,construct,f"{ion_solution}/prodmd/applied_field",voltage,replicate)

            ############################################################
            # CALCULATING SF AREA #
            ############################################################
            if calculate_sf_area == True:
                print("\nREADING DATA FROM %s" %dir_stem)
                GRO = os.path.join(dir_stem,"prerun_confout.gro")
                XTC = GRO.replace("prerun_confout.gro", "traj_analysis.xtc")
                print("STARTING TO LOAD TRAJECTORY")
                u = mda.Universe(GRO,XTC)
                print("TRAJECTORY LOADED\n")

                for sf_res, n in zip(sf_res_list,range(len(sf_res_list))):
                    # Selecting the top residue in the SF....
                    if sf_res[0] == "Q":
                        sf = u.select_atoms(f"protein and resid {sf_res[1:]} and (name NE2 or name OE1)")
                    else:
                        sf = u.select_atoms(f"protein and resid {sf_res[1:]} and name O")


                    timestep = []
                    sf_area = []

                    for frame in u.trajectory[::analysis_ts]:
                        # Making a list for our timesteps...
                        timestep.append(frame.time/1000)
                        possible_sf_areas = []
                        # I iterate over all chain orders so that I don't have to manually state the order of the chains. This is to facilitate a wider use of this code on all TRPV structures...
                        for permutation in list(itertools.permutations([0,1,2,3])):
                            ordered_cog = []
                            for i in permutation:
                                ordered_cog.append(sf.residues[i].atoms.center_of_geometry())
                            possible_sf_areas.append(poly_area(ordered_cog))
                            #print(ordered_cog)
                        sf_area.append(max(possible_sf_areas))


                    d = {"timestep":timestep,"sf_area":sf_area}
                    df_area = pd.DataFrame(data=d)
                    df_area.to_csv(os.path.join(result_dir,f"sf_area_timecourse_{protein}_{construct}_{ion_solution}_{replicate}_{sf_res_names[n]}.csv"), index=False)
                    df_area_list[n].append(df_area)

            ############################################################
            # CALCULATING BACKBONE FLEXIBILITY #
            ############################################################
            if calculate_bb_flexibility == True:
                bb_data_repl = xvg_to_list(os.path.join(dir_stem,"rmsf_sf.xvg"))
                d = {"resnum":bb_data_repl[0],"rmsf":bb_data_repl[1]}
                df_bb = pd.DataFrame(data=d)
                # Converting from nm to A...
                df_bb['rmsf'] = df_bb['rmsf'] * 10
                df_bb['resnum'] = df_bb['resnum'].astype('int64')
                df_bb_list.append(df_bb)

        ############################################################
        # CALCULATING AVERAGES #
        ############################################################
        if calculate_sf_area == True:
            for n in range(len(df_area_list)):
                df_area_concat = pd.concat(df_area_list[n])
                average_area = np.mean(df_area_concat['sf_area'].tolist())
                sem_area = stats.sem(df_area_concat['sf_area'].tolist())
                print(f"Average area of SF {sf_res_names[n]}-residue gate for {protein.upper()} is {round(average_area,2)} +/- {round(sem_area,3)} A^2\n")

                f.write(f"Average area of SF {sf_res_names[n]}-residue gate for {protein.upper()} is {round(average_area,2)} +/- {round(sem_area,3)} A^2\n\n")

        if calculate_bb_flexibility == True:
            df_bb_concat = pd.concat(df_bb_list)
            df_bb_concat_sem = df_bb_concat.copy()
            df_bb_concat = df_bb_concat.groupby("resnum").mean()
            df_bb_concat_sem = df_bb_concat_sem.groupby("resnum").sem()
            df_bb_concat["rmsf_sem"] = df_bb_concat_sem["rmsf"]
            df_bb_concat = df_bb_concat.sort_values("resnum",ascending=False)
            df_bb_concat.to_csv(os.path.join(result_dir,f"sfbb_rmsf_{protein}_{construct}_{ion_solution}.csv"))

    f.close()


    #####################################################################################################
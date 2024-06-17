#!/usr/bin/env python3
#$ -cwd
#$ -V
#$ -N m_na_chap
#$ -e /cluster/uz_lab/cmives/TRP_channels/MDS/permeation_analysis/pore_architecture/m_na_chap.err
#$ -o /cluster/uz_lab/cmives/TRP_channels/MDS/permeation_analysis/pore_architecture/m_na_chap.err
#$ -pe smp 2
#$ -q c6100.q

import os
import sys
import glob
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 

#################################################################################
# SPECIFYING THE PROTEIN OF INTEREST AND THE PATH'S TO THE RELEVANT DIRECTORIES #
#################################################################################

proteins = ["trpm5"]
forcefield = "charmm"
trpmds = "/cluster/uz_lab/cmives/TRP_channels/MDS/"
result_dir = "/cluster/uz_lab/cmives/TRP_channels/MDS/permeation_analysis/pore_architecture"
v_tm = 0.35
perform_chap = True


################################################################################
# DEFINING METHODS AND OTHER VARIABLES #
################################################################################

voltage = str(int(v_tm*1000)) + "mv"


################################################################################
# STARTING THE ANALYSIS OF THE TRAJECTORIES #
################################################################################

comparison_list = [ [], [], [] ]

for protein in proteins:

    if protein == "trpm5":
        construct = "7MBS"


    #####################################################################################################

    ion_solutions = ["nacl","cacl2"]


    for ion_solution, ion_solution_index in zip(ion_solutions,range(len(ion_solutions))):

        if perform_chap == True:

            ###############################################################################
            # RUNNING CHAP ANALYSIS #
            ###############################################################################
            
            os.chdir(os.path.join(trpmds,protein,forcefield,construct,ion_solution,"prodmd/applied_field",voltage,"density"))

            # Running chap on overlapping windows...
            subprocess.run(f"chap -s chap.tpr -f traj_analysis_concat.xtc -n chap.ndx -tu ns -sel-pathway Protein -sel-solvent ion_chap -pf-max-free-dist 0.7 -pf-sel-ipp gates_chap -de-method histogram -de-res 0.025 -out-filename chap_profile_{protein}_{construct}_{ion_solution}_concat", shell=True, check=True)

            # Converting the json output to a csv using the provided chap_json2csv.py script...                  
            os.system(f"/cluster/uz_lab/cmives/miniconda3/envs/chap/chap/scripts/plotting/Python/chap_json2csv_cmi_loop.py -i chap_profile_{protein}_{construct}_{ion_solution}_concat.json")

            # Renaming the csv and pdb file, and deleting the unwanted output files...
            for parameter in ["reproducibilityInformation","pathwaySummary","pathwayProfile","pathwayScalarTimeSeries","pathwayProfileTimeSeries","residueSummary"]:
                os.rename(f"{parameter}.csv", os.path.join(result_dir,f"chap_{parameter}_{protein}_{construct}_{ion_solution}_concat.csv"))

            os.rename(f"chap_profile_{protein}_{construct}_{ion_solution}_concat.pdb", os.path.join(result_dir,f"chap_profile_{protein}_{construct}_{ion_solution}_concat.pdb"))
            os.rename(f"chap_profile_{protein}_{construct}_{ion_solution}_concat.obj", os.path.join(result_dir,f"chap_profile_{protein}_{construct}_{ion_solution}_concat.obj"))
            os.rename(f"chap_profile_{protein}_{construct}_{ion_solution}_concat.json", os.path.join(result_dir,f"chap_profile_{protein}_{construct}_{ion_solution}_concat.json"))
            os.rename(f"chap_profile_{protein}_{construct}_{ion_solution}_concat.mtl", os.path.join(result_dir,f"chap_profile_{protein}_{construct}_{ion_solution}_concat.mtl"))

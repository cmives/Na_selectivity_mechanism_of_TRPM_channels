#!/usr/bin/env python3
#$ -cwd
#$ -V
#$ -N trpm_perm
#$ -e /cluster/uz_lab/cmives/TRP_channels/MDS/permeation_analysis/permeation_identification/trpm_perm.err
#$ -o /cluster/uz_lab/cmives/TRP_channels/MDS/permeation_analysis/permeation_identification/trpm_perm.out
#$ -pe smp 2
#$ -q c6100.q

import os
import sys
import time
import numpy as np
import pandas as pd
import MDAnalysis as mda

#################################################################################
# SPECIFYING THE PROTEIN OF INTEREST AND THE PATH'S TO THE RELEVANT DIRECTORIES #
#################################################################################

protein = ["trpm5"]
forcefield = "charmm"
trpmds = "/cluster/uz_lab/cmives/TRP_channels/MDS/"
result_dir = "/cluster/uz_lab/cmives/TRP_channels/MDS/permeation_analysis/permeation_identification"
v_tm = 0.2
desired_time_step = 1    # in ns
analysis_time_step = 0.2 # in ns

################################################################################
# DEFINING METHODS #
################################################################################

def count_permeations(lst, seq):
     count = 0
     len_seq = len(seq)
     upper_bound = len(lst)-len_seq+1
     for i in range(upper_bound):
         if lst[i:i+len_seq] == seq:
             count += 1
     return count

################################################################################
# STARTING THE ANALYSIS OF THE TRAJECTORIES #
################################################################################
for protein in protein:
    if protein == "trpm5":
        construct = "7MBS"

    ion_solution = ["cacl2"]
    for ion_solution in ion_solution:
        replicate = ["repl_a","repl_b","repl_c"]
        for replicate in replicate:

            print("############################################################################################\nSCRIPT STARTED ON %s FOR %s %s %s %s\n############################################################################################\n" %(time.asctime().upper(),protein.upper(), construct.upper(), ion_solution.upper(), replicate.upper()))

            ############################################################
            #  DEFINING VARIABLES BASED UPON THE INPUT VARIABLES #
            ############################################################
            # Specifiying the names of the ionic species...
            if ion_solution == "cacl2":
                ionic_species = ["CAM","CL"]
            elif ion_solution == "nacl":
                ionic_species = ["NA","CL"]
            elif ion_solution == "kcl":
                ionic_species = ["K","CL"]
            elif ion_solution == "mixed": 
                ionic_species = ["CAM","NA","CL"]

            voltage = str(int(v_tm*1000)) + "mv"

            dir_stem = os.path.join(trpmds,protein,forcefield,construct,ion_solution,"prodmd/applied_field",voltage,replicate)

            traj_ts = int(desired_time_step/0.2)
            analysis_ts = int(analysis_time_step/0.2)

            ############################################################
            # IMPORTING THE COORDINATE AND TRAJECTORY FILES #
            ############################################################
            print("READING DATA FROM %s" %dir_stem)
            GRO = os.path.join(dir_stem,"prerun_confout.gro")
            XTC = GRO.replace("prerun_confout.gro", "traj_analysis.xtc")
            print("STARTING TO LOAD TRAJECTORY")
            u = mda.Universe(GRO,XTC)
            print("TRAJECTORY LOADED\n")              

            traj_step = u.trajectory[::traj_ts]           
            traj_analysis = u.trajectory[::analysis_ts]   

            # Adding a check to make sure that the timestep is as wanted...
            #if (traj_analysis[1].time - traj_analysis[0].time)*1000 != analysis_time_step:
                #sys.exit("INCORRECT TIMESTEP!!!")

            ############################################################
            # IDENTIFYING IONS THAT ENTER THE PORE #
            ############################################################
            #Specifying the pore region, so we can identify ions that are within the pore. This selection is set to update, so will update for every frame in the trajectory. I'm going to define the pore as the cylinder between the upper and lower gate.

            if protein == 'trpm5':
                upper_gate_resid = 906  # Gln
                lower_gate_resid = 966  # Ile


            pore = u.select_atoms(f"cyzone 10 15 -15 (backbone and (resid {upper_gate_resid} or resid {lower_gate_resid})) and not (resname TIP3 or resname POPC)", updating=True) 
            lower_gate = u.select_atoms(f"backbone and resid {lower_gate_resid}")
            upper_gate = u.select_atoms(f"backbone and resid {upper_gate_resid}")

            # Iterating through the trajectory, and identifying the cation/anion residue ID's found within the pore at each timestep...
            pore_ionic_species = [ [] for i in range(len(ionic_species)) ]
            for frame in traj_step:
                print("Scanning frame for %s ns" %(frame.time / 1000))
                for ion, s in zip(ionic_species,range(len(ionic_species))):
                    for residue in pore.residues:
                        if residue.resname == ion:
                            if residue.resid not in pore_ionic_species[s]:
                                pore_ionic_species[s].append(str(residue.resid))

            pore_ionic_species = [list(set(i)) for i in pore_ionic_species]
            print("\n")
            for s in range(len(ionic_species)):
                print(f"Number of pore-interacting {ionic_species[s]} ions = {len(pore_ionic_species[s])}")
            print("\n")

            ############################################################
            # IDENTIFYING IONS THAT ENTER THROUGH THE PORE #
            ############################################################
            # Iterating through the trajectory to get the z-position information for every pore-interacting cation and anion... 
            time_steps = []
            ion_traj_x = [ [] for i in range(len(ionic_species)) ]
            ion_traj_y = [ [] for i in range(len(ionic_species)) ]
            ion_traj_z = [ [] for i in range(len(ionic_species)) ]
            lower_gate_x = []
            lower_gate_y = []
            lower_gate_z = []
            upper_gate_x = []
            upper_gate_y = []
            upper_gate_z = []
            for frame in traj_analysis:
                print(f"Analysing frame for {frame.time / 1000} ns") #RENAME THIS!!!
                time_steps.append(frame.time/1000) # This is to use in making a Pandas DataFrame later in the script...
                pore_CoG = pore.center_of_geometry()
                for s in range(len(ionic_species)):
                    ion_frame_x = []
                    ion_frame_y = []
                    ion_frame_z = []
                    # Getting the x, y, and z coordinates for the pore-interacting cations at a given frame...
                    for ion in pore_ionic_species[s]:
                        ion_resid = ion
                        ion = u.select_atoms(f"resid {ion_resid}")
                        ion_x = ion.center_of_geometry()[0] - pore_CoG[0]
                        ion_frame_x.append(ion_x)
                        ion_y = ion.center_of_geometry()[1] - pore_CoG[1]
                        ion_frame_y.append(ion_y)
                        ion_z = ion.center_of_geometry()[2] - pore_CoG[2]
                        ion_frame_z.append(ion_z)
                    ion_traj_x[s].append(ion_frame_x)
                    ion_traj_y[s].append(ion_frame_y)
                    ion_traj_z[s].append(ion_frame_z)


                # Getting the coordinates for references of the gates...
                upper_x =  upper_gate.center_of_geometry()[0] - pore_CoG[0]
                upper_gate_x.append(upper_x)
                upper_y =  upper_gate.center_of_geometry()[1] - pore_CoG[1]
                upper_gate_y.append(upper_y)
                upper_z =  upper_gate.center_of_geometry()[2] - pore_CoG[2]
                upper_gate_z.append(upper_z)
                lower_x =  lower_gate.center_of_geometry()[0] - pore_CoG[0]
                lower_gate_x.append(lower_x)
                lower_y =  lower_gate.center_of_geometry()[1] - pore_CoG[1]
                lower_gate_y.append(lower_y)
                lower_z =  lower_gate.center_of_geometry()[2] - pore_CoG[2]
                lower_gate_z.append(lower_z)

            # To make a Pandas DataFrame of the permeation events, I want the trajectory timesteps for the x-axis...
            df_template = pd.DataFrame(time_steps,columns=['Timestep'])

            # Going to make and save the gates DataFrame here...
            d = {'Timestep':time_steps, 'Upper_x':upper_gate_x, 'Upper_y':upper_gate_y, 'Upper_z':upper_gate_z, 'Lower_x':lower_gate_x, 'Lower_y':lower_gate_y, 'Lower_z':lower_gate_z}
            df_gates = pd.DataFrame(data=d)
            cols = ['Timestep', 'Upper_x', 'Upper_y', 'Upper_z', 'Lower_x', 'Lower_y', 'Lower_z']
            df_gates = df_gates[cols]
            df_name = f"gates_relative_coordinates_{protein}_{construct}_{ion_solution}_{replicate}_200mv.csv"
            df_path = os.path.join(result_dir,df_name)
            df_gates.to_csv(df_path,index=False)

            ###############################################################
            # IDENTIFYING IONS THAT ENTER FULLY PERMEATE THROUGH THE PORE #
            ###############################################################

            print("\n")

            # The permeation events are calculated based on identifying a sequence of the permeation state. If the ion is above the upper gate it is assigned a state of 1, inside the pore a state of 0, and below the pore a state of -1. By then removing adjacent duplicates, identification of a sequence of 1,0,-1 equates to a successful permeation event.
            permeating_ion = [ [] for i in range(len(ionic_species)) ]
            permeating_time_ion = [ [] for i in range(len(ionic_species)) ]
            permeating_length = [ [] for i in range(len(ionic_species)) ]
            permeating_ion_x_position = [ [] for i in range(len(ionic_species)) ]
            permeating_ion_y_position = [ [] for i in range(len(ionic_species)) ]
            permeating_ion_z_position = [ [] for i in range(len(ionic_species)) ]
            df_list = []
            for ion, s in zip(ionic_species,range(len(ionic_species))):
                for i, n in zip(range(len(time_steps)),range(len(pore_ionic_species[s]))):
                    ion_z_position = [value[i] for value in ion_traj_z[s]]
                    ion_resid = pore_ionic_species[s][n]
                    perm_state = []              
                    for a in range(len(time_steps)):
                        if ion_z_position[a] >= upper_gate_z[a]:
                            perm_state.append(1)
                        elif upper_gate_z[a] >= ion_z_position[a] >= lower_gate_z[a]:
                            perm_state.append(0)
                        elif lower_gate_z[a] >= ion_z_position[a]:
                            perm_state.append(-1) 
                    d = {'Timesteps':time_steps,'State':perm_state}
                    df = pd.DataFrame(data=d)
                    # Removing neighbouring duplicates...
                    df2 = df.loc[df['State'].shift() != df['State']]
                    perm_state = df2['State'].to_list()
                    time_steps_mod = df2['Timesteps'].tolist()
                    if perm_state[0] == 0:
                        perm_count_init = count_permeations(perm_state[:2], [0,-1]) # This is to include ions that entered the pore during equilibration...
                        perm_count = count_permeations(perm_state, [1,0,-1]) + perm_count_init
                    elif perm_state[0] != 0:
                        perm_count = count_permeations(perm_state, [1,0,-1])
                    if perm_count != 0:

                        permeating_ion_z_position[s].append(ion_z_position)
                        ion_x_position = [value[i] for value in ion_traj_x[s]]
                        permeating_ion_x_position[s].append(ion_x_position)
                        ion_y_position = [value[i] for value in ion_traj_y[s]]
                        permeating_ion_y_position[s].append(ion_y_position)

                        df_ion = df_template
                        column_name = str(ion_resid) + "_x"
                        df_ion[column_name] = ion_x_position
                        column_name = str(ion_resid) + "_y"
                        df_ion[column_name] = ion_y_position
                        column_name = str(ion_resid) + "_z"
                        df_ion[column_name] = ion_z_position
                        df_list.append(df_ion)   

                        for p in range(len(perm_state)-2):
                            if p == 0:
                                if perm_state[p] == 0 and perm_state[p+1] ==-1:
                                    permeating_time_ion[s].append(time_steps_mod[p+1])
                                    permeating_ion[s].append(ion_resid)
                                    permeating_length[s].append(np.NaN)
                            else:
                                if perm_state[p] == 1 and perm_state[p+1] == 0 and perm_state[p+2] == -1:
                                    permeating_time_ion[s].append(time_steps_mod[p+2])
                                    permeating_ion[s].append(ion_resid)
                                    permeating_length[s].append(round((time_steps_mod[p+2]-time_steps_mod[p+1]),2))


                # Writing the residue ID's of the permeating ions to a file so they can be easily used in other scripts...
                ion_resid_file = f"permeating_ion_resids_{protein}_{construct}_{ion_solution}_{ionic_species[s]}_{replicate}_200mv.csv"
                ion_resid_path = os.path.join(result_dir,ion_resid_file)
                d = {'permeating_ion':permeating_ion[s],'permeation_time':permeating_time_ion[s]}
                permeating_ion_df = pd.DataFrame(data=d)
                permeating_ion_df = permeating_ion_df.sort_values(by='permeation_time')
                permeating_ion_df.to_csv(os.path.join(result_dir,ion_resid_path), index=False)
                # Saving the permeation event times for calculation of permeation frequency and conductances...
                np.savetxt(f"{result_dir}/permeation_times_{protein}_{construct}_{ion_solution}_{ionic_species[s]}_{replicate}_200mv.csv", sorted(permeating_time_ion[s]), delimiter=",", fmt='%s', header='Permeation_times') 

                print(f"Number of permeating {ionic_species[s]} ions = {len(permeating_time_ion[s])}")

                # Saving the permeation duration times...
                np.savetxt(f"{result_dir}/permeation_duration_{protein}_{construct}_{ion_solution}_{ionic_species[s]}_{replicate}_200mv.csv", sorted(permeating_length[s]), delimiter=",", fmt='%s', header='Permeation_duration') 


            ############################################################
            # ACCOUNTING FOR PBC EFFECTS #
            ############################################################
            # Due to the species moving from one side of the box to another due to PBC, this can make plotting the data diffcult. Therefore if a pbc is spotted, I correct for this by replacing with a None value. This DataFrame is saved separately however...
            for df, s in zip(df_list,range(len(ionic_species))):
                df2_ion = df.iloc[:-1]
                all_ion_z_pbc = []
                for ion in permeating_ion[s]:
                    column_name = str(ion) +"_z"
                    t = df_ion[column_name].tolist()
                    ion_z_pbc = []
                    # Identifying if the ion has jumped from one "side" of the protein to the other...
                    for n in range(len(t)-1):
                        if t[n] < -30 and t[n+1] > 30:
                            # If so, we replace that value with "None"...
                            ion_z_pbc.append(None)
                        elif t[n] > 30 and t[n+1] < -30:
                            ion_z_pbc.append(None)
                        else:
                            ion_z_pbc.append(float(t[n]))
                    all_ion_z_pbc.append(ion_z_pbc)
                    column_name = str(ion) + "_z"
                    df2_ion.loc[:,column_name] = ion_z_pbc
                df_name = f"permeation_identification_pbc_{protein}_{construct}_{ion_solution}_{ionic_species[s]}_{replicate}_200mv.csv"
                df_path = os.path.join(result_dir,df_name)
                df2_ion.to_csv(df_path,index=False)


            ############################################################
            # 
            ############################################################

            print("############################################################################################\nSCRIPT FINISHED ON %s FOR %s %s %s %s\n############################################################################################\n" %(time.asctime().upper(),protein.upper(), construct.upper(), ion_solution.upper(), replicate.upper()))


print("############################################################################################\n############################################################################################\nSCRIPT COMPLETELY FINISHED ON %s\n############################################################################################\n############################################################################################\n" %time.asctime().upper())

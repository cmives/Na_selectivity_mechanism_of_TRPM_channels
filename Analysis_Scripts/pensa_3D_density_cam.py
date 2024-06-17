#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#$ -N density_cam 
#$ -e proderr.err
#$ -o prodout.out
#$ -l gpu=1
#$ -pe smp 4
#$ -cwd
#$ -R y
#$ -jc long
#$ -mods l_hard m_mem_free 16G
#$ -mods l_hard h_rt 72:00:00
#$ -notify

"""
Created on Fri Sep 23 14:28:27 2022

@author: neil
"""

import os
import sys
sys.path.append('/homes/njthomson/pensa/')

from pensa import *

import numpy as np
from scipy import ndimage as ndi
import os
from gridData import Grid
import MDAnalysis as mda
from MDAnalysis.analysis.density import DensityAnalysis
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import align
import biotite.structure as struc
import biotite.structure.io as strucio
from pensa.features.processing import *
from tqdm import tqdm
            


def ca_grid(structure_input, xtc_input, top_atoms=35, out_name="CAM"):
    
    """
    Write out water pockets for the top X most probable atoms (top_atoms).

    Parameters
    ----------
    structure_input : str
        File name for the reference file (PDB or GRO format).
    xtc_input : str
        File name for the trajectory (xtc format).
    atomgroup : str
        Atomgroup selection to calculate the density for (atom name in structure_input).
    top_atoms : int, optional
        Number of atoms to featurize. The default is 35.
    grid_input : str, optional
        File name for the density grid input. The default is None, and a grid is automatically generated.
    write : bool, optional
        If true, a reference pdb will be written out. The default is None.
    write_grid_as : str, optional
        If you choose to write out the grid, you must specify the water model 
        to convert the density into. The default is None. Options are suggested if default.
    out_name : str, optional
        Prefix for all written filenames. The default is None.

    Returns
    -------
        feature_names : list of str
            Names of all features
        features_data : numpy array
            Data for all features

    """

        
    u = mda.Universe(structure_input, xtc_input)
    
    if not os.path.exists('dens/'):
        os.makedirs('dens/')
    p = u.select_atoms("protein")
    pdb_outname = 'dens/' + out_name +"Sites.pdb"
    p_avg = np.zeros_like(p.positions)
    # do a quick average of the protein (in reality you probably want to remove PBC and RMSD-superpose)
    for ts in u.trajectory:
        p_avg += p.positions
    p_avg /= len(u.trajectory)
    # temporarily replace positions with the average
    # p.load_new(p_avg)
    p.positions = p_avg
    # write average protein coordinates
    p.write(pdb_outname)
    # just make sure that we have clean original coordinates again (start at the beginning)
    u.trajectory.rewind()    
    
    density_atomgroup = u.select_atoms("resname CAM and around 10 protein", updating=True)
    # a resolution of delta=1.0 ensures the coordinates of the maxima match the coordinates of the simulation box
    D = DensityAnalysis(density_atomgroup, delta=1.0)
    D.run(verbose=True)
    g = D.density
    D.density.convert_density("A^{-3}")
    D.density.export('dens/' + out_name + "_density.dx", type="double")

    xyz, val = local_maxima_3D(g.grid)
    ## Negate the array to get probabilities in descending order
    val_sort = np.argsort(-1*val.copy())
    newvals = [val[max_val] for max_val in val_sort]  
    coords = [xyz[max_val] for max_val in val_sort]    
    maxdens_coord_str = [str(item)[1:-1] for item in coords]
    atom_information=[]
    atom_dists=[]

    if top_atoms > len(coords):
        top_atoms = len(coords)  


    print('\n')
    print('Featurizing ',top_atoms,' Atoms')
    for at_no in tqdm(range(top_atoms)):
        print('\n')
        print('Atom no: ',at_no+1)
        print('\n')

        ## Find all water atoms within 3.5 Angstroms of density maxima
        # Shifting the coordinates of the maxima by the grid origin to match 
        # the simulation box coordinates
        shifted_coords=coords[at_no]+g.origin
        point_str = str(shifted_coords)[1:-1]
        densval = newvals[at_no]

        atom_ID = "A" + str(at_no+1)
        atom_location = shifted_coords

        atom_information.append([atom_ID,list(atom_location),densval])
        
        ## Write data out and visualize water sites in pdb           
        write_CA_to_pdb(pdb_outname, atom_location, atom_ID)
        u_pdb = mda.Universe(pdb_outname)
        u_pdb.add_TopologyAttr('tempfactors')
        # Write values as beta-factors ("tempfactors") to a PDB file
        for res in range(len(atom_information)):
            #scale the atom resid by the starting resid
            atom_resid = len(u_pdb.residues) - at_no-1 + res
            u_pdb.residues[atom_resid].atoms.tempfactors = atom_information[res][-1]
        u_pdb.atoms.write(pdb_outname)
    
    # Return the dictionaries.
    return print('Pdb file completed.')



def write_CA_to_pdb(pdb_outname, atom_location, atom_ID):
    """
    Write a new atom to a reference structure to visualise conserved non-protein atom sites.

    Parameters
    ----------
    pdb_outname : str
        Filename of reference structure.
    atom_location : array
        (x,y,z) coordinates of the atom location with respect to the reference structure.
    atom_ID : str
        A unique ID for the atom.
    atomgroup : str
        MDAnalysis atomgroup to describe the atom.

    """
    
    ##PDB_VISUALISATION     
    ##rescursively add waters to the pdb file one by one as they are processed           
    # # Read the file into Biotite's structure object (atom array)
    atom_array = strucio.load_structure(pdb_outname)
    res_id = atom_array.res_id[-1] + 1
    # Add an HETATM
    atom = struc.Atom(
        coord = atom_location,
        chain_id = "X",
        # The residue ID is the last ID in the file +1
        res_id = res_id,
        res_name = atom_ID,
        hetero = True,
        atom_name = "CA",
        element = "CA"
        )
    atom_array += struc.array([atom])
    # Save edited structure
    strucio.save_structure(pdb_outname, atom_array)


"""
##############################################################################
##############################################################################
"""
input_gro_1 = 'step6.6_equilibration.gro'
input_xtc_1 = 'traj_analysis_concat.xtc'


# # # Then we featurize the waters common to both simulations
# # # We can do the same analysis for ions using the get_atom_features featurizer. 
atom_feat_a, atom_data_a = ca_grid(structure_input = input_gro_1, 
                                   xtc_input = input_xtc_1)

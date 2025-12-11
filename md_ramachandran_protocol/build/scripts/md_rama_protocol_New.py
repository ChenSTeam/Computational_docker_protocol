from rdkit import Chem
import pandas as pd
from tqdm import tqdm
import peptide_generation
import md_simulation
import md_analysis
import os
import argparse
import traceback
import datetime
from Bio.PDB import PDBParser
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--uaa_file', type=str, help='input uaa file')
    parser.add_argument('--output', type=str, default='.', help='output folder')
    parser.add_argument('--timestep', type=float, default=0.002, help='* unit of picosecond')
    parser.add_argument('--ministep', type=int, default=500000)
    parser.add_argument('--steps', type=int, default=50000000)
    parser.add_argument('--savesteps', type=int, default=50000)

    # Add
    parser.add_argument('--structure_file', type = str, help = 'input file with peptide structure')

    
    return parser.parse_args()

def main():
    # get inputs
    args = parse_args()

    uaa_file = args.uaa_file
    out_folder = args.output
    log_file = f'{out_folder}/log.txt'
    error_output = f'{out_folder}/error.txt'
    ministep = args.ministep
    stepnum = args.steps
    timestep = args.timestep
    savestep = args.savesteps
    # Add
    structure  = args.structure_file
    if (uaa_file && structure):
        print('Please input either uaa_file or structure_file')
        break
    if (uaa_file):
        uaa_matrix = pd.read_csv(uaa_file)
        uaa_name = uaa_matrix.name.tolist()
        uaa_smiles = uaa_matrix.smiles.tolist()
        uaa_list = list(zip(uaa_name, uaa_smiles))

        if not os.path.exists(f'{out_folder}/md'):
            os.makedirs(f'{out_folder}/md')
        if not os.path.exists(f'{out_folder}/rama'):
            os.makedirs(f'{out_folder}/rama')
        with open(log_file, 'w+') as f:
            f.write('uaa_name,uaa_smiles,status\n')

        for uaa_name, uaa_smiles in tqdm(uaa_list, total=len(uaa_list)):

            ligand_path = f'{out_folder}/md/{uaa_name}.sdf'
            topology_file = f'{out_folder}/md/{uaa_name}.cif'
            traj_output = f'{out_folder}/md/{uaa_name}.dcd'
            log_output = f'{out_folder}/md/{uaa_name}_log.txt'
            phi_psi_out = f'{out_folder}/rama/{uaa_name}.txt'

            try:
                #generate tripeptide
                peptide = peptide_generation.uaa_tripeptide_generate(uaa_smiles)
                peptide = peptide_generation.generate_3d_struct(peptide)
                Chem.MolToMolFile(peptide, ligand_path)
            
                #simulation of the tripeptide    
                md_simulation.md_run(ligand_path, topology_file, traj_output, log_output, ministep, stepnum , timestep, savestep)
                md_simulation.recenter_molecule(traj_output, topology_file)
                md_simulation.remove_solvent(traj_output, topology_file)
            
                #analysis of ramachandran properties
                md_analysis.phi_psi_calc(traj_output,topology_file,out=phi_psi_out)
                md_analysis.rama_plot(uaa_name, phi_psi_out)

            except Exception as e:
                if not os.path.exists(error_output):
                    with open(error_output, 'w+') as f:
                         f.write('\n')

                with open(log_file, 'a+') as f:  
                    f.write(f'{uaa_name},{uaa_smiles},failed\n')
            
                timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                error_traceback = traceback.format_exc()
                error_content = (
                    "========================================\n"
                    f"time: {timestamp}\n"
                    f"error type: {type(e).__name__}\n"
                    f"error messages: {e}\n"
                    "------\n"
                    f"{error_traceback}"
                    "========================================\n\n"
                    )

                with open(error_output, 'a+') as f:
                     f.write(f'{uaa_name},{uaa_smiles},failed\n')
                     f.write(error_content)

                continue
        
            with open(log_file, 'a+') as f:
                    f.write(f'{uaa_name},{uaa_smiles},completed\n')

    if (structure):
        structure_matrix = pd.read_csv(structure)
        structure_name = structure_matrix.name.tolist()
        parser = PDBParser(QUIET=True)
        
        if not os.path.exists(f'{out_folder}/md'):
            os.makedirs(f'{out_folder}/md')
        if not os.path.exists(f'{out_folder}/rama'):
            os.makedirs(f'{out_folder}/rama')
        with open(log_file, 'w+') as f:
            f.write('structure_name,status\n')

        for structure_name in tqdm(structure_name, total=len(structure_name)):
            
            ligand_path = f'{out_folder}/md/{structure_name}.sdf'
            topology_file = f'{out_folder}/md/{structure_name}.cif'
            traj_output = f'{out_folder}/md/{structure_name}.dcd'
            log_output = f'{out_folder}/md/{structure_name}_log.txt'
            phi_psi_out = f'{out_folder}/rama/{structure_name}.txt'

            try:
                file_path = f'{structure}/{structure_name}/{structure_name}.pdb'
                protein =  Chem.MolFromPDBFile(file_path, removeHs=True, sanitize=True)
                Chem.MolToMolFile(protein, ligand_path)

                #simulation of the protein  
                md_simulation.md_run(ligand_path, topology_file, traj_output, log_output, ministep, stepnum , timestep, savestep)
                md_simulation.recenter_molecule(traj_output, topology_file)
                md_simulation.remove_solvent(traj_output, topology_file)
            
                #analysis of ramachandran properties
                md_analysis.phi_psi_calc(traj_output,topology_file,out=phi_psi_out)
                md_analysis.rama_plot_structure(structure_name, phi_psi_out)

            except Exception as e:
                if not os.path.exists(error_output):
                    with open(error_output, 'w+') as f:
                         f.write('\n')

                with open(log_file, 'a+') as f:  
                    f.write(f'{structure_name},failed\n')
            
                timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                error_traceback = traceback.format_exc()
                error_content = (
                    "========================================\n"
                    f"time: {timestamp}\n"
                    f"error type: {type(e).__name__}\n"
                    f"error messages: {e}\n"
                    "------\n"
                    f"{error_traceback}"
                    "========================================\n\n"
                    )

                with open(error_output, 'a+') as f:
                     f.write(f'{structure_name},failed\n')
                     f.write(error_content)

                continue
        
            with open(log_file, 'a+') as f:
                    f.write(f'{structure_name},completed\n')
            
if __name__ == '__main__':
    main()

from rdkit import Chem
import pandas as pd
from tqdm import tqdm

import peptide_generation
import md_simulation
import md_analysis
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--uaa_file', type=str, help='input uaa file')
    parser.add_argument('--output', type=str, default='.', help='output folder')
    parser.add_argument('--timestep', type=float, default=0.002, help='* unit of picosecond')
    parser.add_argument('--ministep', type=int, default=50000)
    parser.add_argument('--steps', type=int, default=50000000)
    parser.add_argument('--savesteps', type=int, default=50000)
    
    return parser.parse_args()

def main():
    # get inputs
    args = parse_args()

    uaa_file = args.uaa_file
    out_folder = args.output
    log_file = f'{out_folder}/log.txt'
    ministep = args.ministep
    stepnum = args.steps
    timestep = args.timestep
    savestep = args.savesteps

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

        except:
            with open(log_file, 'a+') as f:
                f.write(f'{uaa_name},{uaa_smiles},failed\n')
            continue
        
        with open(log_file, 'a+') as f:
                    f.write(f'{uaa_name},{uaa_smiles},completed\n')

if __name__ == '__main__':
    main()

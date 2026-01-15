import mdtraj as md

from sys import stdout

# OpenMM imports
import openmm.app as app
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

# OpenFF-toolkit imports
from openff.toolkit import Molecule
from openff.toolkit import Topology as offTopology
from openff.units.openmm import to_openmm as offquantity_to_openmm

import parmed as pmd

import subprocess

import argparse

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--traj_file', type=str, help='input traj file')
    parser.add_argument('--top_file', type=str, help='input topology file')
    parser.add_argument('--ligand_path', type=str, help='input ligand file')
    parser.add_argument('--ligand_resname', type=str, help='ligand resname in topology file')
    parser.add_argument('--outputfolder', type=str, default='.', help='output_folder')
    parser.add_argument('--title', type=str, default='result', help='title for the output files')
    parser.add_argument('--start_frame', type=int, default=1, help='start frame in traj')
    parser.add_argument('--end_frame', type=int, default=100, help='end frame in traj')
    parser.add_argument('--inter_frame', type=int, default=5, help='interval')
    
    return parser.parse_args()

def PBSA_calc(traj_file, top_file, ligand_path, ligand_resname, outputfolder, cond, sf, ef, intf):
    traj = md.load(traj_file, top=top_file)

    # remove all other atoms except protein
    complex_atoms = traj.topology.select(f'chainid 0 or resname {ligand_resname}')
    traj_complex = traj.atom_slice(complex_atoms)
    traj_complex.save('traj_complex.dcd')
    traj_complex[0].save('traj_complex.pdb')

    protein_atoms = traj.topology.select('chainid 0')
    traj_protein = traj.atom_slice(protein_atoms)
    traj_protein[0].save('traj_protein.pdb')

    ligand_atoms = traj.topology.select(f'resname {ligand_resname}')
    traj_ligand = traj.atom_slice(ligand_atoms)
    traj_ligand[0].save('traj_ligand.pdb')

    complex_pdb = app.PDBFile('traj_complex.pdb')
    protein_pdb = app.PDBFile('traj_protein.pdb')
    ff = app.ForceField('amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml')

    ligand = Molecule.from_file(ligand_path)
    smirnoff = SMIRNOFFTemplateGenerator(molecules=ligand)
    ff.registerTemplateGenerator(smirnoff.generator)

    ligand_off_topology = offTopology.from_molecules(molecules=[ligand])
    ligand_omm_topology = ligand_off_topology.to_openmm()
    ligand_positions = offquantity_to_openmm(ligand.conformers[0])

    modeller = app.Modeller(protein_pdb.topology, protein_pdb.positions)

    modeller.add(ligand_omm_topology, ligand_positions)

    system = ff.createSystem(modeller.topology)

    structure = pmd.openmm.load_topology(complex_pdb.topology, system=system, xyz=complex_pdb.positions)
    structure.save('complex.prmtop', overwrite=True)

    rec_atoms = traj.topology.select('chainid 0')  # chain A
    rec_traj = traj.atom_slice(rec_atoms)
    rec_traj[0].save('receptor.pdb')
    rec_pdb = app.PDBFile('receptor.pdb')
    modeller = app.Modeller(rec_pdb.topology, rec_pdb.positions)
    system = ff.createSystem(modeller.topology)
    rec_struct = pmd.openmm.load_topology(rec_pdb.topology, system=system, xyz=rec_pdb.positions)
    rec_struct.save('receptor.prmtop', overwrite=True)

    complex = pmd.load_file("complex.prmtop")
    complex.strip(f"!:{ligand_resname}") 
    complex.save("ligand.prmtop", overwrite=True)

    infile = "mmpbsa.in"

    mmpbsa_content = f"""&general
    startframe={sf},
    endframe={ef},
    interval={intf},
    verbose=1,
    /
    &decomp
    idecomp=1,
    dec_verbose=1,
    /
    &pb
    istrng=0.15,
    /
    """

    with open(infile, "w+") as f:
        f.write(mmpbsa_content)

    cmd = [
        "MMPBSA.py",
        "-O",
        "-i", "mmpbsa.in",
        "-cp", "complex.prmtop",
        "-rp", "receptor.prmtop",
        "-lp", "ligand.prmtop",
        "-y", "traj_complex.dcd",
        "-o", f"{outputfolder}/ligand_o_{cond}.dat",
        "-eo", f"{outputfolder}/ligand_eo_{cond}.dat",
        "-do", f"{outputfolder}/ligand_do_{cond}.dat",
        "-deo", f"{outputfolder}/ligand_deo_{cond}.dat",
    ]

    subprocess.run(cmd, capture_output=True, text=True)

def main():
    args = parse_args()
    traj = args.traj_file
    topo = args.top_file
    ligand_path = args.ligand_path
    ligand_resname = args.ligand_resname
    outputfolder = args.outputfolder
    cond = args.title
    sf = args.start_frame
    ef = args.end_frame
    intf = args.inter_frame

    # first chain as receptor
    PBSA_calc(traj, topo, ligand_path, ligand_resname, outputfolder, cond, sf, ef, intf)
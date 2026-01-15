import mdtraj as md

import openmm.app as app
import parmed as pmd

import subprocess

import argparse

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--traj_file', type=str, help='input traj file')
    parser.add_argument('--top_file', type=str, help='input topology file')
    parser.add_argument('--outputfolder', type=str, default='.', help='output_folder')
    parser.add_argument('--title', type=str, default='result', help='title for the output files')
    parser.add_argument('--start_frame', type=int, default=1, help='start frame in traj')
    parser.add_argument('--end_frame', type=int, default=100, help='end frame in traj')
    parser.add_argument('--inter_frame', type=int, default=5, help='interval')
    
    return parser.parse_args()

def PBSA_calc(traj_file, top_file, outputfolder, cond, sf, ef, intf):
    traj = md.load(traj_file, top=top_file)

    # remove all other atoms except protein
    protein_atoms = traj.topology.select('protein and not resname UNK')

    traj_protein = traj.atom_slice(protein_atoms)

    traj_protein.save('traj_protein.dcd')
    traj_protein[0].save('traj_protein.pdb')

    protein_pdb = app.PDBFile('traj_protein.pdb')
    ff = app.ForceField('amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml')
    modeller = app.Modeller(protein_pdb.topology, protein_pdb.positions)
    system = ff.createSystem(modeller.topology)

    structure = pmd.openmm.load_topology(protein_pdb.topology, system=system, xyz=protein_pdb.positions)

    structure.save('complex.prmtop', overwrite=True)
    #structure.save('complex.inpcrd', overwrite=True)
    #structure.save('traj.rst7', overwrite=True)

    rec_atoms = traj.topology.select('chainid 0')  # chain A
    rec_traj = traj_protein.atom_slice(rec_atoms)
    rec_traj[0].save('receptor.pdb')
    rec_pdb = app.PDBFile('receptor.pdb')
    modeller = app.Modeller(rec_pdb.topology, rec_pdb.positions)
    system = ff.createSystem(modeller.topology)
    rec_struct = pmd.openmm.load_topology(rec_pdb.topology, system=system, xyz=rec_pdb.positions)
    rec_struct.save('receptor.prmtop', overwrite=True)

    lig_atoms = traj.topology.select('chainid 1')  # chain B
    lig_traj = traj_protein.atom_slice(lig_atoms)
    lig_traj[0].save('ligand.pdb')
    lig_pdb = app.PDBFile('ligand.pdb')
    modeller = app.Modeller(lig_pdb.topology, lig_pdb.positions)
    system = ff.createSystem(modeller.topology)
    lig_struct = pmd.openmm.load_topology(lig_pdb.topology, system=system, xyz=lig_pdb.positions)
    lig_struct.save('ligand.prmtop', overwrite=True)

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
        "-y", "traj_protein.dcd",
        "-o", f"{outputfolder}/o_{cond}.dat",
        "-eo", f"{outputfolder}/eo_{cond}.dat",
        "-do", f"{outputfolder}/do_{cond}.dat",
        "-deo", f"{outputfolder}/deo_{cond}.dat",
    ]

    subprocess.run(cmd, capture_output=True, text=True)

def main():
    args = parse_args()
    traj = args.traj_file
    topo = args.top_file
    outputfolder = args.outputfolder
    cond = args.title
    sf = args.start_frame
    ef = args.end_frame
    intf = args.inter_frame

    # first chain as receptor, second chain as ligand
    PBSA_calc(traj, topo, outputfolder, cond, sf, ef, intf)
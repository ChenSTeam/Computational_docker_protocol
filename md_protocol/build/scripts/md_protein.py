from sys import stdout
import argparse
import os

# OpenMM imports
import openmm.app as app
import openmm as mm
import openmm.unit as unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

# OpenFF-toolkit imports
from openff.toolkit import Molecule
from openff.toolkit import Topology as offTopology
from openff.units.openmm import to_openmm as offquantity_to_openmm

from pdbfixer import PDBFixer

import mdtraj as md

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np

import torch

if torch.cuda.is_available():
    device = 'CUDA'
else:
    device = 'CPU'

def fix_pdb(protin, protout):
        # Load the PDB file with missing residues
        fixer = PDBFixer(filename=protin)

        # Find and add missing residues and atoms
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        app.PDBFile.writeFile(fixer.topology, fixer.positions, open(protout, 'w'))
        print(f'Fix the pdb and save to {protout}')

def md_run(protein_path, topology_file, traj_output, log_output, state_save, ministep, stepnum, timestep, savestep):
    #read files
    print('Reading files ...')
    protein_pdb = app.PDBFile(protein_path)

    print('Generate force filed and modeller ...')
    #force field
    ff = app.ForceField('amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml')

    #modeller generate
    modeller = app.Modeller(protein_pdb.topology, protein_pdb.positions)
    modeller.addHydrogens(ff)

    print('Solvate the modeller ...')
    # Solvate
    modeller.addSolvent(ff, padding=1.0*unit.nanometer, ionicStrength=0.15*unit.molar)

    print('Create system, integrator and prepare simulation ...')
    # Create the system, define the integrator, and create the simulation
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=1*unit.nanometer, constraints=app.HBonds)
    integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, timestep*unit.picoseconds)
    platform = mm.openmm.Platform.getPlatform(device)
    simulation = app.Simulation(modeller.topology, system, integrator, platform)

    # set the positions
    simulation.context.setPositions(modeller.positions)

    print("Minimizing energy ...")
    simulation.minimizeEnergy(maxIterations=ministep)

    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)

    #save topology files
    ##recommend save as cif file. if choose pdb file, error (ValueError: invalid literal for int() with base 10: ' A') will happen due to the hexadecimal, when the atom number is very high.
    positions = simulation.context.getState(positions=True).getPositions()
    app.PDBxFile.writeFile(simulation.topology, positions, open(topology_file, 'w'))

    simulation.reporters.append(app.DCDReporter(traj_output, savestep))

    simulation.reporters.append(app.StateDataReporter(log_output, savestep, step=True, totalEnergy=True,
            potentialEnergy=True, kineticEnergy=True, temperature=True, speed=True, volume=True, density=True))

    #simulation.reporters.append(app.StateDataReporter(stdout, savestep, step=True,
            #potentialEnergy=True, temperature=True, speed=True))

    simulation.reporters.append(app.CheckpointReporter(traj_output[:-4] + '_checkpnt.chk', savestep))

    print("Running simulation ...")
    try:
        simulation.step(stepnum)
    except:
        if state_save:
            print('Saving the state ...')
            simulation.loadCheckpoint(traj_output[:-4] + '_checkpnt.chk')
            simulation.saveState(state_save)
        print('Complete!')

    #save the state
    if state_save:
        print('Saving the state ...')
        simulation.saveState(state_save)
    print('Complete!')

def recenter_molecule(*args, save_format = '.dcd'):
    if len(args) == 1:
        traj_file = args[0]
        m = md.load(traj_file)
    elif len(args) == 2:
        traj_file = args[0]
        top_file = args[1]
        m = md.load(traj_file, top=top_file)
    else:
        print('Not available!')
        return 0
    protein = m.topology.guess_anchor_molecules()
    n = m.image_molecules(anchor_molecules=protein)
    os.remove(traj_file)
    n.save(traj_file[:-4] + save_format)

def md_analysis(topo, traj, out, duration):
    m = md.load(traj,top = topo)

    # plot RMSD of all protein atoms
    fig = plt.figure(figsize=(8, 3))
    rmsd = md.rmsd(m, m, frame=0, atom_indices=m.topology.select('protein'))
    plt.plot((m.time + 1) * duration, rmsd)
    plt.xlabel('time (ps)')
    plt.ylabel('RMSD (nm)')
    plt.title('Protein RMSD')
    fig.savefig(f"{out}/rmsd_protein.svg", format="svg", bbox_inches = 'tight')
    plt.close(fig)
    np.savetxt(f'{out}/rmsd_protein.txt', rmsd)

    # plot RMSD of all protein CA atoms
    fig = plt.figure(figsize=(8, 3))
    rmsd = md.rmsd(m, m, frame=0, atom_indices=m.topology.select('protein and name CA'))
    plt.plot((m.time + 1) * duration, rmsd)
    plt.xlabel('time (ps)')
    plt.ylabel('RMSD (nm)')
    plt.title('Protein CA RMSD')
    fig.savefig(f"{out}/rmsd_protein_ca.svg", format="svg", bbox_inches = 'tight')
    plt.close(fig)
    np.savetxt(f'{out}/rmsd_protein_ca.txt', rmsd)

    # plot RMSF of all protein CA
    rmsf = md.rmsf(m, m, frame=0, atom_indices=m.topology.select('protein and name CA'))
    fig = plt.figure(figsize=(8, 3))
    plt.plot(range(1,1 + len(rmsf)),rmsf)
    plt.xlabel('Residue')
    plt.ylabel('RMSF (nm)')
    plt.title('Protein CA RMSF')
    #plt.ylim(0.0, 5.0)
    fig.savefig(f"{out}/rmsf_protein_ca.svg", format="svg", bbox_inches = 'tight')
    plt.close(fig)
    np.savetxt(f'{out}/rmsf_protein_ca.txt', rmsf)

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--pdb', type=str, help='input pdb file')
    parser.add_argument('--output', type=str, help='output folder')
    parser.add_argument('--timestep', type=float, default=0.002, help='* unit of picosecond')
    parser.add_argument('--ministep', type=int, default=50000)
    parser.add_argument('--steps', type=int, default=500000)
    parser.add_argument('--savesteps', type=int, default=5000)
    
    return parser.parse_args()

def main():
    # get inputs
    args = parse_args()

    filename = os.path.splitext(os.path.basename(args.pdb))[0]
    output = args.output + '/' + filename

    if not os.path.exists(output):
        os.makedirs(output)

    protein_path = f'{output}/{filename}_fix.pdb'
    topology_file = f'{output}/topology_{filename}.cif'
    traj_output = f'{output}/traj_{filename}.dcd'
    log_output = f'{output}/log_{filename}.txt'
    save_file = f'{output}/state_{filename}.xml'

    print('## 01 run simulation!')
    fix_pdb(args.pdb, protein_path)
    md_run(protein_path, topology_file, traj_output, log_output, save_file, ministep = args.ministep, stepnum = args.steps, timestep = args.timestep, savestep = args.savesteps)

    print('## 02 run processing!')
    recenter_molecule(traj_output, topology_file)

    print('## 03 run simple analysis!')
    md_analysis(topo = topology_file, traj = traj_output, out = output, duration= args.timestep * args.savesteps)

if __name__ == '__main__':
    main()
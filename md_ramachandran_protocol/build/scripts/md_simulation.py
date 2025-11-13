import os
import GPUtil

# OpenMM imports
import openmm.app as app
import openmm as mm
import openmm.unit as unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

# OpenFF-toolkit imports
from openff.toolkit import Molecule
from openff.toolkit import Topology as offTopology
from openff.units.openmm import to_openmm as offquantity_to_openmm

import mdtraj as md

if len(GPUtil.getAvailable()) > 0:
    device = 'CUDA'
else:
    device = 'CPU'

def md_run(ligand_path, topology_file, traj_output, log_output, ministep, stepnum, timestep, savestep):
    #read files
    print('Reading files ...')
    ligand = Molecule.from_file(ligand_path)

    print('Generate force filed and modeller ...')
    #force field
    smirnoff = SMIRNOFFTemplateGenerator(molecules=ligand)
    ff = app.ForceField('amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml')
    ff.registerTemplateGenerator(smirnoff.generator)

    #modeller generate
    ligand_off_topology = offTopology.from_molecules(molecules=[ligand])
    ligand_omm_topology = ligand_off_topology.to_openmm()
    ligand_positions = offquantity_to_openmm(ligand.conformers[0])
    modeller = app.Modeller(ligand_omm_topology, ligand_positions)
    modeller.addHydrogens(ff)

    print('Solvate the modeller ...')
    # Solvate
    modeller.addSolvent(ff, padding=2.0*unit.nanometer, ionicStrength=0.15*unit.molar)

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

    print("Running simulation ...")
    try:
        simulation.step(stepnum)
    except:
        print('Simulation did not run appropriately!')
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

def remove_solvent(*args, save_format = '.dcd'):
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
    selection = m.topology.select('not water and not resname NA and not resname CL')
    m = m.atom_slice(selection)

    if len(args) == 1:
        os.remove(traj_file)
    elif len(args) == 2:
        os.remove(traj_file)
        os.remove(top_file)
    m.save_dcd(traj_file)
    m.save_cif(top_file) 
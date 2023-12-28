import argparse
import random
import ufedmm
import mdtraj
# Importing CV LIB

from ufedmm import cvlib 
import simtk.openmm as mm
from simtk.openmm import app
from simtk import openmm, unit
from sys import stdout

waters = ['spce', 'tip3p', 'tip4pew', 'tip5p']
parser = argparse.ArgumentParser()
parser.add_argument('--water', dest='water', help='the water model', choices=waters, default=None)
parser.add_argument('--ff', dest='ff', help='the pepdide force field', default='amber03')
parser.add_argument('--seed', dest='seed', help='the RNG seed', default=None)
parser.add_argument('--platform', dest='platform', help='the computation platform', default='Reference')
args = parser.parse_args()

seed = random.SystemRandom().randint(0, 2**31) if args.seed is None else args.seed

pdb = app.PDBFile('not3.pdb')

#**************** Force field ***************************

forcefield = app.ForceField('amber99sb.xml','tip3p.xml')


#********************************************************

system = forcefield.createSystem(pdb.topology,nonbondedMethod=app.PME,nonbondedCutoff=1.0*unit.nanometers,constraints=None, rigidWater=False,removeCMMotion=False)

#****************** RG *******************************

# ALPHA CARBONS (resid=1 to 20) 

alpha_carbons = [atom.index for atom in pdb.topology.atoms() if atom.name == 'CA']

system.Rg= ufedmm.CollectiveVariable('Rg', cvlib.RadiusOfGyration(alpha_carbons))


#****************** RMSD *******************************
print("Alpha_carbons",alpha_carbons)



system.RMSD = ufedmm.CollectiveVariable('RMSD',openmm.RMSDForce(pdb.positions,alpha_carbons))


#***************** RMSD Hydrophobic Core ********************

# Tryptophan(resid=6),Proline resid(12,17,18,19) taken to define RMSD Hydrophibic core 

hcore_atm = [atom.index for atom in pdb.topology.atoms() if atom.residue.name == 'PRO' or atom.residue.name == 'TRP']


system.RMSDHpcore = ufedmm.CollectiveVariable('RMSDHpcore',openmm.RMSDForce(pdb.positions,hcore_atm))

print("core_atm",hcore_atm)


#***************************************************************************

#******************* RMSD Helix ********************************************

# RMSD Alpha Helix "CA" (resid=2 to 8)

rmsd_helix_atoms = [atom.index for atom in pdb.topology.atoms() if atom.residue.index in range(1,8) and atom.name =='CA']
print("RMSD_Helix_atoms",rmsd_helix_atoms)
system.RMSDH = ufedmm.CollectiveVariable('RMSDH',openmm.RMSDForce(pdb.positions,rmsd_helix_atoms))

#*************** Extracting phi psi atom index *******************************

traj = mdtraj.load('not3.pdb')
phi_atoms = [list(map(int, atoms)) for atoms in mdtraj.compute_phi(traj)[0]]
psi_atoms = [list(map(int, atoms)) for atoms in mdtraj.compute_psi(traj)[0]]
dihedral_atoms_list     = phi_atoms[0:18]+psi_atoms[1:19]

#*************** Alpha Helicity *********************************************

reference_angles_list1  = list()
reference_angles_list1[0:18]   = [-1.0]*18
reference_angles_list1[18:36]  = [-0.82]*18

force = openmm.PeriodicTorsionForce()
for atoms, phi_ref in zip(dihedral_atoms_list, reference_angles_list1):
    force.addTorsion(*atoms, 1, phi_ref, 0.500)  
system.AlphaH= ufedmm.CollectiveVariable('AlphaH',force) 

#*********** Beta similarity *****************************************************

reference_angles_list2  = list()
reference_angles_list2[0:18]   = [-1.396]*18
reference_angles_list2[18:36]  = [2.618]*18
force1 = openmm.PeriodicTorsionForce()
for atoms1, phi_ref1 in zip(dihedral_atoms_list, reference_angles_list2):
    force1.addTorsion(*atoms1, 1, phi_ref1, 0.500)  # atoms is a list or tuple, n=1, k=0.5
system.BetaS= ufedmm.CollectiveVariable('BetaS',force1)


#************* Dihcor Cv **********************************************************

traj2 = mdtraj.load('not3.pdb')

phi2_atoms = [list(map(int, atoms)) for atoms in mdtraj.compute_phi(traj)[0]]
psi2_atoms = [list(map(int, atoms)) for atoms in mdtraj.compute_psi(traj)[0]]
phi2_atoms=phi_atoms[0:18]
psi2_atoms=psi_atoms[1:19]

psi2_atoms += psi2_atoms[:-1]
phi2_atoms += phi2_atoms[1:]

energy = '0.5*(1 + cos(dihedral(p1, p2, p3, p4) - dihedral(p5, p6, p7, p8)))'
dihcor_force = openmm.CustomCompoundBondForce(8, energy)

for phi, psi in zip(phi2_atoms,psi2_atoms):
    dihcor_force.addBond(phi + psi, [])

system.dihcor = ufedmm.CollectiveVariable('dihcor', dihcor_force)


#*************** Distance b/w donor and acceptor in salt bridge **************


atm1           = 164
atm2           = 241
dist_d_A       = "distance(p1,p2)"
distance_force = openmm.CustomCompoundBondForce(2, dist_d_A)
distance_force.addBond([atm1,atm2])
system.sbridge = ufedmm.CollectiveVariable('sbridge', distance_force)

#*************** End to End distance Cv *******************************

atm3          = 4
atm4          = 294

dist_E_to_E   = "distance(p1,p2)"
distance_force1 = openmm.CustomCompoundBondForce(2, dist_E_to_E)
distance_force1.addBond([atm3,atm4])
system.dist_EE = ufedmm.CollectiveVariable('dist_EE', distance_force1)

#*****************Effective Mass of Cv's ************************************
#dihcor_M     = system.dihcor.effective_mass(system, pdb.positions,cv_unit=unit.amu*(unit.nanometers/unit.radians)**2)
#print("dihcor_M",dihcor_M)

#dist_EE_M    = system.dist_EE.effective_mass(system,pdb.positions,cv_unit=unit.amu*(unit.nanometers/unit.angstrom)**2)
#print("dist_EE_M",dist_EE_M)

#sbridge_M    = system.sbridge.effective_mass(system,pdb.positions,cv_unit=unit.amu*(unit.nanometers/unit.angstrom)**2)
#print("sbridge_M",sbridge_M)

#AlphaH_M     = system.AlphaH.effective_mass(system,pdb.positions,cv_unit=unit.amu*(unit.nanometers/unit.radians)**2)
#print("AlphaH_M",AlphaH_M)


#BetaS_M     = system.BetaS.effective_mass(system,pdb.positions,cv_unit=unit.amu*(unit.nanometers/unit.radians)**2)
#print("BetaS_M",BetaS_M)

#RG_M         = system.Rg.effective_mass(system,pdb.positions,cv_unit=unit.amu*(unit.nanometers/unit.angstrom)**2)
#print("RG_M",RG_M)


#RMSD_M       =  system.RMSD.effective_mass(system,pdb.positions,cv_unit=unit.amu*(unit.nanometers/unit.angstrom)**2)
#print("RMSD_M",RMSD_M)


#RMSDH_M      = system.RMSDH.effective_mass(system,pdb.positions,cv_unit=unit.amu*(unit.nanometers/unit.angstrom)**2)
#print("RMSDH_M",RMSDH_M)


#RMSDHpcore_M = system.RMSDHpcore.effective_mass(system,pdb.positions,cv_unit=unit.amu*(unit.nanometers/unit.angstrom)**2)
#print("RMSDHpcore_M",RMSDHpcore_M)
#******************************************************************************



#*******************************
temp = 300*unit.kelvin
gamma = 10/unit.picoseconds
dt = 0.5*unit.femtoseconds

nsteps =100000000
#nsteps = 100000
mass1 =2.343*unit.amu*(unit.angstrom/unit.angstrom)**2
mass2 = 0.05*unit.amu*(unit.angstrom/unit.angstrom)**2
mass3 = 2.0*unit.amu*(unit.angstrom/unit.radians)**2
mass4 = 0.05*unit.amu*(unit.angstrom/unit.radians)**2
mass5 = 40.0*unit.amu*(unit.angstrom/unit.angstroms)**2 
mass6 = 1.0*unit.amu*(unit.angstrom)**2
mass7 = 2.0*unit.amu*(unit.angstrom/unit.angstrom)**2
mass8 = 1.5*unit.amu*(unit.angstrom/unit.angstrom)**2
mass9 =0.01*unit.amu*(unit.angstrom/unit.angstrom)**2
mass10=0.01*unit.amu*(unit.angstrom/unit.angstrom)**2
Ks1   = 38.0*unit.kilojoules_per_mole/unit.angstrom**2
Ks2   = 5.0*unit.kilojoules_per_mole/unit.angstroms**2
Ks3   = 36.0*unit.kilojoules_per_mole/unit.radian**2
Ks4   = 2000*unit.kilojoules_per_mole/unit.radian**2
Ks5   = 5000*unit.kilojoules_per_mole/unit.angstrom**2 
Ks6   = 20.0*unit.kilojoules_per_mole 
Ks7   = 8.0*unit.kilojoules_per_mole/unit.angstrom**2
Ks8   = 10.0*unit.kilojoules_per_mole/unit.angstrom**2
Ks9   = 10.0*unit.kilojoules_per_mole/unit.angstrom**2
Ks10  = 10.0*unit.kilojoules_per_mole/unit.angstrom**2
Kh    = 0.01*unit.kilojoules_per_mole/unit.angstrom**2
Ts    = 3000*unit.kelvin
limit=180*unit.degrees
RG_max = 20000*unit.angstrom
RG_min = 0*unit.angstrom
RMSD_max = 20000*unit.angstrom
RMSD_min = 0*unit.angstrom
hyd_min=0*unit.angstrom
hyd_max=10*unit.angstrom
sigma = 0.1
frequency=500
height = 2.0*unit.kilojoules_per_mole # 2.0 inival
bias_factor=2
enforce_gridless=False
deposition_period = 500
val=1.00*unit.angstrom
#val=0.40



#************ Defining Aux CV RMSD and RG *******************
potential1 ='0.5*Ks1*(S_Rg - Rg*100.0)^2'
constants1 ={'Ks1':Ks1}

#***********  Harmonic restraint ****************************
potential2 ='0.5*Ks2*(S_RMSD - (RMSD*500.0))^2 +0.5*Kh*(S_RMSD -(val*500.0))^2'
constants2 = {'Ks2': Ks2,'Kh':Kh,'val':val}


potential3 = '0.5*Ks3*(S_AlphaH - (AlphaH*10.0))^2'
constants3 ={'Ks3':Ks3}

potential4 = '0.5*Ks4*(S_BetaS - BetaS)^2'
constants4 ={'Ks4':Ks4}

potential6 = '0.5*Ks6*(S_dihcor - (dihcor*5.0))^2'
constants6 ={'Ks6':Ks6}

potential7 ='0.5*Ks7*(S_RMSDH - (RMSDH*200.0))^2'
constants7  = {'Ks7': Ks7}

potential8  = '0.5*Ks8*(S_RMSDHpcore - (RMSDHpcore*200.0))^2'
constants8  = {'Ks8': Ks8}

potential9  = '0.5*Ks9*(S_sbridge - (sbridge*20.0))^2'
constants9  = {'Ks9': Ks9}

potential10 = '0.5*Ks10*(S_dist_EE - (dist_EE*20.0))^2'
constants10 = {'Ks10': Ks10}
#*************************************************************
S_Rg         = ufedmm.DynamicalVariable('S_Rg',RG_min,RG_max, mass1,Ts,system.Rg,potential1,sigma=sigma,periodic=False,**constants1)

S_RMSD       = ufedmm.DynamicalVariable('S_RMSD',RMSD_min,RMSD_max, mass2,Ts,system.RMSD,potential2,sigma=None,periodic=False,**constants2)

S_AlphaH     = ufedmm.DynamicalVariable('S_AlphaH',-1800,1800, mass3,Ts,system.AlphaH,potential3,sigma=None,periodic=True,**constants3)

S_BetaS      = ufedmm.DynamicalVariable('S_BetaS',-1800,1800, mass4,Ts,system.BetaS,potential4,sigma=None,periodic=True,**constants4)

S_dihcor     = ufedmm.DynamicalVariable('S_dihcor',-1800,1800, mass6,Ts,system.dihcor,potential6,sigma=None,periodic=False,**constants6)

S_RMSDH      = ufedmm.DynamicalVariable('S_RMSDH',RMSD_min,RMSD_max, mass7,Ts,system.RMSDH,potential7,sigma=None,periodic=False,**constants7)

S_RMSDHpcore = ufedmm.DynamicalVariable('S_RMSDHpcore',RMSD_min,RMSD_max, mass8,Ts,system.RMSDHpcore,potential8,sigma=None,periodic=False,**constants8)

S_sbridge    = ufedmm.DynamicalVariable('S_sbridge',0.0, 10000.0 , mass9,Ts,system.sbridge,potential9,sigma=None,periodic=False,**constants9)

S_dist_EE    = ufedmm.DynamicalVariable('S_dist_EE',0.0, 10000.0 , mass10,Ts,system.dist_EE,potential10,sigma=None,periodic=False,**constants10)
#*****************************************************************************

ufed=ufedmm.UnifiedFreeEnergyDynamics([S_Rg,S_RMSD,S_RMSDH,S_RMSDHpcore,S_AlphaH,S_sbridge,S_dist_EE,S_dihcor], temp, height, frequency, bias_factor,enforce_gridless=enforce_gridless)
#ufed=ufedmm.UnifiedFreeEnergyDynamics([S_Rg,S_RMSD],temp,height, frequency, bias_factor,enforce_gridless=enforce_gridless)


ufedmm.serialize(ufed, 'ufed_object.yml')
for force in system.getForces():
   if not isinstance(force, (openmm.HarmonicBondForce, openmm.HarmonicAngleForce)):
       force.setForceGroup(1)
integrator = ufedmm.MiddleMassiveGGMTIntegrator(temp,40.0*unit.femtoseconds, dt, scheme='VV-Middle',respa_loops=[2,1],bath_loops=4,unroll_loops=False)
#integrator.set_extended_space_time_constants([250.0*unit.femtoseconds,10.0*unit.femtoseconds,80.0*unit.femtoseconds,80.0*unit.femtoseconds,60.0*unit.femtoseconds,40.0*unit.femtoseconds,40.0*unit.femtoseconds,310.0*unit.femtoseconds])
#integrator = ufedmm.RegulatedNHLIntegrator(temp, tau, friction, dt, respa_loops=[1, 1], regulation_parameter=1)
#integrator = ufedmm.GeodesicLangevinIntegrator(temp, gamma, dt)
#integrator.setRandomNumberSeed(seed)

### GPU ####### 
#platform = openmm.Platform.getPlatformByName('OpenCL')
#properties = {'CudaPrecision': 'double'}

platform = openmm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}

##############
simulation = ufed.simulation(pdb.topology,system, integrator, platform)
simulation.context.setPositions(pdb.positions)
#simulation.context.setVelocitiesToTemperature(temp, seed)
print('Minimization stats')
simulation.minimizeEnergy(tolerance=0.1*unit.kilojoules_per_mole, maxIterations=0)
print('Minimization ends')

simulation.context.setVelocitiesToTemperature(temp, seed)
output1 = ufedmm.Tee(stdout, 'COLVAR')
reporter1 = ufedmm.StateDataReporter(output1,10, step=True, multipleTemperatures=False,hillHeights=False ,variables=True,speed=False,separator='\t')

output2 = ufedmm.Tee(stdout, 'HILL')
reporter2 = ufedmm.StateDataReporter(output2,500, step=True, multipleTemperatures=False,hillHeights=True ,variables=True,speed=True,separator='\t')
simulation.reporters.append(reporter1)
simulation.reporters.append(reporter2)


simulation.reporters.append(app.DCDReporter('output.dcd',100000))

simulation.step(nsteps)


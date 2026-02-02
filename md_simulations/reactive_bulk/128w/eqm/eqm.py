import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np
from datetime import datetime as dtt
import evb.utils as utils
import evb.check_forces as cf
from evb.evb import EVBSystem
from myreporter import *
from openmmplumed import *
import json

starttime = dtt.now();print(starttime)

#####################Set steps
box=np.array([1.566,1.566,1.566])
box_wu = box*unit.nanometer
nbcut = 0.78 * unit.nanometer
switchdist = 0.68 * unit.nanometer

path = "inits"
path_to_xml = "inits"
jobname = 'nvtrun'
prename = 'restart' 
restartfile='restart.rst'

with open("%s/para.json"%(path),"r") as inf:
    V=json.load(inf)["VL"][0]

temp = 300 #in K
pressure = 1 #in bar
dt = 0.5 * unit.femtoseconds
ewaldtol = 1e-5
isNVE = False
isPBC = True
nprint_restart = 20000 ## dump restart.rst file every 10 ps
nprint = 100000 # dump trajectory every 50 ps
restart = False
nstep = int(6100000) # in total 3050 ps (25 ps eqm at with weak restrain at CV=-3.0, 25 ps eqm at CV=-3.0 kappa changing from 500 to 10000, 1 ns eqm at CV=-3.0 with kappa=10000, 2 ns then to collect rst files into config_bank, every 10ps one config, in total 200 configs) 


print("###reading topology###")
psf1 = app.CharmmPsfFile("./inits/m1.psf")
psf2 = app.CharmmPsfFile("./inits/m2.psf")
psf3 = app.CharmmPsfFile("./inits/off.psf")
ff1 = app.ForceField("%s/toppar.xml"%(path_to_xml))
ff2 = app.ForceField("%s/toppar.xml"%(path_to_xml))
ff3 = app.ForceField("%s/toppar.xml"%(path_to_xml), "%s/off.xml"%(path_to_xml))

sub1=utils.create_Mysystem(psf1,ff1,box_wu,nbcut,switchdist,ewaldtol,isPBC,isNVE,isoff=False)
sub2=utils.create_Mysystem(psf2,ff2,box_wu,nbcut,switchdist,ewaldtol,isPBC,isNVE,isoff=False)
suboff=utils.create_Mysystem(psf3,ff3,box_wu,nbcut,switchdist,ewaldtol,isPBC,isNVE,isoff=True)

evbsys = EVBSystem(sub1, sub2, suboff,isPBC)
system = evbsys.mergeSystem(V)

fgrps = utils.forcegroupify(system) ## for EVB energy decom
fgrps1 = utils.forcegroupify(sub1)  ## for state 1 energy decom
fgrps2 = utils.forcegroupify(sub2)  ## for state 2 energy decom

restraint = mm.CustomExternalForce('k * periodicdistance(x, y, z, x0, y0, z0)^2')
system.addForce(restraint) # for PBC system, add this to the system
restraint.addPerParticleParameter("k")
restraint.addPerParticleParameter("x0")
restraint.addPerParticleParameter("y0")
restraint.addPerParticleParameter("z0")
restraint.addParticle(5,(100.0 * unit.kilojoules_per_mole / unit.nanometer**2, box_wu[0] / 2, box_wu[1] / 2, box_wu[2] / 2))

script = """
#RESTART
UNITS LENGTH=A TIME=fs
HtoC: DISTANCE ATOMS=5,6
HtoO: DISTANCE ATOMS=5,1
diff: CUSTOM ARG=HtoC,HtoO FUNC=x-y PERIODIC=NO

restraint: MOVINGRESTRAINT ...
    ARG=diff
    STEP0=0          AT0=-3.0 KAPPA0=0.0
    STEP1=50000      AT1=-3.0 KAPPA1=500.0
    STEP2=100000      AT2=-3.0 KAPPA2=10000.0
    STEP3=2100000      AT3=-3.0  KAPPA3=10000.0
    STEP4=6100000      AT4=-3.0  KAPPA4=10000.0
...
PRINT ARG=HtoC,HtoO,diff,restraint.diff_cntr,restraint.diff_work STRIDE=20000 FILE=SMD.cv
#FLUSH STRIDE=10
"""


system.addForce(PlumedForce(script)) 



platform = mm.Platform.getPlatformByName('CUDA')
prop = dict(CudaPrecision='double')
platform.loadPluginsFromDirectory('/home/yutong/software/myOpenMM-plumed-plugin/openmm-plumed1.0_install/lib/plugins/')


simulation = utils.make_simulation(psf1,system,temp, dt, isNVE)
simulation1 = utils.make_simulation(psf1,sub1,temp,dt, isNVE)
simulation2  = utils.make_simulation(psf2,sub2,temp,dt, isNVE)

if not restart:
    crd = app.CharmmCrdFile('inits/m1.crd') # coord file of MeCl + OHR + 128w
    simulation.context.setPositions(crd.positions)
    simulation.context.computeVirtualSites()
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(temp*unit.kelvin)
    simulation.context.setTime(0)
else:
    # start from the equilibrated positions and velocities
    print ('###reading restart file###')
    restartFname="%s"%restartfile
    print('#Startting from %s'%restartFname)
    simulation.loadState('%s'%restartFname)
    simulation1.loadState('%s'%restartFname)
    simulation2.loadState('%s'%restartFname)
    simulation.context.setTime(0)




print ("######setting up reporters######")
dcd=app.DCDReporter(jobname+'.dcd', nprint)
simulation.reporters.append(dcd)
simulation.reporters.append(app.StateDataReporter(jobname+'.log',nprint,time=True, kineticEnergy=True, potentialEnergy=True,totalEnergy=True, temperature=True, speed=True, separator=','))
simulation.reporters.append(RestartReporter(nprint_restart,'./restart/' + jobname + '-restart.log',jobname))


print ("###initial energies###")
state = simulation.context.getState(getEnergy=True)


print('PE', state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole),'kJ/mol')
for k, v in utils.getEnergyDecomp(simulation.context, fgrps).items():
    print(k.__class__.__name__, v.value_in_unit(unit.kilojoules_per_mole),'kJ/mol')


for step in range(nstep):                                                                                            
    simulation.step(1) ## run without bias                                                                        
    ######################################################                                                        
    ###state1, 2 energy decomposition, use fgrps created for state 1,2 systems                                    
    ## load the current coordinates into state 1 and 2 simulations                                                
    #evbpos = state.getPositions(asNumpy=True);                                                                    
    #simulation1.context.setPositions(evbpos);simulation2.context.setPositions(evbpos)                                                
    #state1 = simulation1.context.getState(getPositions=True, getEnergy=True)                                      
    #state2 = simulation2.context.getState(getPositions=True, getEnergy=True)                                      
                                                                                                                  
    #utils.printDecomposedEnergies(state,simulation,fgrps,printHeaders=True)                                      
    #cf.get_coul_pair_decomposition(sub2,simulation2)                                                             
    #cf.nbfix_type_decomposition(sub2,simulation2)                                                                
    #if (step+1)%100000==0:                                      
        #print ("Step %d  Time %f ps"%(step+1,state.getTime()._value))
        #print("c1sq %.6f    c2sq:  %.6f   E1:  %.6f  E2V:  %.6f       E12: %.6f   EVB PE: %.6f\n"%(c[0]*c[0],c[1]*c[1],E11,E2V,E12,EEVB))
        
        #utils.printDecomposedEnergies(state1,simulation1,fgrps1,printHeaders=True)
        #utils.printDecomposedEnergies(state2,simulation2,fgrps2,printHeaders=False)                                   
       # print('\n')
        ### print type decomposed CustomNonbonded energies
        # cf.nbfix_type_decomposition(sub1,simulation1)
        # cf.nbfix_type_decomposition(sub2,simulation2)
        #print("\n============================================================")                                                                                                                 
        #cf.check_harmonic_bond_energy(simulation1,sub1,box,isPBC=False)                                              
        #cf.check_harmonic_bond_energy(simulation2,sub2,box,isPBC=False)                                              
                                                                                                                  
        #cf.nbfix_pair_decomposition(sub1,simulation1)                                                                
        #utils.check_EVB(simulation,simulation1,simulation2,V,box,isPBC,0,4,5)                                         
    #####################################################     
    


# write final state
state = simulation.context.getState(getEnergy=True)
fgrps = utils.forcegroupify(system)
print('\n\nFinal PE', state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole),'kJ/mol')
for k, v in utils.getEnergyDecomp(simulation.context, fgrps).items():
    print(k.__class__.__name__, v.value_in_unit(unit.kilojoules_per_mole),'kJ/mol')

# write final state
state = simulation.context.getState(getPositions=True,getVelocities=True,getEnergy=True)
with open(prename+'-current.rst', 'w') as f:
    f.write(mm.XmlSerializer.serialize(state))

endtime = dtt.now()
print(endtime-starttime)

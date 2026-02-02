
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import os

#sys.path.append('/scratch3/usr/chungchi/myOpenMM_classes/')
#from myreporter import *

class RestartReporter(object):
    def __init__(self, reportInterval, fname, jobname):
        self._reportInterval = reportInterval
        self.file = open(fname, "w")
        self.nrestart = 0
        self.jobname = jobname
        isexist = os.path.exists("./")
        if isexist == "False":
            os.system("mkdir restart")

    def __del__(self):
        self.file.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, True, True, False, True, False)

    def report(self, simulation, state):
        self.nrestart += 1
        time = state.getTime()
        with open('./restart/' + self.jobname + '-' + '%d'%(self.nrestart) + '.rst', 'w') as f:
            f.write(mm.XmlSerializer.serialize(state))
        self.file.write("%f %s \n"%(time.value_in_unit(unit.picosecond), './restart/' + self.jobname + '-' + '%d'%(self.nrestart) + '.rst'))
        self.file.flush()

class VelocityReporter(object):
    def __init__(self, reportInterval, fname):
        self._reportInterval = reportInterval
        self.file = open(fname, "w")

    def __del__(self):
        self.file.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, True, False, False, False)

    def report(self, simulation, state):
        time = state.getTime()
        self.file.write("%f \n"%(time.value_in_unit(picosecond)))
        vels = state.getVelocities().value_in_unit(nanometer/picosecond)
        for v in vels:
            self.file.write('%g %g %g\n' % (v[0], v[1], v[2]))

class ForceReporter(object):
    def __init__(self, reportInterval, fname):
        self._reportInterval = reportInterval
        self.file = open(fname, "w")

    def __del__(self):
        self.file.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, False, True, False, False)

    def report(self, simulation, state):
        time = state.getTime()
        self.file.write("%f \n"%(time.value_in_unit(picosecond)))
        forces = state.getForces().value_in_unit(kilojoules/mole/nanometer)
        for f in forces:
            self.file.write('%g %g %g\n' % (f[0], f[1], f[2]))

##usuage example:
##simulation.reporters.append(Force_Reporter(print_interval_force,'./restart/' + jobname + '-restart.log', jobname)) ##print frequency (steps), outputfile name, jobname(for restartreporter)

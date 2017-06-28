###############################################################
# MRSA Transmission Model                                     #
# Universal Contact Precaution Induced Change in Contact Rate #
# Parameter Sensitivity Sweep                                 #
###############################################################

# Module Imports
import os
import stochpy
import numpy as numpy
import random

workingdir = os.getcwd()

# General simulation parameters
random.seed(12567)
start_time = 0.0
end_time = 8760
n_runs = 1000

# Model output storage arrays
acquisitions = numpy.empty([n_runs,3])

# All parameters are uniform random draws +/- 10% of the value used in the main text

# Baseline Model Run
def Baseline(iteration):
    model = stochpy.SSA()
    model.Model(model_file='MRSA_Colonization.psc', dir=workingdir)
    model.ChangeParameter('rho',random.uniform(3.3232,4.9848))
    model.ChangeParameter('sigma',random.uniform(0.0432,0.0648))
    model.ChangeParameter('psi',random.uniform(0.0745144,0.1117716))
    model.ChangeParameter('theta',random.uniform(0.007592,0.011388))
    colonized_admits = random.uniform(0.06232,0.09348)
    model.ChangeParameter('nu_C',colonized_admits)
    model.ChangeParameter('nu_U',1-colonized_admits)
    model.ChangeParameter('iota',random.uniform(4.592,6.888))
    model.ChangeParameter('tau',random.uniform(1.9112,2.8668))
    model.Endtime(end_time)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    Incident = outcomes[4][0][-1]
    acquisitions[iteration,0] = Incident

# Reduced Contact
def Reduced(iteration):
    model = stochpy.SSA()
    model.Model(model_file='MRSA_Colonization.psc', dir=workingdir)
    model.ChangeParameter('rho',random.uniform(2.734994,4.10249))
    model.ChangeParameter('sigma',random.uniform(0.0432,0.0648))
    model.ChangeParameter('psi',random.uniform(0.0745144,0.1117716))
    model.ChangeParameter('theta',random.uniform(0.007592,0.011388))
    colonized_admits = random.uniform(0.06232,0.09348)
    model.ChangeParameter('nu_C',colonized_admits)
    model.ChangeParameter('nu_U',1-colonized_admits)
    model.ChangeParameter('iota',random.uniform(4.592,6.888))
    model.ChangeParameter('tau',random.uniform(1.9112,2.8668))
    model.Endtime(end_time)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    Incident = outcomes[4][0][-1]
    acquisitions[iteration,1] = Incident

# Efficient Contact
def Efficient(iteration):
    model = stochpy.SSA()
    model.Model(model_file='MRSA_Colonization.psc', dir=workingdir)
    model.ChangeParameter('rho',random.uniform(3.3232,4.9848))
    model.ChangeParameter('sigma',random.uniform(0.0432,0.0648))
    model.ChangeParameter('psi',random.uniform(0.0745144,0.1117716))
    model.ChangeParameter('theta',random.uniform(0.007592,0.011388))
    colonized_admits = random.uniform(0.06232,0.09348)
    model.ChangeParameter('nu_C',colonized_admits)
    model.ChangeParameter('nu_U',1-colonized_admits)
    model.ChangeParameter('iota',random.uniform(4.592,6.888))
    model.ChangeParameter('tau',random.uniform(1.572918,2.359376))
    model.Endtime(end_time)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    Incident = outcomes[4][0][-1]
    acquisitions[iteration,2] = Incident

for i in range(0,n_runs):
	print("*** Iteration %i of %i ***" % (i+1,n_runs))
	Baseline(i)
	Reduced(i)
	Efficient(i)
	
numpy.savetxt('ContactChangeOutcomes_MRSA_Sensitivity.csv',acquisitions,delimiter=','
,header="Baseline,Reduced,Efficient",comments='')

print("*************************")
print("***** Runs Complete *****")
print("*************************")
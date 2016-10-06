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
    model.ChangeParameter('rho',random.uniform(3.7386,4.5694))
    model.ChangeParameter('sigma',random.uniform(0.0486,0.0594))
    model.ChangeParameter('psi',random.uniform(0.0838287,0.1024573))
    model.ChangeParameter('theta',random.uniform(0.008541,0.010439))
    colonized_admits = random.uniform(0.07011,0.08569)
    model.ChangeParameter('nu_C',colonized_admits)
    model.ChangeParameter('nu_U',1-colonized_admits)
    model.ChangeParameter('iota',random.uniform(5.166,6.314))
    model.ChangeParameter('tau',random.uniform(2.1501,2.6279))
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
    model.ChangeParameter('rho',random.uniform(3.076868,3.760616))
    model.ChangeParameter('sigma',random.uniform(0.0486,0.0594))
    model.ChangeParameter('psi',random.uniform(0.0838287,0.1024573))
    model.ChangeParameter('theta',random.uniform(0.008541,0.010439))
    colonized_admits = random.uniform(0.07011,0.08569)
    model.ChangeParameter('nu_C',colonized_admits)
    model.ChangeParameter('nu_U',1-colonized_admits)
    model.ChangeParameter('iota',random.uniform(5.166,6.314))
    model.ChangeParameter('tau',random.uniform(2.1501,2.6279))
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
    model.ChangeParameter('rho',random.uniform(3.7386,4.5694))
    model.ChangeParameter('sigma',random.uniform(0.0486,0.0594))
    model.ChangeParameter('psi',random.uniform(0.0838287,0.1024573))
    model.ChangeParameter('theta',random.uniform(0.008541,0.010439))
    colonized_admits = random.uniform(0.07011,0.08569)
    model.ChangeParameter('nu_C',colonized_admits)
    model.ChangeParameter('nu_U',1-colonized_admits)
    model.ChangeParameter('iota',random.uniform(5.166,6.314))
    model.ChangeParameter('tau',random.uniform(1.769532,2.162762))
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
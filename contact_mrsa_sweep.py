###############################################################
# MRSA Transmission Model                                     #
# Universal Contact Precaution Induced Change in Contact Rate #
###############################################################

# Module Imports
import os
import stochpy
import numpy as numpy
import random

workingdir = os.getcwd()

# General simulation parameters
random.seed(23408)
start_time = 0.0
end_time = 8760
n_runs = 5000

# Model output storage arrays
acquisitions = numpy.empty([n_runs,4])

# Reduced Contact
# Here we assume that the 18.3% reduction in worker visits is simply a flat reduction
def Reduced(iteration):
    model = stochpy.SSA()
    model.Model(model_file='MRSA_Colonization.psc', dir=workingdir)
    percentchange = random.uniform(0.50,1)
    model.ChangeParameter('rho',(4.154*percentchange))
    model.Endtime(end_time)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    Incident = outcomes[4][0][-1]
    acquisitions[iteration,0] = (1-percentchange)*100
    acquisitions[iteration,1] = Incident

# Efficient Contact
# Here we assume the same *number* of patient care tasks, but done more efficiently
# This reduces the opportunities to wash hands on exit (more tasks/visit)
# Need to change iota

def Efficient(iteration):
    model = stochpy.SSA()
    model.Model(model_file='MRSA_Colonization.psc', dir=workingdir)
    percentchange = random.uniform(0.50,1)
    model.ChangeParameter('tau',(2.389*percentchange))
    model.Endtime(end_time)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    Incident = outcomes[4][0][-1]
    acquisitions[iteration,2] = (1-percentchange)*100
    acquisitions[iteration,3] = Incident
    

for i in range(0,n_runs):
	print("*** Iteration %i of %i ***" % (i+1,n_runs))
	Reduced(i)
	Efficient(i)
	
numpy.savetxt('ContactChangeOutcomes_MRSA_Sweep.csv',acquisitions,delimiter=','
,header="ReducedPer,ReducedCases,EfficientPer,EfficientCases",comments='')

print("*************************")
print("***** Runs Complete *****")
print("*************************")
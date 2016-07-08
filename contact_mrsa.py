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
random.seed(80209)
start_time = 0.0
end_time = 8760
n_runs = 1000

# Model output storage arrays
acquisitions = numpy.empty([n_runs,5])

# Baseline Model Run
def Baseline(iteration):
    model = stochpy.SSA()
    model.Model(model_file='MRSA_Colonization.psc', dir=workingdir)
    model.Endtime(end_time)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    Incident = outcomes[4][0][-1]
    acquisitions[iteration,0] = Incident

# Reduced Contact
# Here we assume that the 18.3% reduction in worker visits is simply a flat reduction
def Reduced(iteration):
    model = stochpy.SSA()
    model.Model(model_file='MRSA_Colonization.psc', dir=workingdir)
    model.ChangeParameter('rho',(4.154*0.817))
    model.Endtime(end_time)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    Incident = outcomes[4][0][-1]
    acquisitions[iteration,1] = Incident

# Efficient Contact
# Here we assume the same *number* of patient care tasks, but done more efficiently
# This reduces the opportunities to wash hands on exit (more tasks/visit)
def Efficient(iteration):
    model = stochpy.SSA()
    model.Model(model_file='MRSA_Colonization.psc', dir=workingdir)
    model.ChangeParameter('tau',(2.389*0.817))
    model.Endtime(end_time)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    Incident = outcomes[4][0][-1]
    acquisitions[iteration,2] = Incident
    
# Effect of Hand Hygiene Compliance
# Increase in outgoing hand-hygiene compliance in intervention arm
# No difference in entry, +15.4% on exit. Overall compliance = 0.6425

def RedHH(iteration):
    model = stochpy.SSA()
    model.Model(model_file='MRSA_Colonization.psc', dir=workingdir)
    model.ChangeParameter('rho',(4.154*0.817))
    model.ChangeParameter('iota',(6.520))
    model.Endtime(end_time)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    Incident = outcomes[4][0][-1]
    acquisitions[iteration,3] = Incident

def EffHH(iteration):
    model = stochpy.SSA()
    model.Model(model_file='MRSA_Colonization.psc', dir=workingdir)
    model.ChangeParameter('iota',6.520)
    model.ChangeParameter('tau',(2.389*0.817))
    model.Endtime(end_time)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    Incident = outcomes[4][0][-1]
    acquisitions[iteration,4] = Incident	

for i in range(0,n_runs):
	print("*** Iteration %i of %i ***" % (i+1,n_runs))
	Baseline(i)
	Reduced(i)
	Efficient(i)
	RedHH(i)
	EffHH(i)
	
numpy.savetxt('ContactChangeOutcomes_MRSA.csv',acquisitions,delimiter=','
,header="Baseline,Reduced,Efficient,RedHH,EffHH",comments='')

print("*************************")
print("***** Runs Complete *****")
print("*************************")
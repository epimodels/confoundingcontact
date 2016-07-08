###############################################################
# MRSA Transmission Model                                     #
# Universal Contact Precaution Induced Change in Contact Rate #
###############################################################

# Module Imports
import os
import stochpy
import numpy as numpy
import pylab as pl
import matplotlib as mpl
import random

workingdir = os.getcwd()

# General simulation parameters
random.seed(96820)
start_time = 0.0
end_time = 8760
n_runs = 1

# Load model file - currently absolute path
baseline = stochpy.SSA()
baseline.Model(model_file='MRSA_Colonization.psc', dir=workingdir)
baseline.Endtime(end_time)
baseline.DoStochSim()

trajectories = baseline.data_stochsim.getSpecies()
time = trajectories[:,0]
H = trajectories[:,1]
S = trajectories[:,2]
U = trajectories[:,3]
C = trajectories[:,4]
Acquisition = trajectories[:,5]

fig = pl.figure("PatientTrajectories", figsize=(15,10))
pl.ax = fig.add_subplot(211)
pl.plot(time,S,'black', linewidth=2)
pl.plot(time,H,'0.75',linewidth=2)
pl.axis([0,end_time,0,9])
pl.xlabel("Time (Hours)")
pl.ylabel("Healthcare Workers")
pl.legend(("Uncontaminated","Contaminated"), shadow=False, fancybox=False, ncol=2)

patient_axix_max = max(max(U),max(C)) * 1.6

pl.ax = fig.add_subplot(212)
pl.plot(time,U,'black', linewidth=2) 
pl.plot(time,C,'0.75', linewidth=2)
pl.axis([0,end_time,-0.5,patient_axix_max])
pl.xlabel("Time (Hours)")	
pl.ylabel("Patients")
pl.legend(("Uncolonized", "Colonized"), shadow=False, fancybox=False, ncol=2)

pl.show()

a = raw_input("Save Figure (Y/N)?")

if a == "Y":
	pl.savefig('stoch_patients.pdf')

else:
	pass






import os
import stochpy
import random
import numpy as numpy
from scipy import stats
import matplotlib.pyplot as plt

workingdir = os.getcwd()

# Simulation parameters
random.seed(83222)
start_time = 0.0
end_time = 8760
n_runs = 25

# Run is a single run of the model that returns the number of incident acquisitions
def Run(filename,k):
    model = stochpy.SSA()
    model.Model(model_file=filename, dir=workingdir)
    model.Endtime(end_time)
    model.ChangeParameter('psi',k)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    acquisitions = outcomes[4][0][-1]
    return acquisitions

# Batch collects the average outcome of n_runs of the model
def Batch(filename,k):
    batch_cases = numpy.empty([n_runs,1])
    for i in range(0,n_runs):
        batch_cases[i,0] = Run(filename,k)
        batch_avg = numpy.mean(batch_cases)
    return batch_avg

# Basic accept-reject algorithm, converts model outcomes to incidence rate
def AcceptReject(objective,result,tolerance):
    rate = result/(365*18)
    lnrate = numpy.log(rate)
    lnobjective = numpy.log(objective/(365*18))
    distance = abs(lnobjective - lnrate)
    if distance <= tolerance:
        accept = 1
    else:
        accept = 0
    return accept

# Main fitting function
def ABCFit(config,objective,tol,iterations,priorhi,priorlow):
    results = numpy.zeros([iterations,3])
    for i in range(0,iterations):
        draw = random.uniform(priorhi,priorlow)
        results[i,0] = draw
        sim_avg = Batch(config,draw)
        results[i,1] = sim_avg
        decision = AcceptReject(objective=objective,result=sim_avg,tolerance=tol)
        results[i,2] = decision
        print("*** Iteration %i of %i ***" % (i+1,iterations))
    return results

# Fitting Function - Tweak prior to improve acceptance rate    
MRSA_Fit = ABCFit(config='MRSA_Colonization.psc',objective=39.03,tol=0.25,iterations=1000,priorhi=0.20,priorlow=0)

# ABC Diagnostics and Outcome
accept_rate = numpy.mean(MRSA_Fit[:,2])
print("ABC Acceptance Rate: %f" % accept_rate)

psi_estimate = MRSA_Fit[:,0][MRSA_Fit[:,2]==1]
psi_central = numpy.median(psi_estimate)

psi_upper = numpy.percentile(psi_estimate,97.5)
psi_lower = numpy.percentile(psi_estimate,2.5)

print("psi Median (95%% CI): %f (%f, %f)" % (psi_central,psi_lower,psi_upper))

# Plotting the posterior estimate
density = stats.kde.gaussian_kde(psi_estimate)

x = numpy.arange(0,0.20,0.0001)
high_density = max(density(x))
plt.plot(x,density(x))
plt.plot((psi_central, psi_central), (0,high_density), 'r--')
a = plt.Rectangle((0, 0), 0.25, 0.25, fc="red", alpha=1.00)
b = plt.Rectangle((0, 0), 0.25, 0.25, fc="blue",alpha=1.00)
plt.legend([a,b], ["Median","Density Estimate"],loc=1,fontsize='small')
plt.xlabel("psi",fontsize=16)
plt.ylabel("Density",fontsize=16)

plt.savefig('psiEstimateDistribution.pdf')
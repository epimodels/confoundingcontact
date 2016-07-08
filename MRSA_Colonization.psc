###################################################
# Dynamic Transmission Model of MRSA		      #
# Queue-based Steady State Population             #
# Author: Eric Lofgren (Eric.Lofgren@gmail.com)   #
###################################################

# Descriptive Information for PML File
Modelname: MRSA
Description: PML Implementation of MRSA transmission model 

# Set model to run with numbers of individuals
Species_In_Conc: False
Output_In_Conc: False

### Model Reactions ###

# Reactions Governing Movement of Healthcare Workers #
R1:
	H > S
	H*iota

R2:
	H > S
	H*tau*(C/(C+U))
	
R3:
	S > H
	rho*sigma*C*(S/(S + H + C + U))

# Reactions Involving Uncontaminated, Low Risk Patients (U) #
R4:
	U > C + Acquisition
	rho*psi*U*(H/(S + H + C + U))
	
R5:
	U > U
	theta*U*nu_U

R6:	
	U > C
	theta*U*nu_C

	
# Reactions Involving Contaminated, Low Risk Patients (C) #
	
R7:
	C > C
	theta*C*nu_C

R8:
	C > U
	theta*C*nu_U

### Parameter Values ###
## Time Values are in HOURS ##
# Compartments #
S = 6
H = 1
U = 17
C = 1
Acquisition = 0

# Contact Rates and Contamination Probabilities
rho = 4.154 # direct care tasks per patient per hour
sigma = 0.054
psi = 0.093143 

# Fitting Results
#ABC Acceptance Rate: 0.189000
#psi Median (95% CI): 0.093143 (0.076011, 0.115808)

# Exit (death/discharge) rates
theta = 0.00949

# Admission Proportions
nu_C = 0.0779
nu_U = 0.9221

# Handwashing and Gown/Glove Change Rates
iota = 5.74 #10.682 direct care tasks per hour with 56.55% compliance and 95% efficacy
tau = 2.389 #2.89 gown/glove chances per hour with 82.66% compliance

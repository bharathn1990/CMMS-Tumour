'''File that stores all global variables and constants'''
import os	
import numpy as np
import random

#
# ******** GENERAL PARAMETERS **********
#

nSites = 101 		# no. of lattice sites in each dimension

#
# Determine Initial Concentration of O2
#

maxOxygen = 0.082
minOxygen = 0.07
ratio = round(random.random(),3)						# Generate random ratio.
initOxygen = minOxygen + (maxOxygen - minOxygen)*ratio 	# Initial concentration of oxygen everywhere
initOxygen = round(initOxygen,5)
print(initOxygen)

#
# Create the directory for the results.	
#

if not os.path.exists("results"):
	os.mkdir("results")
dirName = 'results/' + str(initOxygen)
if not os.path.exists(dirName):	
	os.mkdir(dirName)


#
# ******* MONTE CARLO PARAMETERS **********
#

J = 12 				# Adhesive energy constant
kbT = 4.183365e-21 	# constant that determines probability of change (Boltzmann x Temperature)


gammaP = 10			# volume energy constants for each type of cell 
gammaQ = 30
gammaM = 0
gammaN = 50


#
# ******* DIFFUSION PARAMETERS **********
#

dt = 9e-08
nIterDiff = 600

oxygenP, oxygenQ, oxygenN = 50,20,0	# Oxygen consumption for each type of cell

thresholdP2Q = 0.06					# Thresholds for cell conversion
thresholdQ2N = 0.045
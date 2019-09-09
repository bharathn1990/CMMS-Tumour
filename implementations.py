''' All the functions called by the main file run.py are defined here.'''
import random
from config import *
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import os
import imageio

def calculate_hamiltonian(cellInfo,cellIDs) :
	''' Function that calculates the Hamiltonian '''
	
	totalEnergy, adhesiveEnergy, volumeEnergy =0.0, 0 , 0
	
	# ************ ADHESIVE ENERGY ******************* #
		
	for i in range(1,nSites-1):
		for j in range(1,nSites-1):
			neighbours = cellIDs[i-1:i+2,j-1:j+2]
			nDifferent = 9 - len(neighbours[neighbours == cellIDs[i,j]])
			adhesiveEnergy += (nDifferent * J)/2 # We divide by 2 due to double counting
	
	# ************ VOLUME ENERGY ******************* #
	
	for key in cellInfo.keys():
		currentVolume = len(cellIDs[cellIDs == key])
		volumeEnergy += (cellInfo[key]['gamma']) * ((cellInfo[key]['targetVol'] - currentVolume))**2
	
	totalEnergy = adhesiveEnergy + volumeEnergy
	return totalEnergy;



def calculate_deltaH(cellInfo,cellIDs,posX,posY,oldCellID,newCellID) :
	''' Function that calculates the change in Hamiltonian '''	
	#
	# First calculate old energy contribution
	#
	
	# ************ ADHESIVE ENERGY ******************* #
	cellIDs[posX,posY] = oldCellID    
	oldAdhesiveEnergyContribution = 0
	for i in range(max(1,posX-1),min(posX+2,nSites-1)):		# loop only through immediate neighbours
		for j in range(max(1,posY-1),min(posY+2,nSites-1)):
			neighbours = cellIDs[i-1:i+2,j-1:j+2]
			
			iDifferent = neighbours[neighbours != cellIDs[i,j]]
			nDifferent = len(iDifferent)
			oldAdhesiveEnergyContribution += nDifferent*J/2

	# ************ VOLUME ENERGY ******************* #
	volumeOld = len(cellIDs[cellIDs == oldCellID])
	volumeNew = len(cellIDs[cellIDs == newCellID])
	oldVolumeEnergyContribution = (cellInfo[oldCellID]['gamma']) * ((cellInfo[oldCellID]['targetVol'] - volumeOld))**2
	oldVolumeEnergyContribution += (cellInfo[newCellID]['gamma']) * ((cellInfo[newCellID]['targetVol'] - volumeNew))**2
	
	#
	# Calculate new energy contribution
	#
	
	# ************ ADHESIVE ENERGY ******************* #
	cellIDs[posX,posY] = newCellID
	newAdhesiveEnergyContribution = 0
	for i in range(max(1,posX-1),min(posX+2,nSites-1)):
		for j in range(max(1,posY-1),min(posY+2,nSites-1)):
			neighbours = cellIDs[i-1:i+2,j-1:j+2]
			iDifferent = neighbours[neighbours != cellIDs[i,j]]
			nDifferent = len(iDifferent)
			newAdhesiveEnergyContribution += nDifferent*J/2

	# ************ VOLUME ENERGY ******************* #
	volumeOld = len(cellIDs[cellIDs == oldCellID])
	volumeNew = len(cellIDs[cellIDs == newCellID])
	
	newVolumeEnergyContribution = (cellInfo[oldCellID]['gamma']) * ((cellInfo[oldCellID]['targetVol'] - volumeOld))**2
	newVolumeEnergyContribution += (cellInfo[newCellID]['gamma']) * ((cellInfo[newCellID]['targetVol'] - volumeNew))**2
	
	deltaH = newAdhesiveEnergyContribution + newVolumeEnergyContribution - oldAdhesiveEnergyContribution - oldVolumeEnergyContribution

	return deltaH	
	

	

def update_cellInfo(cellInfo,cellIDs):
	''' Update the cell volumes and information 
	- For proliferating cells, the target volume increases to twice the initial volume
	- For Quiescent and Necrotic cells, the target volume is the current volume 
	- For Medium cells, there is no target volume.
	'''	
	
	for key in cellInfo.keys():
		cellInfo[key]['vol'] = len(cellIDs[cellIDs == key])		# update the vol. of the current cellID
	
		if(cellInfo[key]['type'] == 'P'): 						# update the vol energy constant & target volume
			cellInfo[key]['targetVol'] = cellInfo[key]['targetVol'] + 2
			cellInfo[key]['gamma'] = gammaP
		if(cellInfo[key]['type'] == 'N'): 
			cellInfo[key]['targetVol'] = cellInfo[key]['vol']
			cellInfo[key]['gamma'] = gammaN
		if(cellInfo[key]['type'] == 'Q'): 
			cellInfo[key]['targetVol'] = cellInfo[key]['vol']
			cellInfo[key]['gamma'] = gammaQ
		if(cellInfo[key]['type'] == 'M'): 
			cellInfo[key]['targetVol'] = 0
			cellInfo[key]['gamma'] = gammaM
	
	return cellInfo


	

def perform_montecarlo(cellInfo,cellIDs):
	''' Perform Monte Carlo Iterations by just calculating deltaH'''	
	
	iter = 0
	H0 = calculate_hamiltonian(cellInfo,cellIDs)

	listOfKeys = np.unique(cellIDs) 		# list of unique cellIDs

	nIter = len(cellIDs[cellIDs != 0])		# no. of iterations is proportional to no. of non-zero entries

	while iter <= int(300*nIter):            
		
		#
		# Pick a position on the lattice and get its cell ID
		#
		
		posX,posY = random.randint(0,nSites-1), random.randint(0,nSites-1)
		oldCellID = cellIDs[posX,posY]
		
		#
		# Make sure it is along a boundary, otherwise continue to next iteration
		#
		
		neighbourCells = cellIDs[posX-1:posX+2,posY-1:posY+2]
		nSame = len(neighbourCells[neighbourCells == cellIDs[posX,posY]])
		if(nSame == np.prod(np.shape(neighbourCells))): continue
		
		#
		# Generate a new cell ID to switch with
		#
		
		tempList = listOfKeys[listOfKeys != oldCellID]		
		newCellID = tempList[random.randint(0,len(tempList)-1)]
				
		#
		# Assign new cell ID to the current location
		#
		
		cellIDs[posX,posY] = newCellID
		iter +=1
		
		#
		# Calculate the change in Hamiltonian with this new configuration
		#
		
		dH = calculate_deltaH(cellInfo,cellIDs,posX,posY,oldCellID,newCellID)
		
		#
		# Decide which configuration to keep based on deltaH
		#
		
		if(dH < 0): 					# If Hamiltonian is lowered, keep the current state
			H0 = H0 + dH			
		else: 							# Otherwise, determine the probability of change
			probability = math.exp(-1*(dH)/kbT) 
			randProb = random.random() 				
			if(randProb<=probability): 
				H0 = H0 + dH 
			else: 
				cellIDs[posX,posY] = oldCellID 	# If a change does not occur, just use the old cell ID				
			
	return cellInfo, cellIDs, H0	



	
def perform_diffusion(cellInfo, cellIDs, oxygen):
	''' Perform pseudo-diffusion for oxygen only '''		
	#
	# Get the oxygen consumption rate at each site based on the type of cell
	#
	
	oxygenConsumptionRate = get_oxygen_consumption(cellInfo,cellIDs)
	
	#
	# Keep reducing the level of oxygen based on the no. of iterations
	#	
		
	oxygen =  oxygen - dt*nIterDiff*oxygenConsumptionRate
	
	return oxygen	
	
	

	
def update_cell_type(cellInfo,cellIDs,oxygen):
	''' Update the cell types at each lattice site based on the concentration of factors. 
	The routine is run right AFTER the DIFFUSION of factors is conducted. '''
		
	for key in cellInfo.keys():		
		
		cellOxygen = np.average(oxygen[cellIDs == key]) # Average oxygen in the cell
		
		if(cellInfo[key]['type'] == 'P'): 				# If concentration is below threshold, 'demote' from P to Q
			if(cellOxygen < thresholdP2Q):
				cellInfo[key]['type'] = 'Q'				# Only the cell type changes, ID remains the same
						
		if(cellInfo[key]['type'] == 'Q') :				# for Q to N , we need to 'move' all the cellIDs to -1 
			if(cellOxygen < thresholdQ2N):				
				cellIDs[cellIDs == key] = -1
				cellInfo[key]['targetVol'] = 0

	return cellInfo, cellIDs

	


def split_cells(cellInfo,cellIDs):
	''' Split a cell into 2 if it is a P type and has hit its target volume'''	
	
	cellInfoNew = cellInfo.copy()
	nCells = len(cellInfo) - 2
	for key in cellInfo.keys():
		if(cellInfo[key]['type']  == 'P') and (len(cellIDs[cellIDs == key]) >= 8):
		
			#
			# Split cell into two IDs
			#
			
			latticeSites = np.where(cellIDs == key)
			splitLatticeSites = [latticeSites[0][:4], latticeSites[1][:4]] 	# Choose one half of the lattice sites
			cellIDs[splitLatticeSites] = nCells + 1 						# Switch to new cell IDs 
			
			#
			# Create new dictionary item for the new cell
			#
			
			cellInfoNew[nCells+1] = cellInfoNew[key].copy()				# Update the information in the cellInfo dictionary
			cellInfoNew[nCells+1]['vol'] = 4
			cellInfoNew[nCells+1]['targetVol'] = 4
			
			#
			# Update the old cell ID to reflect the halved volume
			#
			
			cellInfoNew[key]['vol'] = 4
			cellInfoNew[key]['targetVol'] = 4
			nCells += 1

	return cellInfoNew, cellIDs

	
	

def get_oxygen_consumption(cellInfo,cellIDs):
	''' This function determines the oxygen consumption rate at each lattice site based on the cell type. '''
	
	oxygenConsumption = np.zeros((nSites,nSites))
	
	for key in cellInfo.keys():
		if(cellInfo[key]['type'] == 'N'):	oxygenConsumption[cellIDs == key] = oxygenN
			
		if(cellInfo[key]['type'] == 'Q'):	oxygenConsumption[cellIDs == key] = oxygenQ
	
		if(cellInfo[key]['type'] == 'P'):	oxygenConsumption[cellIDs == key] = oxygenP
				
	return oxygenConsumption
	

	

def print_to_file(cellInfo,cellIDs,i):
	''' Print cellIDs to file and plot the cell growth. '''	
	
	displayIDs = cellIDs.copy()		# to show only type of cell and not all cellIDs
	
	for key in cellInfo.keys():

		if(cellInfo[key]['type'] == 'P'):		
			displayIDs[cellIDs == key] = 1
	
		if(cellInfo[key]['type'] == 'Q'):		
			displayIDs[cellIDs == key] = 2

	#
	# Plot the tumour growth
	#
	
	colors = 'black blue red yellow'.split()
	cmap = matplotlib.colors.ListedColormap(colors, name='colors', N=None)
	plt.imshow(displayIDs, cmap=cmap, vmin=-1, vmax=2)
	plt.savefig(dirName + "/display_"+ str(initOxygen) + "_" + str(i) + ".png")
	

	#
	# Save cellIDs and displayIDs
	#
	
	filename = dirName + "/display_"+ str(initOxygen) + "_" + str(i) + ".txt"
	filename1 = dirName + "/cellIDs_"+ str(initOxygen) + "_" + str(i) + ".txt"
	np.savetxt(filename,displayIDs,delimiter=' ',newline='\n ', fmt = '%d')
	np.savetxt(filename1,cellIDs,delimiter=' ',newline='\n ', fmt = '%d')



	
def get_prognosis(cellIDs):
	''' Display the prognosis for the patient '''
	
	#
	# Print the ratio of tumour sites (all types) to the total number of sites
	#
	
	print('Size of tumour is : ', len(cellIDs[cellIDs != 0])/(nSites)**2)
	if((len(cellIDs[cellIDs != 0])/(nSites)**2) > 0.0225): 
		print('Your tumour is Metastatic :-( ')
	else:
		print('Your tumour is Benign :-)' )



		
def generate_gif():
	''' Code to generate the gifs'''
	# Load each file into a list
	frames = []
	for i in range(0,21):
		frames.append(imageio.imread(dirName + "/display_" + str(initOxygen) + "_" + str(i) + ".png"))
	
	# Save them as frames into a gif 
	exportName = dirName + "/output.gif"
	kargs = { 'duration': 1 }
	imageio.mimsave(exportName, frames, 'GIF', **kargs)			
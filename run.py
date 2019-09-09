import numpy as np
from implementations import *
from config import *
import os

if __name__ == '__main__':
	''' Main function for the pseudo oxygen diffusion '''	
		
	#
	# ********************* PROBLEM SETUP ******************************
	#
	
	
	# Initialize matrices for cell Latice sites
	cellIDs = np.zeros((nSites,nSites)) 		# Cell ID	
	oxygen = np.ones((nSites,nSites))*initOxygen			# Oxygen concentration

		
	# Initialize a single proliferating cell at the centre
	cellCentre = int((nSites-1)/2)
	cellIDs[cellCentre-1:cellCentre+1,cellCentre-1:cellCentre+1] = 1
		

	# Create a dictionary that stores cell information 
	cellInfo = {}
	
	# Initialize with the medium, necrotic and a single proliferating cell
	cellInfo[1] = {'type' : 'P', 'vol' : 4, 'targetVol' : 6, 'gamma' : 10}
	cellInfo[0] = {'type' : 'M', 'vol' : nSites*nSites - 4, 'targetVol' : 0, 'gamma' : 0} 
	cellInfo[-1] = {'type' : 'N', 'vol' : 0, 'targetVol' : 0, 'gamma' : 50} # Necrotic cells will all be added to the ID -1

	
	#
	# ************************* SIMULATION LOOP *******************************
	#
	
	H0 = calculate_hamiltonian(cellInfo,cellIDs) # base energy

	for i in range(21):
		print('iteration = ', i)
		
		# print current cellIDs to file
		print_to_file(cellInfo,cellIDs,i)
		
		
		# update cell volumes and constants 
		cellInfo = update_cellInfo(cellInfo,cellIDs)

		
		# montecarlo steps to ensure minimum energy physical state
		cellInfo,cellIDs,H0 = perform_montecarlo(cellInfo,cellIDs)
		
		
		# consumption of oxygen by proliferating and quiescent cells
		oxygen = perform_diffusion(cellInfo,cellIDs, oxygen)
		
		
		# update cell IDs and type based on oxygen concentration
		cellInfo, cellIDs = update_cell_type(cellInfo,cellIDs,oxygen)
		
		
		# split proliferating cells into 2 if required
		cellInfo, cellIDs = split_cells(cellInfo,cellIDs)
		
	# get the final prognosis
	get_prognosis(cellIDs)	
	
	# generate the GIF animation
	generate_gif()
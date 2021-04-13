###### Tiny little script to convert the output of the Spacal Summary simulations from root file to numpy.
###### USAGE:
### 		python3 rootToNumpy.py fileToConvert.root
###### the output will be a numpy file:
###		fileToConvert.root_numpy.npy
###### Dec2019 L.M.

import os
import sys
import ROOT
import numpy as np


datafilename 	= sys.argv[1]			### .root simulation output file.
dataTTreename  	= "summary"             ### Name of the ttree containing the data.





### Load data in a dataframe
print("Converting:", datafilename)
print ("... creating the dataframe...")
df = ROOT.RDataFrame(dataTTreename, datafilename)
columns = df.GetColumnNames()

### Convert dataframe to dictionary
print ("... converting it to numpy array...")
data_dict = df.AsNumpy(columns=df.GetColumnNames())
#print(data_dict)

### Dictionary to list of arrays
list = [array for ind, array in data_dict.items()]

### List of arrays to array of arrays and transpose it
matrix = np.stack(list)
matrix = matrix.T
#print("### stack: ", matrix)
#print("Stack format: ", np.shape(matrix))


### Save data
np.save(datafilename+"_numpy.npy", matrix)


print ("... done.")

# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 22:07:43 2021

@author:
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This script extracts the information from the shape model created
    from gias2.
    
    Use gias2-py27 environment for gias2 functionality.
    
    TODO: add more notes...
    
"""

# %% Import packages

import os
# import gias2
import numpy as np
from gias2.learning import PCA
from gias2.fieldwork.field import geometric_field

# %% Set-up

#Navigate to shape model folder
os.chdir('..\\ShapeModel\\shape_model')

# %% Load in shape model

#Load the pca file
pc = PCA.loadPrincipalComponents('tibia.pc.npz')

#Print the % variance for each PC
#Get the components
componentVar = pc.getNormSpectrum()
#Set the number to print
nPrintPC = 10
for ii in range(nPrintPC):
   print('PC%d: %4.2f%%' % (ii + 1, componentVar[ii] * 100))
   
#Print cumulative variance
currCumulative = 0
for ii in range(nPrintPC):
   print('PC%d: %4.2f%%' % (ii + 1,
                            (currCumulative+componentVar[ii]) * 100))
   #Update cumulative variance
   currCumulative += componentVar[ii]
   
#Get the mean points and reshape to n x 3
meanPoints = pc.getMean()
meanPoints = meanPoints.reshape(int(len(meanPoints)/3),3)

#Export reconstructed points to file
np.savetxt('pointCloud_mean.csv', meanPoints, delimiter = ',')

#Reconstruct the first 10 PCs and get the points
for cc in range(10):

    #Set PC to reconstruct points
    reconPC = [cc, ]
    
    #Reconstruct points +/- 3SD's
    sdPlusPoints = pc.reconstruct(pc.getWeightsBySD(reconPC, [+3.0, ]), reconPC)
    sdMinusPoints = pc.reconstruct(pc.getWeightsBySD(reconPC, [-3.0, ]), reconPC)
    
    #Reshape to n x 3
    sdPlusPoints = sdPlusPoints.reshape(int(len(sdPlusPoints)/3),3)
    sdMinusPoints = sdMinusPoints.reshape(int(len(sdMinusPoints)/3),3)

    #Export reconstructed points to file
    np.savetxt('pointCloud_pc'+str(cc+1)+'_p3.csv', sdPlusPoints, delimiter = ',')
    np.savetxt('pointCloud_pc'+str(cc+1)+'_m3.csv', sdMinusPoints, delimiter = ',')


# %%

# %% Old below - raw data variant

# %% Import packages

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# %% Set-up

#Navigate to shape model folder
os.chdir('..\\ShapeModel\\shape_model')

# %% Load shape model data

modelMean = np.load('mean.npy', allow_pickle = True)
modelMean = modelMean.reshape(int(len(modelMean)/3),3)
modelModes = np.load('modes.npy', allow_pickle = True)
modelSD = np.load('SD.npy', allow_pickle = True)
modelSizes = np.load('sizes.npy', allow_pickle = True)
modelWeights = np.load('weights.npy', allow_pickle = True)
modelProjWeights = np.load('projectedWeights.npy', allow_pickle = True)



# %% test code...

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(modelMean[:,0],
           modelMean[:,1], 
           modelMean[:,2])
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
import gias2

# %% Set-up

#Navigate to shape model folder
os.chdir('..\\ShapeModel\\shape_model')

# %% Load in shape model







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
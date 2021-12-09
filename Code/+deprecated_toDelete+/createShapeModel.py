# -*- coding: utf-8 -*-
"""

@authors:
    Meghan Keast
    Centre for Sport Research
    Deakin University
    mfkeast@deakin.edu.au
    
    Aaron Fox
    Centre for Sport Research
    Deakin University
    aaron.f@deakin.edu.au
    
    This script is a fairly simple storage file that holds to commands to run
    the statistical shape model processing workflow in gias2. The commands encased
    within the 'os.system()' commands can actually be run from the command line,
    but we included them within this script for ease of use within Python IDE's
    such as Spyder.
    
"""

# %% Import packages

import os
from gias2.learning import PCA

# %% Set-up

#Set home directory
homeDir = os.getcwd()

# %% Generate tibia shape model

#Navigate to tibia shape model folder
os.chdir('..\\ShapeModels\\tibia')

#Run mesh fitting step
#This step non-rigidly registers the surface meshes to a 'master' mesh (which
#is the first in the .txt list), thereby representing all of the segmentations
#by correspondent meshes. The gias-rbfreg uses radial basis functions to register
#non-correspondent point clouds.
os.system('gias-rbfreg -b rbfreg_list.txt -d fitted_meshes --outext .stl')

#Run rigid alignment step
#This step rigidly aligns the meshes to remove rotational and translational variations.
#Here is where there is the option to not scale or scale the meshes to remove size
#variation too - hence why two variants are provided (i.e. corr_r vs. corr_rs argument)
os.system('gias-rigidreg corr_r -b rigidreg_list.txt -d aligned_meshes --outext .stl')
os.system('gias-rigidreg corr_rs -b rigidreg_list.txt -d aligned_meshes_scaled --outext .stl')

#Run PCA to generate shape model
#The last step is to run principal component analysis on the surface meshes to generate
#the shape model. Here we also generate meshes of the mean model and the principal
#components at iterative standard deviations to the mean. As above, we can generate
#shape models on unscaled and scaled meshes. Here we also provide the function how many
#principal components to create, and this is the number of meshes - 1 to create a full model.
#We however only reconstruct meshes for the first 7 components, as this accounts for 97% of
#variance in the unscaled model.
os.system('gias-trainpcashapemodel pca_list.txt -n 29 -o shape_model/tibia')
os.system('gias-trainpcashapemodel pca_list_scaled.txt -n 29 -o shape_model_scaled/tibia')

#Export shape model for use in Matlab
#This is to ensure the principal component model can be easily loaded and used
#in Matlab in subsequent scripts.
tibiaPC = PCA.loadPrincipalComponents('shape_model\\tibia.pc.npz')
tibiaPC.savemat('shape_model\\tibia.pc.mat')
tibiaPC_scaled = PCA.loadPrincipalComponents('shape_model_scaled\\tibia.pc.npz')
tibiaPC_scaled.savemat('shape_model_scaled\\tibia.pc.mat')

# %% Generate tibia-fibula shape model

#Note that all of the comments in the above section relating to the tibia model, apply
#equally in this section for creating the tibia-fibula model.

#Navigate to tibia shape model folder
os.chdir('..\\tibia-fibula')

#Run mesh fitting step
os.system('gias-rbfreg -b rbfreg_list.txt -d fitted_meshes --outext .stl')

#Run rigid alignment step
os.system('gias-rigidreg corr_r -b rigidreg_list.txt -d aligned_meshes --outext .stl')
os.system('gias-rigidreg corr_rs -b rigidreg_list.txt -d aligned_meshes_scaled --outext .stl')

#Run PCA to generate shape model
os.system('gias-trainpcashapemodel pca_list.txt -n 28 -o shape_model/tibia-fibula')
os.system('gias-trainpcashapemodel pca_list_scaled.txt -n 28 -o shape_model_scaled/tibia-fibula')

#Export shape model for use in Matlab
tibFibPC = PCA.loadPrincipalComponents('shape_model\\tibia-fibula.pc.npz')
tibFibPC.savemat('shape_model\\tibia-fibula.pc.mat')
tibFibPC_scaled = PCA.loadPrincipalComponents('shape_model_scaled\\tibia-fibula.pc.npz')
tibFibPC_scaled.savemat('shape_model_scaled\\tibia-fibula.pc.mat')

# %% Generate trabecular shape model

#Note that all of the comments in the above section relating to the tibia model, apply
#equally in this section for creating the trabecular model.

#Navigate to tibia shape model folder
os.chdir('..\\trabecular')

#Run mesh fitting step
os.system('gias-rbfreg -b rbfreg_list.txt -d fitted_meshes --outext .stl')

#Run rigid alignment step
os.system('gias-rigidreg corr_r -b rigidreg_list.txt -d aligned_meshes --outext .stl')
os.system('gias-rigidreg corr_rs -b rigidreg_list.txt -d aligned_meshes_scaled --outext .stl')

#Run PCA to generate shape model
os.system('gias-trainpcashapemodel pca_list.txt -n 29 -o shape_model/trabecular')
os.system('gias-trainpcashapemodel pca_list_scaled.txt -n 29 -o shape_model_scaled/trabecular')

#Export shape model for use in Matlab
trabecularPC = PCA.loadPrincipalComponents('shape_model\\trabecular.pc.npz')
trabecularPC.savemat('shape_model\\trabecular.pc.mat')
trabecularPC_scaled = PCA.loadPrincipalComponents('shape_model_scaled\\trabecular.pc.npz')
trabecularPC_scaled.savemat('shape_model_scaled\\trabecular.pc.mat')

# %% Finish up

#Return to home directory
os.chdir(homeDir)

# %%% ----- End of createShapeModel.py -----
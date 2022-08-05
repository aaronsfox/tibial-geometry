# Shape Models

*We recommend downloading a Markdown reader like [Typora](https://typora.io/) to view the readme files within this dataset.*

This folder contains the shape model data for the three models included with our paper and dataset, those being — (1) tibia; (2) tibia and fibula; and (3) tibial trabecular. Below outlines the details of the files and folders included with each shape model.

## {model}ShapeModel.mat

This .mat file contains all of the data associated with the shape model (e.g. mean values, loadings etc.) and can be loaded into MATLAB to be used and manipulated as desired. The example cases presented with this repository provide multiple demonstrations of loading and using the shape model.

```matlab
%Example MATLAB code to load shape model data
load('tibiaShapeModel.mat')
```

## {model}ShapeModel_mean.stl

The .stl file within each shape model folder represents the mean points of the model as a surface file. Note that visually this will look the same as the mean points included in the shape model .mat file — however due to the way in which the .stl was exported, the ordering of the points is changed.

## figures

The figures folder contains two sub-folders:

- **errorSummary:** The figure included in this folder is a summary of mean and peak position errors, and the Jaccard index for reconstructions of the data included in the shape model using the reduced principal component set. These figures for each model are included in our preprint/manuscript associated with these data.
- **reconstructionErrors:** The figures included in this folder are comparisons of surfaces reconstructed using the reduced principal component set of the model to the original surfaces. A heatmap approach is used to distinguish the magnitude of point-to-point difference between the reconstructed and original surface (shown in light grey).

## results

The results folder contains a .csv file that tabulates the reconstruction error data for the original surfaces using the reduced principal component set. These data are presented in the 'errorSummary' folder in figure format (see above in **figures**). Animations of the effect of each principal component on the mean model are also included here in .gif format.

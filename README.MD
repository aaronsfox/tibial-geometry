# Statistical Shape Model of the Tibia-Fibula

*We recommend downloading a Markdown reader like [Typora](https://typora.io/) to view the readme files within this dataset.*

This README file aims to provide an overview of key considerations and elements included within this dataset.

## Dependencies

A series of MATLAB dependencies are required to run the code associated with this dataset. Our code has been tested on MATLAB version 2019b. Some specific code dependencies are included in the **supplementary** directories in various folders, and these are added to the MATLAB path at the beginning of the functions where they are used. Other external major dependencies that aren't included directly with the dataset are:

- [GIBBON](https://www.gibboncode.org/): used to prepare and analyse the 3D surface data
- [geom3d](https://au.mathworks.com/matlabcentral/fileexchange/24484-geom3d): used to simplify the process of applying 3D geometric calculations
- [export_fig](https://au.mathworks.com/matlabcentral/fileexchange/23629-export_fig): used to export high resolution figures from MATLAB

## Segmentation

Lower limb computed tomography (CT) scans from the right tibia and fibula of 30 cadavers (male n = 20, female n =10) were selected for use in this dataset. CT scans were provided by the [New Mexico Decedent Image Database](https://nmdid.unm.edu/) (NMDID). The NMDID dataset includes whole body CT scans of over 15,000 New Mexicans who died between 2010-2017. CT scan Slice thickness is 1 mm with 0.5 mm overlap. The images included for use in creation of the dataset were from individuals between 18 to 40 years who participated in impact-based physical activities (i.e. team sports, dancing, recreational running and walking). Images were not selected for inclusion if there was any noticeable damage to the right tibia-fibula or the inclusion of medical devices (plates, rods, or screws). Images were segmented using mimics innovation suite.

Tibias were segmented into two surfaces representing the outer boundaries of the trabecular and cortical bone, while fibulas where segmented as one surface representing the outer shape of the entire bone. Cortical bone was classified as all visible hard bone, while trabecular was all soft bone and the medulla cavity. The CT files and segmentation files will not be provided as a part of this dataset. STL files of the segmented images are available in the [segmentation folder](Segmentation).

## Code

The following points outline the functions and processes involved of the various scripts included with this dataset. The order of these dot points also corresponds to the general order in which these scripts were each run to generate and analyse the shape models.

- **createShapeModel.m:** This script contains the required commands to (i) remesh the surfaces to an even number of points; (ii) non-rigidly and rigidly register the surfaces to a target mesh; and (iii) create the statistical shape models. The process of registering the surfaces is time consuming (~1-2 hours), and hence an option is provided to do this or skip and load the existing data (see **Data** section below). This code also saves the shape model data to be used at a later stage.

The **supplementary** folder includes a series of external and custom functions used within the main script for various steps. Where using these external code sources for your own work, please cite the original authors of these works.

## Shape Models

The [shape models](ShapeModels) folder contains all data associated with the statistical shape models from our work. See the [README](ShapeModels/README.MD) file associated with this folder for further details on what is included here. 

## Data

The [data](data) folder contains one file — *registeredSurfaceDataPoints.mat* — a file that contains the three-dimensional point data of the bone surfaces used to create the statistical shape models. The process of registering the points from the original surfaces is time-consuming, hence we have saved and included this file here for those who wish to re-create the shape models but avoid this computationally expensive step. This can be loaded into MATLAB using the basic command:

```matlab
load('registerdSurfaceDataPoints.mat')
```

## Examples

This dataset and associated statistical shape models can be used in several ways to assist with skeletal focused research of the tibia-fibula. We do not have the scope to highlight each and every potential application, however have provided a series of example cases of where and how the shape models may be used. Our hope is that these examples can be directly used, or assist in guiding other uses.

- **[Case 1: Generating Surface Samples](Examples/case1)** — this example case demonstrates how to use the shape model data to reconstruct a randomly sampled 'population' of surfaces.
- **[Case 2: Predicting and Generating Trabecular Volumes](Examples/case2)** — this example case demonstrates how to combine the tibia and trabecular shape models to predict and generate the trabecular volume from a tibial surface.
- **[Case 3: Generating Tibia-Fibula Surfaces from Landmarks](Examples/case3)** — this example case demonstrates how to use the tibia-fibula shape model to estimate and reconstruct surfaces from palpable landmarks on the tibia and fibula.

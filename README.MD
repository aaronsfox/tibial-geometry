# Statistical Shape Model of the Tibia-Fibula

*We recommend downloading a Markdown reader like [Typora](https://typora.io/) to view the readme files within this dataset.*



## Dependencies

A series of MATLAB dependencies are required to run the code associated with this dataset. Our code has been tested on MATLAB version 2019b. Certain code dependencies are included in the 'Code' or 'Supplementary' directories in various folders, and these are typically added to the MATLAB path at the beginning of the code. Other major dependencies include:

- [GIBBON](https://www.gibboncode.org/): used to prepare and analyse the 3D surface data
- [geom3d](https://au.mathworks.com/matlabcentral/fileexchange/24484-geom3d): used to simplify the process of applying 3D geometric calculations
- [export_fig](https://au.mathworks.com/matlabcentral/fileexchange/23629-export_fig): used to export high resolution figures from MATLAB



## Segmentation

Lower limb computed tomography (CT) scans from the right tibia and fibula of 30 cadavers (male n = 20, female n =10) were selected for use in this dataset. CT scans were provided by the [New Mexico Decedent Image Database](https://nmdid.unm.edu/) (NMDID). The NMDID dataset includes whole body CT scans of over 15,000 New Mexicans who died between 2010-2017. CT scan Slice thickness is 1 mm with 0.5 mm overlap. The images included for use in creation of the dataset were from individuals between 18 to 40 years who participated in impact-based physical activities (i.e. team sports, dancing, recreational running and walking). Images were not selected for inclusion if there was any noticeable damage to the right tibia-fibula or the inclusion of medical devices (plates, rods, or screws). Images were segmented using mimics innovation suite.
Tibias were segmented into two surfaces representing the outer boundaries of the trabecular and cortical bone, while fibulas where segmented as one surface representing the outer shape of the entire bone. Cortical bone was classified as all visible hard bone, while trabecular was all soft bone and the medulla cavity. 

****The CT files and segmentation files will not be provided as apart of this dataset. STL files of the segmented images will be made avalible.

## Code

The following points outline the functions and processes involved of the various scripts included with this dataset. The order of these dot points also corresponds to the general order in which these scripts were each run to generate and analyse the shape models.

- **createShapeModel.m:** This script contains the required commands to (i) remesh the surfaces to an even number of points; (ii) non-rigidly and rigidly register the surfaces to a target mesh; and (iii) create the statistical shape models. The process of registering the surfaces is time consuming (~1-2 hours), and hence an option is provided to do this or skip and load the existing data. This code also saves the shape model data to be used at a later stage.



## Examples

***TODO: add example list and links to documentation***

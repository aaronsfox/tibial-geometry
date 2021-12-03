# Statistical Shape Model of the Tibia



NOTES:

- Shape model developed using gias2 in Python 2.7 environment
- GIBBON toolbox used to prepare and analyse data
- Any Python 3 used?
- Need geom3d toolbox
- case-134065 has incomplete fibula due to damage
- Consider use of GIBBON anim8 to demonstrate shape features





STEPS:

- Original segmentations placed into starting 'ShapeModels' folder by ***setUpShapeModel.m*** code



## Code

The following points outline the functions and processes involved of the various scripts included with this dataset. The order of these dot points also corresponds to the order in which these scripts were each run to generate and analyse the shape models. Both MATLAB and Python code were used in this project — mainly due to the usefulness of both MATLAB (i.e. GIBBON, geom3d) and Python (i.e. gias2) packages for processing these type of data.

- **setUpShapeModel.m:** This script was run on the raw segmented tibia and fibula surfaces, along with the digitised landmarks to — (i) smooth and remesh the surfaces; and (ii) align the surfaces to an appropriate tibial coordinate system. This step was a precursory step to the statistical shape model and served to speed-up this process by reducing the total number of vertices on the surfaces and roughly aligning them to the same point in 3D space.
- **createShapeModel.py:** This script contains the required commands to pass through *gias2* to generate the statistical shape models for the tibia and tibia-fibula approaches. The commands within this can be passed straight into the command line, but are housed in this .py script as well for easier use in Python IDEs, such as Spyder. The last step in this script is to export the created shape models to a .mat file for use in later MATLAB scripts.


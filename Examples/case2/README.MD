# Case 2: Generate Trabecular from Surface Mesh

*We recommend downloading a Markdown reader like [Typora](https://typora.io/) to view the readme files within this dataset.*

## Background

This example case uses both the tibia and trabecular shape model created from this dataset to form a predictive model of the trabecular surface from the tibia surface. The theoretical basis for this is that the shape of the tibia can be effective in predicting the shape of the internal trabecular surface. The practical basis for this is that there are certain scenarios where generating actual models of the trabecular surface are difficult or improbable. For example:

- Segmenting the trabecular surface from CT scans is possible, but is more time-consuming and requires more manual correction than the overall external surface (i.e. the cortical bone)
- MRI scans that aren't optimised for detecting the trabecular surface do not provide great information for segmenting the trabecular surface
- Bone surface reconstruction methods from incomplete scans (e.g. using something like the [MAP-Client](https://simtk.org/projects/map)) only provide an estimate of the external bone surface — not the internal trabecular surface

There are various scenarios where an adequate model of the cortical and trabecular surfaces is required. Finite element simulations of tibial stress during running is one example, where including the different material properties of these different bone areas is a part of the simulation process (e.g. Edwards et al. 2010; Haider et al. 2020).





NOTES:

- Approach here is to demonstrate how to generate a volumetric trabecular mesh of the tibia from a tibia or tibia-fibula surface mesh and the relevant surface shape models
- 





DEPENDENCIES:

- GIBBON toolbox
- geom3d
- nonrigidICP from https://au.mathworks.com/matlabcentral/fileexchange/41396-nonrigidicp (with some edits to suppress the visualisations)
- intriangulation from https://au.mathworks.com/matlabcentral/fileexchange/43381-intriangulation-vertices-faces-testp-heavytest



## Code

- **TODO.m:** ...





## References

Edwards et al. (2010). Effects of running speed on a probabilistic stress fracture model. *Clin Biomech*, 25: 372-377.

Haider et al. (2020). Subject-specific finite element models of the tibia with realistic boundary conditions predict bending deformations consistent with in vivo measurement. *J Biomech Eng*, 142: 021010. 

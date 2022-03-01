# Case 3: Generating Tibia-Fibula Surfaces from a Minimal Landmark Set

*We recommend downloading a Markdown reader like [Typora](https://typora.io/) to view the readme files within this dataset.*

## Background

This example case uses landmarks extracted from a static motion capture trial to form a predictive model of the tibia-fibula surface using the shape model as a basis for the prediction. The theoretical basis for this is that the relative positioning of landmarks on the tibia-fibula surface can be used to effectively reconstruct the surface (e.g. Bruce et al. 2021). The practical basis for this is that there are many scenarios in biomechanical laboratories where motion capture data of this nature can be collected — however, collecting the additional imaging (e.g. CT scans) are unavailable due to prohibitive access or costs. Providing a method to develop surface models of the tibia-fibula complex from a more common set of data can increase access to research using structural models.





DEPENDENCIES:

- GIBBON toolbox
- geom3d
- nonrigidICP from https://au.mathworks.com/matlabcentral/fileexchange/41396-nonrigidicp (with some edits to suppress the visualisations)
- BTK - https://code.google.com/archive/p/b-tk/downloads



## Code

- **TODO.m:** ...





## References

Bruce et al. (2021). A statistical shape model of the tibia-fibula complex: sexual dimorphism and effects of age on reconstruction accuracy from anatomical landmarks. *Comput Meth Biomech Biomed Eng*, doi: 10.1080/10255842.2021.1985111.

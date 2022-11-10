# Case 2: Generate Trabecular from Surface Mesh

*We recommend downloading a Markdown reader like [Typora](https://typora.io/) to view the readme files within this dataset.*

## Background

This example case uses the combined tibia and trabecular shape model created from this dataset to make predictions of the trabecular surface from the outer tibial surface data. The theoretical basis for this is that the shape of the outer cortical tibia can be effective in predicting the shape of the internal trabecular surface. The practical basis for this is that there are certain scenarios where generating actual models of the trabecular surface are difficult or improbable. For example:

- Segmenting the trabecular surface from CT scans is possible, but is more time-consuming and requires more manual correction than the overall external surface (i.e. the cortical bone)
- MRI scans that aren't optimised for detecting the trabecular surface do not provide great information for segmenting the trabecular surface
- Bone surface reconstruction methods from incomplete scans (e.g. using something like the [MAP-Client](https://simtk.org/projects/map)) only provide an estimate of the external bone surface — not the internal trabecular surface

There are various scenarios where an adequate model of the cortical and trabecular surfaces is required. Finite element simulations of tibial stress during running is one example, where including the different material properties of these different bone areas is a part of the simulation process (e.g. Edwards et al. 2010; Haider et al. 2020).

## Code

- **generateTrabecularFromSurface.m:** This is the main script associated with this example case, and runs through the processes of: (i) comparing predicted reconstructions with the original surfaces as a means of evaluating the accuracy of this process; and (iii) predicting the trabecular shape of a new surface from Nolte et al. (2016).

The **supplementary** folder includes a series of external and custom functions used within the main script for various steps. Where using these external code sources for your own work, please cite the original authors of these works.

## References

Edwards et al. (2010). Effects of running speed on a probabilistic stress fracture model. *Clin Biomech*, 25: 372-377.

Haider et al. (2020). Subject-specific finite element models of the tibia with realistic boundary conditions predict bending deformations consistent with in vivo measurement. *J Biomech Eng*, 142: 021010. 

Nolte et al. (2016). Non-linear scaling of a musculoskeletal model of the lower limb using statistical shape models. *J Biomech*, 49: 3576-3581.

# Case 3: Generating Tibia-Fibula Surfaces from a Minimal Landmark Set

*We recommend downloading a Markdown reader like [Typora](https://typora.io/) to view the readme files within this dataset.*

## Background

This example case uses an optimisation process combined with the tibia-fibula shape model to reconstruct surfaces from Nolte et al. (2016) from a set of landmarks digitised on the surfaces. The theoretical basis for this is that the relative positioning of landmarks on the tibia-fibula can be used to effectively reconstruct the surface (e.g. Bruce et al. 2021). The practical basis for this is that there are many scenarios in biomechanical laboratories where landmark data can be readily collected (e.g. during static motion capture trials) — however, collecting additional imaging (e.g. CT scans) are unavailable due to prohibitive access or costs. Providing a method to develop surface models of the tibia-fibula from a more common set of data can increase access to research using three-dimensional models.

## Code

- **generateSurfaceFromLandmarks.m:** This is the main script associated with this example case, and runs an optimisation process that fits the tibia-fibula statistical shape model to the digitised landmarks on the sample surfaces from Nolte et al. (2016). The reconstructions are then evaluated for error against the original surfaces.

The **supplementary** folder includes a series of external and custom functions used within the main script for various steps. Where using these external code sources for your own work, please cite the original authors of these works.

## References

Bruce et al. (2021). A statistical shape model of the tibia-fibula complex: sexual dimorphism and effects of age on reconstruction accuracy from anatomical landmarks. *Comput Meth Biomech Biomed Eng*, doi: 10.1080/10255842.2021.1985111.

Nolte et al. (2016). Non-linear scaling of a musculoskeletal model of the lower limb using statistical shape models. *J Biomech*, 49: 3576-3581.

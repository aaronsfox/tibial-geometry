# Case 1: Generating Surface Samples

*We recommend downloading a Markdown reader like [Typora](https://typora.io/) to view the readme files within this dataset.*

## Background

This example case uses uses the SSM of the tibia-fibula to create a sample set of tibia-fibula surfaces representing standardized variation across all and/or a selection of the model components. The practical basis for this is that there are many scenarios where collecting imaging (e.g. CT scans) to develop a large sample of bone models is unattainable due to prohibitive access or costs. Providing a method to develop surface models of the tibia-fibula increases access to research using three-dimensional models. This code can be used to produce a representative sample of surfaces or a ‘simulated population.’ Further, this example provides the ability to manipulate specific components of the shape models allowing samples to be generated with isolated geometric variations of interest.

## Code

- **simulatedPopulations.m:** This is the main script associated with this example case. It allows users to generate a select number of surfaces from a pre-loaded shape model. While the example provided uses the tibia-fibula shape model, this process could be applied to any of the shape models included in our dataset. The function is structured to allow users to: (i) load the desired shape model; (ii) set the desired number of samples (i.e. n) to be generated; (iii) select the components of the model to be perturbed in creating the surface variations; (iv) set the magnitude of variation to perturb the model components by (in standard deviations); and (v) set the output options. This function takes user inputs and randomly samples values for the chosen components within the standard deviation bounds provided to reconstruct n samples. Users can also select to output a visual representation of the surfaces with heatmaps applied to highlight the areas of variation, as well as outputting the three-dimensional surfaces in STL format. 

## Note

It is important to be wary of the limitations associated with this process when using the tibia-plus-trabecular model. This shape model includes the trabecular surface inside the outer cortical surface — and without cautious manipulation of the principal components can result in projection of the trabecular surface through the cortical surface (i.e. an unrealistic scenario). When using this shape model, we therefore recommend:

1. Inspect any surfaces created for projection of the trabecular through the cortical surface
2. Consider the standard deviation bounds to place on the various principal components. It is likely that only small deviations can be used on principal components which impact the trabecular more than the cortical surface.
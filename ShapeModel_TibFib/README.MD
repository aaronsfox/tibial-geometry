***Used gias2-py37 environment***

gias2 input commands:

- gias-rbfreg -b rbfreg_list.txt -d fitted_meshes --outext .stl
- gias-rigidreg corr_r -b rigidreg_list.txt -d aligned_meshes --outext .stl
- gias-trainpcashapemodel pca_list.txt -n 3 -r 0 1 2 -o shape_model/tibia-fibula --plot_pcs 3 -v
  - NOTE: fix up number of components


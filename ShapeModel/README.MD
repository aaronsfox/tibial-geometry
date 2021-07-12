***Used gias2-py37 environment***

gias2 input commands:

- gias-rbfreg -b rbfreg_list.txt -d fitted_meshes --outext .stl
- gias-rigidreg corr_r -b rigidreg_list.txt -d aligned_meshes --outext .stl
- gias-trainpcashapemodel pca_list.txt -n 10 -r 0 1 2 3 4 -o shape_model/tibia --plot_pcs 10 -v
  - Note: 141682 removed for now from shape model given lack of landmarks and registration issues
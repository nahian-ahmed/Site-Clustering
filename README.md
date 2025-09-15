<!----------------------------------------------------------->
<!----------------------------------------------------------->
<!-- June 27, 2025                                         -->
<!----------------------------------------------------------->
<!----------------------------------------------------------->


# A comparison of clustering approaches to create sites for occupancy models from opportunistic biodiversity surveys

<!-- <p align="center">
  <img src="site_clustering_process_flow.png"  width="800">
</p> -->

## 1. Data 

Data available at <https://doi.org/10.5281/zenodo.13118082>.

`occupancy_feature_raster` contains raster file of occupancy/site features, and `checklist_data` contains eBird checklists of 31 species. 

Please place both directories in root directory before executing code.


## 2. Requirements

Code was last run using `R` version 4.4.2, `gdal` version 3.10.0, `geos` version 3.13.0, and `proj` version 9.5.0.

Required `R` packages are listed in `dependencies.txt`.

## 3. Instructions

### 3.1. Experiments on simulation data

Run the following command:

`Rscript run_simulation_experiments.R`


### 3.2. Experiments on species data
Run the following command:

`Rscript run_species_experiments.R`

Replace `species_names` in `run_species_experiments.R` to run experiments for specific species. Abbreviations of species names are in `checklist_data` (e.g, `NOFL` = Northern Flicker, `COHA` = Cooper's Hawk, etc.).
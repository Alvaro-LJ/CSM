# <img src="man/figures/Sticker.png" width="60" style = "align-items: center"> Comprehensive Spatial Methods (CSM) 
**Summary:** CSM is an R package designed as a toolbox to analyze spatially resolved tissue data. To see a demo of CSM capabilities please see related [*publication*](publicationURL).


## What you will find in this Repository (Folder description)
[Requirements and installation](#Main-features-and-requirements) 
[Main outline of CSM](#Main-outline-of-CSMs) 
[Publication and associated datasets](#Publication) 
[Citation](#Images) 

## Main features and requirements
* CSM was developed using R version > 4.3
* CSM depends on many packages to work. However not all suggested packages are required to perform specific tasks.
* CSM includes >100 functions organized in 6+1 modules.
* CSM provides several datasets to test function capabilities.
* You can install CSM in R using `devtools::install_github(Alvaro-LJ/CSM)`.

## Main outline of CSM
<img src="man/figures/CSM_OUTLINE.png" width="750" style = "align-items: center">

- MODULE 0: Data generation and formatting. Includes functions to perform the following tasks:
  - Divide large images into tiles. Useful to work with whole slide tissue samples.
  - Perform pixel thresholding and quantification. Used to work with extra-cellular biomolecules.
  - Extract color channels from RGB images. May be required to work with HE, IHC or other histochemical staining techniques.
  - Perform cell segmentation and feature extraction from images.
  - Arrange cell feature matrices into an adequate format.
  - Arrange image metadata into an adequate format.
  - Perform general quality checks and check computational resources available.
- MODULE 1: Data normalization. Includes functions to perform the following tasks:
  - Normalize cell feature expression data according to image and slide belonging.
- MODULE 2: Data thresholding. Includes functions to perform the following tasks:
  - Define positive thresholds for features of interest. Thresholded data can be used to assign cell labels.

## Publication
This folder contains the CSM scripts and datasets to replicate the results of our [*publication*](publicationURL).





## Citation
Please cite this paper in case our method or parts of it were helpful in your work.
```diff
@article{XXX,
  title={XXX},
  author={XXX},
  journal={XXX},
  year={XXX}
}
```

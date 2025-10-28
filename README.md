# <img src="man/figures/Sticker.png" width="60" style = "align-items: center"> Comprehensive Spatial Methods (CSM) 
**Summary:** CSM is an R package designed as a toolbox to analyze spatially resolved tissue data. To see a demo of CSM capabilities please see related [*publication*](publicationURL).


## What you will find in this Repository (Folder description)
+ [Requirements and installation](#Main-features-and-requirements) 
+ [Main outline of CSM](#Main-outline-of-CSMs) 
+ [Publication and associated datasets](#Publication) 
+ [Citation](#Images) 
+ [Reporting Bugs](#🐛-Reporting-Bugs)

## Main features and requirements
* CSM was developed using R version > 4.3
* CSM depends on many packages to work. However not all suggested packages are required to perform specific tasks.
* CSM includes >100 functions organized in 6+1 modules.
* CSM provides several datasets to test function capabilities.
* You can install CSM in R using `devtools::install_github(Alvaro-LJ/CSM)`.

## Main outline of CSM
<img src="man/figures/CSM_OUTLINE.png" width="750" style = "align-items: center">

- MODULE 0: Data generation and formatting
  - Divide large images into tiles. Useful to work with whole slide tissue samples.
  - Perform pixel thresholding and quantification. Used to work with extra-cellular biomolecules.
  - Extract color channels from RGB images. May be required to work with HE, IHC or other histochemical staining techniques.
  - Perform cell segmentation and feature extraction from images.
  - Arrange cell feature matrices into an adequate format.
  - Arrange image metadata into an adequate format.
  - Perform general quality checks and check computational resources available.
  
- MODULE 1: Data normalization
  - Normalize cell feature expression data according to image and slide belonging.
  
- MODULE 2: Data thresholding
  - Define positive thresholds for features of interest. Thresholded data can be used to assign cell labels.
  
- MODULE 3: Cell phenotype labeling
  - Define cell phenotype labels based on thresholded features (see MODULE 2)
  - Define cell phenotype labels based on unsupervised clustering.
  - Define cell phenotype labels based on semi-supervised clustering.
  - Define cell phenotype labels based on user trained algorithms.
  - Check concordance between different cell phenotyping methods.
  - Assign cell phenotype labels by merging various methods.
  - Quantify cell phenotypes and find associations with image metadata.

- MODULE 4: Heterogeneity assessment
  - Calculate global cell composition heterogeneity
  - Calculate spatial heterogeneity using tiling approaches as well as texture features
  
- MODULE 5: Cell-Cell spatial associations
  - Calculate spatial associations between pairs of cell types and thriads.
  - Calculate spatial association occurring by chance.
  
- MODULE 6: Neighborhood analysis and tissue structures
  - Calculate cellular neighborhoods using various algorithms
  - Divide tissue into compartments according to a single cell type (for example Tumor/stromal compartments).

## Publication
To see examples of use of CSM you can have a look at our associated publication [*publication*](publicationURL).<br />
Also this publication has an associated GitHub [*repository*](CSM repositoryURL) where user can find test datasets and examples of use of CSM.

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

## 🐛 Reporting Bugs

If you encounter any bugs or unexpected behavior, please help us improve by reporting them! <br />

You can report issues directly through the [GitHub Issues page](https://github.com/Alvaro-LJ). When submitting a bug report, please include:

- A clear description of the problem
- Steps to reproduce the issue
- Screenshots or logs if applicable
- Your environment (OS, Rversion)

We appreciate your feedback and contributions!


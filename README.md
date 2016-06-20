# RAVEL
[![Build Status](https://travis-ci.org/Jfortin1/RAVEL.svg?branch=master)](https://travis-ci.org/Jfortin1/RAVEL)

### Imaging suite for the statistical analysis of structural MRIs. 

RAVEL stands for Removal of Artificial Voxel Effect by Linear regression, the main preprocessing function of the package. 

#### Installation

```{r}
library(devtools)
install_github("jfortin1/RAVEL")
```


#### preprocessRAVEL

The function `preprocessRAVEL` applies the RAVEL correction described in XX


#### preprocessWS

The function `preprocessWS` applies the White Stripe intensity normalization described in XX

#### Tissue Segmentation

There are different tissue segmentation algorithms available in R. My favorite is the FSL FAST segmentation via the [`fslr`](https://cran.r-project.org/web/packages/fslr/index.html) package. Note that the `fslr` package requires FSL to be installed on your machine; see the [FSL website](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/).

```{r}
library(fslr)
library(RAVELData)
have.fsl() # Should be TRUE if fsl is correctly installed
img.file <- 
img <- readNIfTI()
```

- fast
- athropos

#### RAVEL for longitudinal data

#### RAVEL for multimodal data





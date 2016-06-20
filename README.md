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

There are different tissue segmentation algorithms available in R. My favorite is the FSL FAST segmentation via the [`fslr`](https://cran.r-project.org/web/packages/fslr/index.html) package. Note that the `fslr` package requires FSL to be installed on your machine; see the [FSL website](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/). Here is an example how to perform segmentation on the JHU-MNI-ss template, after skull removal, that is included in the `RAVELData` package. First, let's make sure we have `fslr` correctly installed:

```{r}
library(fslr)
library(RAVELData)
have.fsl() # Should be TRUE if fsl is correctly installed
```

We first load the skull-stripped image to be segmented:

```{r}
img_path <- system.file(package="RAVELData", "data/JHU_MNI_SS_T1_Brain.nii.gz")
img <- readNIfTI(img_path)
ortho2(img, crosshairs=FALSE, mfrow=c(1,3), add.orient=FALSE, ylim=c(0,400))
```
The last line of code produces via the `ortho2` function from the `fslr` package the following visualization of the template:


<p align="center">
<img src="https://github.com/Jfortin1/RAVEL/blob/master/images/template.png" width="750"/>
</p>

We perform a 3-class tissue segmentation on the T1-w image with the FAST segmentation algorithm:

```{r}
seg <- fast(img, verbose=FALSE, opts="-t 1 -n 3") 
ortho2(seg, crosshairs=FALSE, mfrow=c(1,3), add.orient=FALSE)
```
<p align="center">
<img src="https://github.com/Jfortin1/RAVEL/blob/master/images/seg.png" width="750"/>
</p>


The object `seg` is an image that contains the segmentation labels `0,1,2` and `3` referring to Background, CSF, GM and WM voxels respectively. 


#### Intensity Inhomogeneity correction

To install `ANTsR`, please visit the [package GitHub page](https://github.com/stnava/ANTsR).

```{r}
library(ANTsR)
img_path <- system.file(package="RAVELData", "data/JHU_MNI_SS_T1.nii.gz")
img <- antsImageRead(img_path, 3)
img_n4 <- n4BiasFieldCorrection(img)
```

#### Registration to template

Tp perform a non-linear registration to the JHU-MNI-ss template, one can use the diffeomorphism algorithm via the `ANTsR` package. Note that we perform the registration with the skulls on. Here is an example where we register the scan1 from the `RAVELData` package to the JHU-MNI-ss template:

```{r}
template_path <- system.file(package="RAVELData", "data/JHU_MNI_SS_T1.nii.gz")
template    <- antsImageRead(template_path, 3)
brain_path <- system.file(package="RAVELData", "data/scan1.nii.gz")
brain <- antsImageRead(brain_path,3)
outprefix <- gsub(".nii.gz","",brain_path) # Prefix for the output files
output <- antsRegistration(fixed = template, moving = brain, typeofTransform = "SyN",  outprefix = outprefix)
brain_reg   <- antsImageClone(output$warpedmovout) # Registered brain
```

For visualization:

```{}
ortho2(ants2oro(brain), crosshairs=FALSE, mfrow=c(1,3), add.orient=FALSE)
ortho2(ants2oro(brain_reg), crosshairs=FALSE, mfrow=c(1,3), add.orient=FALSE)
```


  
#### Coregistration (for more than one modality)

#### RAVEL for longitudinal data

#### RAVEL for multimodal data





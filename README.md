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




### Preprocessing images

#### Registration to template

Tp perform a non-linear registration to the JHU-MNI-ss template, one can use the diffeomorphism algorithm via the `ANTsR` package. To install `ANTsR`, please visit the [package GitHub page](https://github.com/stnava/ANTsR). Note that we perform the registration with the skulls on. Here is an example where we register the scan1 from the `RAVELData` package to the JHU-MNI-ss template:

```{r}
library(ANTsR)
template_path <- system.file(package="RAVELData", "data/JHU_MNI_SS_T1.nii.gz")
template    <- antsImageRead(template_path, 3)
scan_path <- system.file(package="RAVELData", "data/scan1.nii.gz")
scan <- antsImageRead(scan_path,3)
outprefix <- gsub(".nii.gz","",scan_path) # Prefix for the output files
output <- antsRegistration(fixed = template, moving = scan, typeofTransform = "SyN",  outprefix = outprefix)
scan_reg   <- antsImageClone(output$warpedmovout) # Registered brain
```
The object `scan_reg` contains the scan registed to the template. Note that the object is in the `ANTsR` format. Since I prefer to work with the `oro.nifti` package, which is compatible with `flsr`, I convert the object to a `nifti` object using the function `ants2oro` as follows:

```{r}
scan_reg <- ants2oro(scan_reg)
```
I can save the registered brain in the NIfTi format using the `writeNIfTI` command:

```{r}
writeNIfTI(scan_reg, "scan_reg")
```
Since `scan_reg` is converted to a `nifti` object, we can use the function `ortho2` from the `fslr` package to visualize the scan: 

```{r}
ortho2(scan_reg, crosshairs=FALSE, mfrow=c(1,3), add.orient=FALSE)
```

#### Intensity inhomogeneity correction

We perform intensity inhomogeneity correction on the registered scan using the N4 Correction from the `ANTsR` package:

```{r}
scan_reg <- oro2ants(scan_reg) # Convert to ANTsR object
scan_reg_n4 <- n4BiasFieldCorrection(scan_reg)
scan_reg_n4 <- ants2oro(scan_reg_n4) # Conversion to nifti object for further processing
```

#### Skull stripping

```{r}
template_brain_mask_path <- system.file(package="RAVELData", "data/JHU_MNI_SS_T1_Brain_Mask.nii.gz")
template_brain_mask <- readNIfTI(template_brain_mask_path, reorient=FALSE)
scan_reg_n4_brain <- niftiarr(scan_reg_n4, scan_reg_n4*template_brain_mask)
```

Visualization:

```{}
ortho2(scan_reg_n4_brain, crosshairs=FALSE, mfrow=c(1,3), add.orient=FALSE)
```
 
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

  
  

  
  
#### Coregistration (for more than one modality)

#### RAVEL for longitudinal data

#### RAVEL for multimodal data





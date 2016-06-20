# RAVEL
### An imaging suite for the statistical analysis of structural MRIs. 
--------

**Creator**: Jean-Philippe Fortin, jeanphi@mail.med.upenn.edu

**Authors**: Jean-Philippe Fortin, John Muschelli, Russell T. Shinohara

##### Software status

| Resource:      | Travis CI     |
| -------------  |  ------------- |
| _Platform:_    | _Linux_       |
| R CMD check    | <a href="https://travis-ci.org/Jfortin1/RAVEL"><img src="https://travis-ci.org/Jfortin1/RAVEL.svg?branch=master" alt="Build status"></a> |



### 0. Introduction

RAVEL stands for Removal of Artificial Voxel Effect by Linear regression, the main preprocessing function of the package. In Section 1, we explain how to use the RAVEL algorithm as well as other intensity normalization techniques. In Section 2, we present a pre-normalization preprocessing pipeline from raw images to processed images ready for intensity normalization. In Section XX, we present different tools for the post-normalizations statistical analysis. In Section XX, we present additional functions that help the visualization of images and statistical results. 

##### Installation

```{r}
library(devtools)
install_github("jfortin1/RAVEL")
```

### 1. Intensity normalization with WhiteStripe and RAVEL



The function `preprocessRAVEL` applies the RAVEL correction described in XX. The function `preprocessWS` applies the White Stripe intensity normalization described in XX

##### Creation of a control region for the RAVEL algorithm




### 2. Preprocessing images

We present a pre-normalization preprocessing pipeline implemented in the R software, from raw images to images ready for intensity normalization and statistical analysis. Once the images are preprocessed, users can apply their favorite intensity normalization and the scan-effect correction tool RAVEL as presented in Section 1 above. We present a preprocessing pipeline that uses the R packages `ANTsR` and `fslr`. While we have chosen to use a specific template space (JHU-MNI-ss), a specific registration (non-linear diffeomorphic registration) and a specific tissue segmentation (FSL FAST), users can choose other algorithms prior to intensity normalization and in order for RAVEL to work. The only requirement is that the images are registered to the same template space. 

##### Prelude

To preprocess the images, we use the packages `fslr` and `ANTsR`. The package `fslr` is available on CRAN, and requires FSL to be installed on your machine; see the [FSL website](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/). To install `ANTsR`, please visit the [package GitHub page](https://github.com/stnava/ANTsR). The JHU-MNI-ss atlas is included in the package `RAVELData`, together with 4 MPRage T1-w scans. The package is available on GitHub at [https://github.com/Jfortin1/RAVELData](https://github.com/Jfortin1/RAVELData). Once the packages are properly installed, we need to load them into R:

```{r}
library(fslr)
library(ANTsR)
library(RAVELData)
have.fsl() # Should be TRUE if fsl is correctly installed
```

and let's specify the path for the different files that we will need:
```{r}
# JHU-MNI-ss template:
template_path <- system.file(package="RAVELData", "data/JHU_MNI_SS_T1.nii.gz") 
# JHU-MNI-ss template brain mask:
template_brain_mask_path <- system.file(package="RAVELData", "data/JHU_MNI_SS_T1_Brain_Mask.nii.gz") 
# Example of T1-w MPRAGE image
scan_path <- system.file(package="RAVELData", "data/scan1.nii.gz")
```

##### Registration to template

Tp perform a non-linear registration to the JHU-MNI-ss template, one can use the diffeomorphism algorithm via the `ANTsR` package.  Note that we perform the registration with the skulls on. Here is an example where we register the scan1 from the `RAVELData` package to the JHU-MNI-ss template:

```{r}
template    <- antsImageRead(template_path, 3)
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

##### Intensity inhomogeneity correction

We perform intensity inhomogeneity correction on the registered scan using the N4 Correction from the `ANTsR` package:

```{r}
scan_reg <- oro2ants(scan_reg) # Convert to ANTsR object
scan_reg_n4 <- n4BiasFieldCorrection(scan_reg)
scan_reg_n4 <- ants2oro(scan_reg_n4) # Conversion to nifti object for further processing
```

##### Skull stripping

```{r}
template_brain_mask <- readNIfTI(template_brain_mask_path, reorient=FALSE)
scan_reg_n4_brain <- niftiarr(scan_reg_n4, scan_reg_n4*template_brain_mask)
```

Visualization:

```{}
ortho2(scan_reg_n4_brain, crosshairs=FALSE, mfrow=c(1,3), add.orient=FALSE)
```
 
##### Tissue Segmentation

There are different tissue segmentation algorithms available in R. My favorite is the FSL FAST segmentation via the [`fslr`](https://cran.r-project.org/web/packages/fslr/index.html) package. 

Let's produce the tissue segmentation for the `scan_reg_n4_brain` scan above:

```{r}
ortho2(scan_reg_n4_brain, crosshairs=FALSE, mfrow=c(1,3), add.orient=FALSE, ylim=c(0,400))
```
The last line of code produces via the `ortho2` function from the `fslr` package the following visualization of the template:


<p align="center">
<img src="https://github.com/Jfortin1/RAVEL/blob/master/images/template.png" width="600"/>
</p>

We perform a 3-class tissue segmentation on the T1-w image with the FAST segmentation algorithm:

```{r}
scan_reg_n4_brain_seg <- fast(scan_reg_n4_brain, verbose=FALSE, opts="-t 1 -n 3") 
ortho2(scan_reg_n4_brain_seg, crosshairs=FALSE, mfrow=c(1,3), add.orient=FALSE)
```
<p align="center">
<img src="https://github.com/Jfortin1/RAVEL/blob/master/images/seg.png" width="600"/>
</p>


The object `scan_reg_n4_brain_seg` is an image that contains the segmentation labels `0,1,2` and `3` referring to Background, CSF, GM and WM voxels respectively. 

  
##### Creation of a tissue mask

Suppose we want to create a mask for CSF.

```{r}
scan_reg_n4_brain_csf_mask <- scan_reg_n4_brain_seg
scan_reg_n4_brain_csf_mask[scan_reg_n4_brain_csf_mask!=1] <- 0
ortho2(scan_reg_n4_brain_csf_mask, crosshairs=FALSE, mfrow=c(1,3), add.orient=FALSE)
```
We use the fact that the file `scan_reg_n4_brain_seg` is equal to 1 for CSF, 2 for GM and 3 for WM. FOr instance, a WM mask could be created as follows:

```{r}
scan_reg_n4_brain_wm_mask <- scan_reg_n4_brain_seg
scan_reg_n4_brain_wm_mask[scan_reg_n4_brain_wm_mask!=3] <- 0
ortho2(scan_reg_n4_brain_wm_mask, crosshairs=FALSE, mfrow=c(1,3), add.orient=FALSE)
```



### 3. Advance stuff and extensions
  
##### Coregistration (for more than one modality)

##### RAVEL for longitudinal data

##### RAVEL for multimodal data





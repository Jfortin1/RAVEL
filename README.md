

RAVEL <img src="sticker.png" width = "150" align="right" />
===========================================================

### Intensity normalizations for structural MRIs

**Creator**: Jean-Philippe Fortin, <fortin946@gmail.com>

**Authors**: Jean-Philippe Fortin, John Muschelli

**License**: GPL-2

##### Software status

| Resource:   | Travis CI                                                                                                                                |
|-------------|------------------------------------------------------------------------------------------------------------------------------------------|
| Platform:   | OSX                                                                                                                                      |
| R CMD check | <a href="https://travis-ci.org/Jfortin1/RAVEL"><img src="https://travis-ci.org/Jfortin1/RAVEL.svg?branch=master" alt="Build status"></a> |

##### References

<table style="width:58%;">
<colgroup>
<col width="19%" />
<col width="19%" />
<col width="19%" />
</colgroup>
<thead>
<tr class="header">
<th>Method</th>
<th>Citation</th>
<th>Paper Link</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>RAVEL</td>
<td>Jean-Philippe Fortin, Elizabeth M Sweeney, John Muschelli, Ciprian M Crainiceanu, Russell T Shinohara, Alzheimer's Disease Neuroimaging Initiative, et al. <strong>Removing inter-subject technical variability in magnetic resonance imaging studies</strong>. NeuroImage, 132:198–212, 2016.</td>
<td><a href="http://www.sciencedirect.com/science/article/pii/S1053811916001452">Link</a></td>
</tr>
<tr class="even">
<td>WhiteStripe</td>
<td>Russell T Shinohara, Elizabeth M Sweeney, Jeff Goldsmith, Navid Shiee, Farrah J Mateen, Peter A Calabresi, Samson Jarso, Dzung L Pham, Daniel S Reich, Ciprian M Crainiceanu, Australian Imaging Biomarkers Lifestyle Flagship Study of Ageing, and Alzheimer’s Disease Neuroimaging Initiative. <strong>Statistical normalization techniques for magnetic resonance imaging</strong>. Neuroimage Clin, 6:9–19, 2014.</td>
<td><a href="http://www.sciencedirect.com/science/article/pii/S221315821400117X">Link</a></td>
</tr>
</tbody>
</table>

Table of content
----------------

-   [1. Introduction](#id-section1)
-   [2. Image preprocessing](#id-section2)
-   [3. Intensity normalization and RAVEL correction](#id-section3)

1. Introduction
---------------

RAVEL is an R package that combines the preprocessing and statistical analysis of magnetic resonance imaging (MRI) datasets within one framework. Users can start with raw images in the NIfTI format, and end up with a variety of statistical results associated with voxels and regions of interest (ROI) in the brain. RAVEL stands for *Removal of Artificial Voxel Effect by Linear regression*, the main preprocessing function of the package that allows an effective removal of between-scan unwanted variation. We have shown in [a recent paper](http://www.sciencedirect.com/science/article/pii/S1053811916001452) that RAVEL improves significantly population-wide statistical inference. RAVEL is now part of the [Neuroconductor project](https://neuroconductor.org/).

##### Installation

You can install RAVEL from github with:

``` r
# install.packages("devtools")
devtools::install_github("jfortin1/RAVEL")
```

2. Image preprocessing
----------------------

We present a pre-normalization preprocessing pipeline implemented in the R software, from raw images to images ready for intensity normalization and statistical analysis. Once the images are preprocessed, users can apply their favorite intensity normalization and the scan-effect correction tool RAVEL as presented in Section 1 above. We present a preprocessing pipeline that uses the R packages `ANTsR` and `fslr`. While we have chosen to use a specific template space (JHU-MNI-ss), a specific registration (non-linear diffeomorphic registration) and a specific tissue segmentation (FSL FAST), users can choose other algorithms prior to intensity normalization and in order for RAVEL to work. The only requirement is that the images are registered to the same template space.

### 2.1. Prelude

To preprocess the images, we use the packages `fslr` and `ANTsR`. The package `fslr` is available on CRAN, and requires FSL to be installed on your machine; see the [FSL website](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) for installation. For `ANTsR`, we recommend to install the latest stable version available at the ANTsR [GitHub page](https://github.com/stnava/ANTsR/releases/). The version used for this vignette was `ANTsR_0.3.2.tgz`. For the template space, we use the JHU-MNI-ss atlas (see Section 1.2) included in the `EveTemplate` package, available on GitHub at <https://github.com/Jfortin1/EveTemplate>. For data examples, we use 4 T1-w scans from the package `RAVELData` available on GitHub at <https://github.com/Jfortin1/RAVELData>. Once the packages are properly installed, we are ready to start our preprocessing of T1-w images. We first load the packages into R:

``` r
library(fslr)
library(ANTsR)
library(RAVELData)
library(EveTemplate)
have.fsl() # Should be TRUE if fsl is correctly installed
```

and let's specify the path for the different files that we will need:

``` r
# JHU-MNI-ss template:
library(EveTemplate)
template_path <- getEvePath("T1")
# JHU-MNI-ss template brain mask:
template_brain_mask_path <- getEvePath("Brain_Mask")
# Example of T1-w MPRAGE image
scan_path <- system.file(package="RAVELData", "data/scan1.nii.gz")
```

### 2.2. JHU-MNI-ss template (*EVE* atlas)

### 2.3. Registration to template

Tp perform a non-linear registration to the JHU-MNI-ss template, one can use the diffeomorphism algorithm via the `ANTsR` package. Note that we perform the registration with the skulls on. Here is an example where we register the scan1 from the `RAVELData` package to the JHU-MNI-ss template:

``` r
library(ANTsRCore)
library(ANTsR)
template    <- antsImageRead(template_path, 3)
scan <- antsImageRead(scan_path,3)
outprefix <- gsub(".nii.gz","",scan_path) # Prefix for the output files
output <- antsRegistration(fixed = template, moving = scan, typeofTransform = "SyN",  outprefix = outprefix)
scan_reg   <- antsImageClone(output$warpedmovout) # Registered brain
```

The object `scan_reg` contains the scan registed to the template. Note that the object is in the `ANTsR` format. Since I prefer to work with the `oro.nifti` package, which is compatible with `flsr`, I convert the object to a `nifti` object using the function `ants2oro` as follows:

``` r
# devtools::install_github("muschellij2/extrantsr")
# or
# source("https://neuroconductor.org/sites/default/files/neurocLite.R")
# neuro_install("extrantsr")
library(extrantsr)
scan_reg <- extrantsr::ants2oro(scan_reg)
```

I can save the registered brain in the NIfTi format using the `writeNIfTI` command:

``` r
writeNIfTI(scan_reg, "scan_reg")
```

Since `scan_reg` is converted to a `nifti` object, we can use the function `ortho2` from the `fslr` package to visualize the scan:

``` r
ortho2(scan_reg, crosshairs=FALSE, mfrow=c(1,3), add.orient=FALSE)
```

### 2.4. Intensity inhomogeneity correction

We perform intensity inhomogeneity correction on the registered scan using the N4 Correction from the `ANTsR` package:

``` r
scan_reg <- extrantsr::oro2ants(scan_reg) # Convert to ANTsR object
scan_reg_n4 <- n4BiasFieldCorrection(scan_reg)
scan_reg_n4 <- extrantsr::ants2oro(scan_reg_n4) # Conversion to nifti object for further processing
```

### 2.5. Skull stripping

``` r
template_brain_mask <- readNIfTI(template_brain_mask_path, reorient=FALSE)
scan_reg_n4_brain <- niftiarr(scan_reg_n4, scan_reg_n4*template_brain_mask)
ortho2(scan_reg_n4_brain, crosshairs=FALSE, mfrow=c(1,3), add.orient=FALSE)
```

### 2.6. Tissue Segmentation

There are different tissue segmentation algorithms available in R. My favorite is the FSL FAST segmentation via the [`fslr`](https://cran.r-project.org/web/packages/fslr/index.html) package. Let's produce the tissue segm
## ----style, eval=TRUE, echo=FALSE, results="asis"---------------------------------------
BiocStyle::latex()

## ---------------------------------------------------------------------------------------
library(RAVEL)
library(RAVELData)
datadir <- system.file("data",package="RAVELData")
files <- list.files(datadir, pattern="csf_mask.nii.gz", full.names=TRUE)
masks <- pblapply(files, readNIfTI)
control.mask<- mask_intersect(masks)

## ---------------------------------------------------------------------------------------
files <- list.files(datadir, pattern="processed.nii.gz", full.names=TRUE)
template    <- file.path(datadir, "JHU_MNI_SS_T1.nii.gz")
brain.mask <- file.path(datadir, "JHU_MNI_SS_T1_Brain_Mask.nii.gz")
RAVEL(input.files = files, brain.mask = brain.mask, control.mask = control.mask,WhiteStripe=FALSE)

## ----sessionInfo, results="asis", echo=FALSE, eval=TRUE---------------------------------
toLatex(sessionInfo())


library(RAVEL)
library(RAVELData)
library(EveTemplate)
dir <- file.path(find.package("RAVELData"), "extdata")
input.files  <- list.files(dir, full.names=TRUE, pattern="processed.nii.gz")
control.mask <- list.files(dir, full.names=TRUE, pattern="mask_intersection.nii.gz")
brain.mask   <- EveTemplate::getEvePath("Brain_Mask")


Y.raw <- normalizeRaw(input.files=input.files,
	brain.mask=brain.mask,
	returnMatrix=TRUE
)

Y.ravel <- normalizeRAVEL(input.files=input.files,
	control.mask=control.mask,
	brain.mask=brain.mask,
	k=1, 
	returnMatrix=TRUE,
	WhiteStripe=FALSE
)

# Example using biological covariates:
x   <- rnorm(4)
mod <- model.matrix(~x)
Y.ravel.mod <- normalizeRAVEL(input.files=input.files,
	control.mask=control.mask,
	brain.mask=brain.mask,
	k=1, 
	mod=mod,
	returnMatrix=TRUE,
	WhiteStripe=FALSE
)
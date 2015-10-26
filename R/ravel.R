
library(oro.nifti)
library(fslr)


# WhiteStripe option not yet implemented
# Assuming images are registered and normalized beforehand
ravel <- function(input.files, output.files=NULL, WhiteStripe=FALSE, brain.mask=NULL, control.mask, k=1, verbose=TRUE){
	if (!is.null(brain.mask)){
		brain.mask <- readNIfTI(brain.mask, reorient=FALSE)
		brain.indices <- brain.mask==1
	} 


	if (is.null(output.files)){
		output.files <- gsub(".nii.gz|.nii","_ravel.nii.gz", input.files)
	}


	if (verbose){
		cat("[ravel] Creating the voxel intensities matrix V \n")
	}
	# Matrix of voxel intensities:
	V <- do.call(cbind,lapply(input.files, function(x){
		brain <- readNIfTI(x, reorient=FALSE)
		if (!is.null(brain.mask)){
			brain <- as.vector(brain[brain.indices])
		}
		brain
	}))

	if (verbose){
		cat("[ravel] Creating the control voxel matrix Vc \n")
	}
	# Submatrix of control voxels:
	control.indices <- as.vector(readNIfTI(brain.mask, reorient=FALSE)==1)
	Vc <- V[control.indices,,drop=FALSE]
	if (verbose){
		cat("[ravel] Estimating the unwanted factors Z \n")
	}
	Z  <- svd(Vc)$v[, 1:k, drop=FALSE] # Unwanted factors



	if (verbose){
		cat("[ravel] Performing RAVEL correction \n")
	}
	V.norm <- .ravel.correction(V,Z)


	if (verbose){
		cat("[ravel] Writing on the disk the RAVEL-corrected images \n")
	}
	template <- readNIfTI(input.files[1], reorient=FALSE)
	if (!is.null(brain.mask)){
		template[!brain.indices] <- 0
		template[brain.indices]  <- 1

	for (i in 1:ncol(V.norm)){
		.write.brain(brain.norm = V.norm[,i], output.file = output.files[i], template=template)
		if (verbose){print(i)}
	}

}


# RAVEL correction procedure:
.ravel.correction <- function(V, Z){
	means <- rowMeans(Y)
	beta   <- solve(t(X) %*% X) %*% t(X) %*% t(Y)
	fitted <- t(X %*% beta)
	res   <- Y - fitted
	res   <- res + means
	res
}


.write.brain <- function(brain.norm, output.file, template){
	if (!is.null(brain.mask)){
		template[brain.indices] <- brain.norm
	} else {
		template <- .niftiarr(template, brain.norm)
	}
	template <- .cal_img(template)
	writeNIfTI(template, output.file)
}


# Function from fslr package (written by John Muschelli)
.niftiarr <-function (img, arr) {
    x = img
    if (!is(arr, "array")) {
        arr = array(arr, dim = dim(img))
    }
    arr = as(arr, "array")
    class(arr) = "numeric"
    stopifnot(all(dim(arr) == dim(img)))
    x@.Data = arr
    x = cal_img(x)
    x = zero_trans(x)
}


# Function from fslr package (written by John Muschelli)
.cal_img <- function (img) {
    cmax = max(img, na.rm = TRUE)
    cmax = ifelse(is.finite(cmax), cmax, 0)
    cmin = min(img, na.rm = TRUE)
    cmin = ifelse(is.finite(cmin), cmin, 0)
    cal.max(img) = cmax
    cal.min(img) = cmin
    img
}






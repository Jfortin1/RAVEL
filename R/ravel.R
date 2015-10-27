library(oro.nifti)
library(fslr)
library(WhiteStripe)


# WhiteStripe option not yet implemented
# Assuming images are registered and normalized beforehand
ravel <- function(input.files, output.files=NULL, brain.mask=NULL, control.mask, WhiteStripe=FALSE, WhiteStripe_Type="T1",  k=1, verbose=TRUE){
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

	if (WhiteStripe & WhiteStripe_Type!="T1"){
		cat("WhiteStripe is not performed -- image modality not supported \n")
	} else if (WhiteStripe & WhiteStripe_Type=="T1"){
		if (verbose){
			cat("[ravel] Performing White Stripe intensity normalization \n")
		}

		# Performing White Stripe normalization: 
		V <- do.call(cbind,lapply(input.files, function(x){
			brain   <- readNIfTI(x, reorient=FALSE)
			indices <- whitestripe(brain, type=WhiteStripe_Type)
			brain    <- whitestripe_norm(brain, indices$whitestripe.ind)
			if (!is.null(brain.mask)){
				brain <- as.vector(brain[brain.indices])
			}
			brain
		}))
	}
	type="T1"
  	
  	
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
		template <- fslr::niftiarr(template, brain.norm)
	}
	template <- oro.nifti::cal_img(template)
	writeNIfTI(template, output.file)
}











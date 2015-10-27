
# WhiteStripe option not yet implemented
# Assuming images are registered and normalized beforehand
RAVEL <- function(input.files, output.files=NULL, brain.mask=NULL, control.mask=NULL, WhiteStripe=FALSE, WhiteStripe_Type="T1",  k=1, verbose=TRUE){
	# RAVEL correction procedure:
	.ravel.correction <- function(V, Z){
		means <- rowMeans(V)
		beta   <- solve(t(Z) %*% Z) %*% t(Z) %*% t(V)
		fitted <- t(Z %*% beta)
		res   <- V - fitted
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
		output.file <- gsub(".nii.gz|.nii", "", output.file)
		writeNIfTI(template, output.file)
	}



	if (!is.null(brain.mask)){
		brain.mask <- readNIfTI(brain.mask, reorient=FALSE)
		brain.indices <- brain.mask==1
	} 

	if (is.null(output.files)){
		output.files <- gsub(".nii.gz|.nii","_RAVEL.nii.gz", input.files)
	}


	
	# Matrix of voxel intensities:
	if (!WhiteStripe){
		if (verbose){
			cat("[RAVEL] Creating the voxel intensities matrix V \n")
		}
		V <- do.call(cbind,lapply(input.files, function(x){
			brain <- readNIfTI(x, reorient=FALSE)
			if (!is.null(brain.mask)){
				brain <- as.vector(brain[brain.indices])
			}
			brain
		}))
	}

	

	if (WhiteStripe & WhiteStripe_Type!="T1"){
		cat("WhiteStripe is not performed -- image modality not supported \n")
	} else if (WhiteStripe & WhiteStripe_Type=="T1"){
		if (verbose){
			cat("[RAVEL] Performing White Stripe intensity normalization \n")
		}
		if (verbose){
			cat("[RAVEL] Creating the voxel intensities matrix V \n")
		}
		# Performing White Stripe normalization: 
		V <- do.call(cbind,lapply(input.files, function(x){
			brain   <- readNIfTI(x, reorient=FALSE)
			indices <- whitestripe(brain, type=WhiteStripe_Type, verbose=FALSE)
			brain    <- whitestripe_norm(brain, indices$whitestripe.ind)
			if (!is.null(brain.mask)){
				brain <- as.vector(brain[brain.indices])
			}
			brain
		}))
	}
	
  	
  	
	if (verbose){
		cat("[RAVEL] Creating the control voxel matrix Vc \n")
	}
	# Submatrix of control voxels:
	control.indices <- readNIfTI(control.mask, reorient=FALSE)==1
	if (!is.null(brain.mask)){
		control.indices <- control.indices[brain.mask==1]
	}


	Vc <- V[control.indices,,drop=FALSE]
	if (verbose){
		cat("[RAVEL] Estimating the unwanted factors Z \n")
	}
	Z  <- svd(Vc)$v[, 1:k, drop=FALSE] # Unwanted factors



	if (verbose){
		cat("[RAVEL] Performing RAVEL correction \n")
	}
	V.norm <- .ravel.correction(V,Z)


	if (verbose){
		cat("[RAVEL] Writing out the corrected images \n")
	}
	template <- readNIfTI(input.files[1], reorient=FALSE)
	if (!is.null(brain.mask)){
		template[!brain.indices] <- 0
		template[brain.indices]  <- 1
	}

	for (i in 1:ncol(V.norm)){
		.write.brain(brain.norm = V.norm[,i], output.file = output.files[i], template=template)
		#if (verbose){print(i)}
	}

}














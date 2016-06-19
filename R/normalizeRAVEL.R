# Assuming images are registered and normalized beforehand
normalizeRAVEL <- function(input.files, output.files=NULL, brain.mask=NULL, control.mask=NULL, WhiteStripe=FALSE, WhiteStripe_Type="T1",  k=1, verbose=TRUE, writeToDisk=FALSE, returnObject=FALSE){
	
	# RAVEL correction procedure:
	if (WhiteStripe & WhiteStripe_Type!="T1") stop("Only image modality T1 is supported at the moment for WhiteStripe")
	if (!verbose) pboptions(type="none") 

	if (!is.null(brain.mask)){
		brain.mask <- readNIfTI(brain.mask, reorient=FALSE)
		brain.indices <- brain.mask==1
	} else {
		stop("brain.mask must be provided.")
	}

	if (is.null(output.files)){
		output.files <- gsub(".nii.gz|.nii","_RAVEL.nii.gz", input.files)
	}

	.ravel_correction <- function(V, Z){
		means <- rowMeans(V)
		beta   <- solve(t(Z) %*% Z) %*% t(Z) %*% t(V)
		fitted <- t(Z %*% beta)
		res   <- V - fitted
		res   <- res + means
		res
	}

	.write_brain <- function(brain.norm, output.file, template){
		if (!is.null(brain.mask)){
			template[brain.indices] <- brain.norm
		} else {
			template <- fslr::niftiarr(template, brain.norm)
		}
		template <- oro.nifti::cal_img(template)
		output.file <- gsub(".nii.gz|.nii", "", output.file)
		writeNIfTI(template, output.file)
	}

	cat("[RAVEL] Creating the voxel intensities matrix V. \n")
	if (WhiteSripe) cat("[RAVEL] WhiteStripe intensity normalization is applied to each scan. \n")
	if (!WhiteSripe) cat("[RAVEL] WhiteStripe intensity normalization not applied.  \n")


	# Matrix of voxel intensities:
	V <- pblapply(input.files, function(x){
		brain <- readNIfTI(x, reorient=FALSE)
		if (WhiteStripe){
			indices <- whitestripe(brain, type=WhiteStripe_Type, verbose=FALSE)
			brain   <- whitestripe_norm(brain, indices$whitestripe.ind)
		}
		if (!is.null(brain.mask)){
			brain <- as.vector(brain[brain.indices])
		}
		brain
	})
	V <- do.call(cbind, V)


  	# Submatrix of control voxels:
	if (verbose) cat("[RAVEL] Creating the control voxel matrix Vc. \n")
	if (class(control.mask)=="nifti"){
		control.indices <- control.mask==1
	} else {
		control.indices <- readNIfTI(control.mask, reorient=FALSE)==1
	}
	control.indices <- control.indices[brain.mask==1]
	Vc <- V[control.indices,,drop=FALSE]




	if (verbose) cat("[RAVEL] Estimating the unwanted factors Z. \n")
	Z  <- svd(Vc)$v[, 1:k, drop=FALSE] # Unwanted factors


	if (verbose) cat("[RAVEL] Performing RAVEL correction \n")
	V.norm <- .ravel_correction(V,Z)


	if (verbose & writeToDisk){
		cat("[RAVEL] Writing out the corrected images \n")
		template <- readNIfTI(input.files[1], reorient=FALSE)
		template[!brain.indices] <- 0
		template[brain.indices]  <- 1
		pblapply(1:ncol(V.norm), function(i){
			.write_brain(brain.norm = V.norm[,i], output.file = output.files[i], template=template)
		})
	} 
	if (returnObject){
		return(V)
	}
}

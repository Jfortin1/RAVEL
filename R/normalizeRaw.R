# Assuming images are registered and normalized beforehand
normalizeRAW <- function(input.files, output.files=NULL, brain.mask=NULL,  verbose=TRUE, writeToDisk=FALSE, returnMatrix=TRUE){
	
	# RAVEL correction procedure:
	if (!verbose) pboptions(type="none") 

	if (!is.null(brain.mask)){
		brain.mask <- readNIfTI(brain.mask, reorient=FALSE)
		brain.indices <- brain.mask==1
	} else {
		stop("brain.mask must be provided.")
	}

	if (is.null(output.files)){
		output.files <- gsub(".nii.gz|.nii","_RAW.nii.gz", input.files)
	}

	

	cat("[preprocessRAVEL] Creating the voxel intensities matrix V. \n")


	# Matrix of voxel intensities:
	V <- pblapply(input.files, function(x){
		brain <- readNIfTI(x, reorient=FALSE)
		if (!is.null(brain.mask)){
			brain <- as.vector(brain[brain.indices])
		}
		brain
	})
	V <- do.call(cbind, V)



	if (verbose & writeToDisk){
		cat("[preprocessRAVEL] Writing out the corrected images \n")
		template <- readNIfTI(input.files[1], reorient=FALSE)
		template[!brain.indices] <- 0
		template[brain.indices]  <- 1
		pblapply(1:ncol(V.norm), function(i){
			.write_brain(brain.norm = V.norm[,i], output.file = output.files[i], template=template)
		})
	} 
	if (returnMatrix){
		return(V)
	}
}

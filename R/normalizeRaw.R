# Assuming images are registered and normalized beforehand
normalizeRaw <- function(input.files, output.files=NULL, brain.mask=NULL, writeToDisk=FALSE, returnMatrix=TRUE, verbose=TRUE){
	
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

	cat("[normalizeRaw] Creating the voxel intensities matrix V. \n")


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
		cat("[normalizeRaw] Writing out the corrected images \n")
		pblapply(1:ncol(V), function(i){
			.write_brain(brain.norm = V[,i], output.file = output.files[i], brain.mask=brain.mask)
		})
	} 
	if (returnMatrix){
		return(V)
	}
}

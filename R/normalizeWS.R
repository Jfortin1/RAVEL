# Assuming images are registered and normalized beforehand
normalizeWS <- function(input.files, output.files=NULL, brain.mask=NULL, 
	WhiteStripe_Type=c("T1", "T2", "FLAIR"), writeToDisk=FALSE, returnMatrix=FALSE, verbose=TRUE){
	
	WhiteStripe_Type <- match.arg(WhiteStripe_Type)
	if (WhiteStripe_Type=="FLAIR") WhiteStripe_Type <- "T2"
	# RAVEL correction procedure:
	if (!verbose) pboptions(type="none") 

	if (!is.null(brain.mask)){
		brain.mask <- readNIfTI(brain.mask, reorient=FALSE)
		brain.indices <- brain.mask==1
	} else {
		stop("brain.mask must be provided.")
	}

	if (is.null(output.files)){
		output.files <- gsub(".nii.gz|.nii","_WS.nii.gz", input.files)
	}
	
    cat("[normalizeWS] WhiteStripe intensity normalization is applied to each scan. \n")
	# Matrix of voxel intensities:
	V <- pblapply(input.files, function(x){
		brain <- readNIfTI(x, reorient=FALSE)
        indices <- whitestripe(brain, type=WhiteStripe_Type, verbose=FALSE)
        brain   <- whitestripe_norm(brain, indices$whitestripe.ind)
        brain <- as.vector(brain[brain.indices])
		brain
	})
	V <- do.call(cbind, V)

	if (verbose & writeToDisk){
		cat("[normalizeWS] Writing out the corrected images \n")
		pblapply(1:ncol(V), function(i){
			.write_brain(brain.norm = V[,i], output.file = output.files[i], brain.mask=brain.mask)
		})
	} 
	if (returnMatrix){
		return(V)
	}	
}
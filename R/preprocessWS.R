# Assuming images are registered and normalized beforehand
normalizeWS <- function(input.files, output.files=NULL, brain.mask=NULL, WhiteStripe=FALSE, WhiteStripe_Type="T1", verbose=TRUE, writeToDisk=FALSE, returnObject=FALSE){
	
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
		output.files <- gsub(".nii.gz|.nii","_WS.nii.gz", input.files)
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

	
	if (WhiteSripe) cat("[RAVEL] WhiteStripe intensity normalization is applied to each scan. \n")
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

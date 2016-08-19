.getOrdinaryT <- function(fit){
		fit$coef/fit$stdev.unscaled/fit$sigma
}

getTMap <- function(V, design, columnIndex = NULL, weights=NULL, brain.mask, eBayes = TRUE, output.file="Tmap.nii.gz",  
	returnObject = TRUE, writeToDisk = FALSE, verbose = TRUE){
	if (is.null(columnIndex)){
		stop("columnIndex must be provided.")
	} else {
		if (columnIndex > ncol(design)){
			stop("columnIndex must be smaller than the number of columns of the design matrix.")
		}
	}
	if (!is.null(brain.mask)){
		brain.mask <- readNIfTI(brain.mask, reorient=FALSE)
		brain.indices <- brain.mask==1
	} else {
		stop("brain.mask must be provided.")
	}
	if (nrow(V) != sum(brain.indices)){
		stop("Number of rows of V does not agree with the brain mask.")
	}

	fit <- lmFit(V, design = design, weights = weights)
	Tmat <- .getOrdinaryT(fit)
	# Need to add some p-values maybe ?
	if (eBayes){
		fit <- eBayes(fit)
	}

	if (writeToDisk){
		if (verbose) cat("[getTMap] Writing out the T map \n")
		.write_brain(brain.norm = Tmat[, columnIndex], output.file = output.file, brain.mask=brain.mask)
	} 

	if (returnObject & !eBayes ){
		return(list(fit=fit, Tmat=Tmat))
	} else {
		return(list(fit=fit))
	}


}
.getOrdinaryT <- function(fit){
		fit$coef/fit$stdev.unscaled/fit$sigma
}



#' Compute t-statistic map with or without Empirical Bayes.
#' 
#' Compute t-statistic map with or without Empirical Bayes.
#' 
#' 
#' @param V Matrix of voxels intensities (rows are voxels, columns are
#' subjects).
#' @param design Design matrix for the covariates of interest, with rows
#' corresponding to subjects and columns to coefficients to be estimated.
#' @param columnIndex Numeric indicating which covariate of the design matrix
#' should be considered.
#' @param weights Non-negative observation weights. Can be a numeric matrix of
#' individual weivhts, of same size as the voxel matrix, or a numeric vector of
#' subject weights with length equal to 'ncol' of the voxel matrix, or a nmeric
#' vector of voxel weights with length equal to 'nrow' of the voxel matrix.
#' @param brain.mask Path of a NIfTI file mask to specify subset of voxels for
#' which the t-statistics should be computed.
#' @param eBayes Should the limma Empirical Bayes method be used to compute
#' moderated t-statistics?
#' @param output.file Filename of the NIfTI file to be saved to disk,
#' containing the computed t-statistics.
#' @param returnObject Should the t-statistics be returned as a matrix?
#' @param writeToDisk Should the t-statistics image be saved to the disk?
#' @param verbose Should messages be printed?
#' @return If \code{returnObject} is \code{TRUE}, a matrix of the t-statitics
#' is returned, and if \code{writeToDisk} is \code{TRUE}, the t-statistics
#' image is saved to disk as a NIfTI file.
#' @author Jean-Philippe Fortin
#' @importFrom limma lmFit eBayes
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

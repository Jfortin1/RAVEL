# Assuming images are registered and normalized beforehand


#' RAVEL intensity normalization
#' 
#' RAVEL intensity normalization.
#' 
#' 
#' @param input.files Vector of filenames for the input images. Must be NIfTI
#' files.
#' @param output.files Optional vector of filenames for the output images.  By
#' default, will be the \code{input.files} with suffix "RAVEL".
#' @param brain.mask Filename for the brain binary mask specifying the template
#' space brain.  Must be a NIfTI file.
#' @param control.mask Filename for the control region binary mask to be used
#' in RAVEL.  Must be a NIfTI file.
#' @param WhiteStripe Should White Stripe intensity normalization be performed
#' prior to RAVEL?.
#' @param WhiteStripe_Type What modality is used for WhiteStripe? Should be one
#' of T1, T2 or FLAIR.
#' @param k Number of unwanted factors to be included in the RAVEL model.
#' @param returnMatrix Should the matrix of normalized intensities be returned?
#' @param writeToDisk Should the normalized images be saved to disk?
#' @param verbose Should messages be printed?
#' @param ... Arguments to be passed to whitestripe
#' @return if \code{returnMatrix} is \code{FALSE}, no value returned, but
#' RAVEL-corrected images are saved. If \code{returnMatrix} is \code{TRUE},
#' RAVEL-corrected images are saved and a matrix of normalized intensities is
#' returned.
#' @author Jean-Philippe Fortin
#' @importFrom pbapply pboptions
#' @importFrom oro.nifti readNIfTI
#' @importFrom WhiteStripe whitestripe whitestripe_norm
#' @export
normalizeRAVEL <- function(input.files,
                           output.files = NULL,
                           brain.mask = NULL,
                           control.mask = NULL,
                           WhiteStripe = TRUE,
                           WhiteStripe_Type = c("T1", "T2", "FLAIR"),
                           k = 1,
                           returnMatrix = TRUE,
                           writeToDisk = FALSE,
                           verbose = TRUE, ...) {
  # RAVEL correction procedure:
  WhiteStripe_Type <- match.arg(WhiteStripe_Type)
  if (WhiteStripe_Type == "FLAIR")
    WhiteStripe_Type <- "T2"
  if (!verbose) {
    pboptions(type = "none")
  }
  
  if (!is.null(brain.mask)) {
    if (is.character(brain.mask)) {
      brain.mask <- readNIfTI(brain.mask, reorient = FALSE)
    }
    brain.indices <- brain.mask == 1
  } else {
    stop("brain.mask must be provided.")
  }
  
  if (is.null(output.files)) {
    output.files <- gsub(".nii.gz|.nii", "_RAVEL.nii.gz", input.files)
  }
  
  .ravel_correction <- function(V, Z) {
    means <- rowMeans(V)
    beta   <- solve(t(Z) %*% Z) %*% t(Z) %*% t(V)
    fitted <- t(Z %*% beta)
    res   <- V - fitted
    res   <- res + means
    res
  }
  
  if (verbose) {
    message("[normalizeRAVEL] Creating the voxel intensities matrix V. \n")
    if (WhiteStripe) {
      message(
        paste0(
          "[normalizeRAVEL] WhiteStripe intensity normalization",
          " is applied to each scan. \n"
        )
      )
    } else {
      message(
        paste0(
          "[normalizeRAVEL] WhiteStripe intensity normalization",
          " not applied.  \n"))
    }
  }
  # Matrix of voxel intensities:
  V <- pblapply(input.files, function(x) {
    brain <- readNIfTI(x, reorient = FALSE)
    if (WhiteStripe) {
      indices <-
        whitestripe(brain, type = WhiteStripe_Type, verbose = FALSE, ...)
      brain   <- whitestripe_norm(brain, indices$whitestripe.ind)
    }
    if (!is.null(brain.mask)) {
      brain <- as.vector(brain[brain.indices])
    }
    brain
  })
  V <- do.call(cbind, V)
  
  
  # Submatrix of control voxels:
  if (verbose)
    message("[normalizeRAVEL] Creating the control voxel matrix Vc. \n")
  if (class(control.mask) == "nifti") {
    control.indices <- control.mask == 1
  } else {
    control.indices <- readNIfTI(control.mask, reorient = FALSE) == 1
  }
  control.indices <- control.indices[brain.mask == 1]
  Vc <- V[control.indices, , drop = FALSE]
  
  
  
  
  if (verbose) {
    message("[normalizeRAVEL] Estimating the unwanted factors Z. \n")
  }
  Z  <- svd(Vc)$v[, 1:k, drop = FALSE] # Unwanted factors
  
  
  if (verbose) {
    message("[normalizeRAVEL] Performing RAVEL correction \n")
  }
  V.norm <- .ravel_correction(V, Z)
  
  
  if (writeToDisk) {
    if (verbose) {
      message("[normalizeRAVEL] Writing out the corrected images \n")
    }
    pblapply(1:ncol(V.norm), function(i) {
      .write_brain(
        brain.norm = V.norm[, i],
        output.file = output.files[i],
        brain.mask = brain.mask
      )
    })
  }
  if (returnMatrix) {
    return(V)
  }
}




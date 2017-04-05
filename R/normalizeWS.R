# Assuming images are registered and normalized beforehand


#' WhiteStripe intensity normalization.
#' 
#' WhiteStripe intensity normalization.
#' 
#' 
#' @param input.files Vector of filenames for the input images. Must be NIfTI
#' files.
#' @param output.files Optional vector of filenames for the output images. By
#' default, will be the \code{input.files} with suffix "WS".
#' @param brain.mask Filename for the brain binary mask specifying the template
#' space brain. Must be a NIfTI file.
#' @param WhiteStripe_Type What modality is used for WhiteStripe? Should be one
#' of T1, T2 or FLAIR.
#' @param writeToDisk Should the normalized scans be saved to the disk?
#' @param returnMatrix Should the matrix of normalized intensities be returned?
#' @param verbose Should messages be printed?
#' @return if \code{returnMatrix} is \code{FALSE}, no value returned, but
#' WhiteStripe-normalized images are saved. If \code{returnMatrix} is
#' \code{TRUE}, WhiteStripe-normalized images are saved and a matrix of
#' normalized intensities is returned.
#' @param ... additional arguments to pass to \code{\link{whitestripe}}
#' @author Jean-Philippe Fortin
#' @importFrom pbapply pboptions pblapply
#' @importFrom oro.nifti readNIfTI
#' @importFrom WhiteStripe whitestripe whitestripe_norm 
#' @export
normalizeWS <-
  function(input.files,
           output.files = NULL,
           brain.mask = NULL,
           WhiteStripe_Type = c("T1", "T2", "FLAIR"),
           writeToDisk = FALSE,
           returnMatrix = TRUE,
           verbose = TRUE,
           ...) {
    WhiteStripe_Type <- match.arg(WhiteStripe_Type)
    if (WhiteStripe_Type == "FLAIR") {
      WhiteStripe_Type <- "T2"
    }
    # RAVEL correction procedure:
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
      output.files <- gsub(".nii.gz|.nii", "_WS.nii.gz", input.files)
    }
    
    if (verbose) {
      message("[normalizeWS] WhiteStripe intensity normalization is applied to each scan. \n")
    }
    
    # Matrix of voxel intensities:
    V <- pblapply(input.files, function(x) {
      brain <- readNIfTI(x, reorient = FALSE)
      indices <- whitestripe(brain, type = WhiteStripe_Type, 
                             verbose = FALSE, ...)
      brain   <- whitestripe_norm(brain, indices$whitestripe.ind)
      brain <- as.vector(brain[brain.indices])
      brain
    })
    V <- do.call(cbind, V)
    
    if (writeToDisk) {
      if (verbose)
        cat("[normalizeWS] Writing out the corrected images \n")
      pblapply(1:ncol(V), function(i) {
        .write_brain(
          brain.norm = V[, i],
          output.file = output.files[i],
          brain.mask = brain.mask
        )
      })
    }
    if (returnMatrix) {
      return(V)
    }
  }

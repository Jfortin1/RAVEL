# Assuming images are registered and normalized beforehand


#' Whole-brain z-score intensity normalization
#'
#' Whole-brain z-score intensity normalization.
#'
#'
#' @param input.files Vector of filenames for the input images. Must be NIfTI
#' files.
#' @param output.files Optional vector of filenames for the output images. By
#' default, will be the \code{input.files} with suffix "ZScore".
#' @param brain.mask Filename for the brain binary mask specifying the template
#' space brain. Must be a NIfTI file.
#' @param writeToDisk Should the scans be saved to the disk? FALSE by default. 
#' @param returnMatrix Should the matrix of intensities be returned? FALSE by default.
#' @param verbose Should messages be printed?
#' @return if \code{returnMatrix} is \code{FALSE}, no value returned, but
#' Zscore-normalized images are saved. If \code{returnMatrix} is \code{TRUE},
#' Zscore-normalized images are saved and a matrix of normalized intensities is
#' returned.
#' @author Jean-Philippe Fortin
#' @importFrom pbapply pblapply
#' @importFrom stats sd
#' @export
normalizeZScore <- function(input.files,
           output.files = NULL,
           brain.mask = NULL,
           writeToDisk = FALSE,
           returnMatrix = TRUE,
           verbose = TRUE
){
    if (!verbose) {
      pboptions(type = "none")
    }

    if (!is.null(brain.mask)) {
      brain.mask = neurobase::check_nifti(brain.mask,
                                          reorient = FALSE,
                                          allow.array = FALSE)
      brain.indices <- brain.mask == 1
    } else {
      stop("brain.mask must be provided.")
    }

    if (verbose) {
      message("[normalizeZScore] Creating the voxel intensities matrix V. \n")
    }

    # Matrix of voxel intensities:
    V <- pblapply(input.files, function(x) {
      brain <- check_nifti(x, reorient = FALSE, allow.array = FALSE)
      if (!is.null(brain.mask)) {
        brain <- as.vector(brain[brain.indices])
      }
      brain
    })
    V <- do.call(cbind, V)
    # Whole brain z-score normalization:
    means <- apply(V,2,mean,na.rm=TRUE)
    sds <- apply(V,2,sd, na.rm=TRUE)
    V <- sweep(V,2,means, "-")
    V <- sweep(V,2,sds, "/")

    input.files = checkimg(input.files)
    if (is.null(output.files)) {
      output.files <- gsub(".nii.gz|.nii", "_ZScore.nii.gz", input.files)
    }

    if (writeToDisk) {
      if (verbose) {
        message("[normalizeZScore] Writing out the corrected images \n")
      }
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

# Assuming images are registered and normalized beforehand


#' Create matrix of voxel intensities without normalization.
#'
#' Create matrix of voxel intensities without normalization.
#'
#'
#' @param input.files Vector of filenames for the input images. Must be NIfTI
#' files.
#' @param output.files Optional vector of filenames for the output images. By
#' default, will be the \code{input.files} with suffix "WS".
#' @param brain.mask Filename for the brain binary mask specifying the template
#' space brain. Must be a NIfTI file.
#' @param writeToDisk Should the scans be saved to the disk?
#' @param returnMatrix Should the matrix of intensities be returned?
#' @param verbose Should messages be printed?
#' @return if \code{returnMatrix} is \code{FALSE}, no value returned, but
#' images are saved. If \code{returnMatrix} is \code{TRUE}, images are saved
#' and a matrix of intensities is returned.
#' @author Jean-Philippe Fortin
#' @importFrom pbapply pblapply
#' @export
normalizeRaw <-
  function(input.files,
           output.files = NULL,
           brain.mask = NULL,
           writeToDisk = FALSE,
           returnMatrix = TRUE,
           verbose = TRUE) {
    # RAVEL correction procedure:
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
      message("[normalizeRaw] Creating the voxel intensities matrix V. \n")
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

    input.files = checkimg(input.files)
    if (is.null(output.files)) {
      output.files <- gsub(".nii.gz|.nii", "_RAW.nii.gz", input.files)
    }

    if (writeToDisk) {
      if (verbose) {
        message("[normalizeRaw] Writing out the corrected images \n")
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

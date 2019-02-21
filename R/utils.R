# Create intersection mask:


#' Intersection of a list of brain masks
#' 
#' Intersection of a list of brain masks.
#' 
#' 
#' @param list List of the mask input files (nifti files) or mask nifti
#' objects.
#' @param output.file Optional filename for the output image.
#' @param prob Probability of intersection. Ex: 0.9 will return a mask for
#' voxels that intersect 90\% of the time.
#' @param reorient Should the input images be reoriented? FALSE by default.
#' @param returnObject Should a nifti object be returned?
#' @param writeToDisk Should the intersection mask be saved to the disk?
#' @param verbose Should messages be printed?
#' @return If \code{returnObject} is \code{TRUE}, a \code{nifti} object is
#' returned, and if \code{writeToDisk} is \code{TRUE}, the intersection mask is
#' saved to disk as a NIfTI file.
#' @author Jean-Philippe Fortin
#' @examples
#' 
#' if (require(EveTemplate)){
#'   path1 <- path2 <- getEvePath()
#'   inter <- maskIntersect(list(path1,path2), writeToDisk=FALSE)
#' }
#' @importFrom oro.nifti writeNIfTI
#' @export
maskIntersect <- function(list,
                          output.file = NULL,
                          prob = 1,
                          reorient = FALSE,
                          returnObject = TRUE,
                          writeToDisk = TRUE,
                          verbose = TRUE) {
  if (!verbose) {
    pboptions(type = "none")
  }
  
  # Checks:
  if (is.atomic(list)) {
    list <- as.list(list)
  }
  if (writeToDisk & is.null(output.file)) {
    stop("output.file must be specified if writeToDisk is true.")
  }
  
  n <- length(list)
  inter  <- list[[1]]
  if (class(inter) == "nifti") {
    inter <- Reduce("+", list)
  } else if (class(inter) == "character") {
    inter <- readNIfTI(list[[1]], reorient = reorient)
    for (i in 2:n) {
      inter <- inter + readNIfTI(list[[i]], reorient = reorient)
    }
  } else {
    stop("list must be either a list of nifti objects or a list of NIfTI file paths.")
  }
  
  # Creating the intersection map:
  cutoff <- floor(prob * n)
  inter[inter < cutoff]  <- 0
  inter[inter >= cutoff] <- 1
  
  # Writing to disk:
  if (writeToDisk) {
    filename <- gsub(".nii.gz|.nii", "", output.file)
    writeNIfTI(inter, filename)
  }
  
  # Returning object:
  if (returnObject) {
    return(inter)
  }
}

.write_brain <- function(brain.norm, output.file, brain.mask) {
  if (is.character(brain.mask)) {
    brain.mask <- readNIfTI(brain.mask, reorient = FALSE)
  }
  brain.mask[brain.mask == 1] <- brain.norm
  output.file <- gsub(".nii.gz|.nii", "", output.file)
  writenii(nim=brain.mask, filename=output.file)
}

# # Function from the extrantsr package (thx to John Muschelli):
# tempants <- function(x, gzipped = TRUE){
#   if (inherits(x, "character")) {
#     return(x)
#   } else {
#     if (inherits(x, "antsImage")) {
#       ext = ".nii"
#       if (gzipped) ext = ".nii.gz"
#       tfile = paste0(tempfile(), ext)
#       antsImageWrite(x, tfile)
#       return(tfile)
#     }  else {
#       stop("x has unknown class - not char or nifti")
#     }
#   }
#   return(FALSE)
# }

# ants2oro <- function(img, reorient = FALSE){
#   if ( is.antsImage(img) | is.character(img) ) {
#     fname = tempants(img)
#     img = readNIfTI(fname, reorient = reorient)
#     return(img)
#   }
#   if ( is.nifti(img) ) {
#     return(img)
#   }
#   stop("img not class nifti or antsImage")
#   return(NULL)
# }

# oro2ants <- function(img, reference = NULL){
#   if (!is.null(reference)) {
#     if (is.antsImage(reference)) {
#       img = as(img, Class = "array")
#       aimg = as.antsImage(img)
#       aimg = antsCopyImageInfo(
#         target = aimg,
#         reference = reference)
#       return(aimg)
#     }
#   }
#   if (  is.nifti(img) | is.character(img) ) {
#     fname = checkimg(img)
#     stopifnot(file.exists(fname))
#     img = antsImageRead(fname)
#     return(img)
#   }
#   if ( is.antsImage(img) ) {
#     return(img)
#   }
#   stop("img not class nifti or antsImage")
#   return(NULL)
# }


# # To clean registration files:
# .clean_reg_files <- function(outprefix){
#     file.remove(paste0(outprefix, "1InverseWarp.nii.gz"))
#     file.remove(paste0(outprefix, "1Warp.nii.gz"))
#     file.remove(paste0(outprefix, "0GenericAffine.mat"))
# }

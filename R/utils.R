# Create intersection mask:
mask_intersect <- function(list){
	n <- length(list)
	inter <- list[[1]]
	for (i in 2:n){
		inter[!(inter==1 & list[[i]]==1)] <- 0
		inter
	}
	inter
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

# Function from the extrantsr package (thx to John Muschelli): 
tempants <- function(x, gzipped = TRUE){
  if (inherits(x, "character")) {
    return(x)
  } else {
    if (inherits(x, "antsImage")) {
      ext = ".nii"
      if (gzipped) ext = ".nii.gz"
      tfile = paste0(tempfile(), ext)
      antsImageWrite(x, tfile)
      return(tfile)
    }  else {
      stop("x has unknown class - not char or nifti")
    }
  }
  return(FALSE)
}

ants2oro <- function(img, reorient = FALSE){
  if ( is.antsImage(img) | is.character(img) ) {
    fname = tempants(img)
    img = readNIfTI(fname, reorient = reorient)
    return(img)
  }
  if ( is.nifti(img) ) {
    return(img)
  }
  stop("img not class nifti or antsImage")
  return(NULL)
}

oro2ants <- function(img, reference = NULL){
  if (!is.null(reference)) {
    if (is.antsImage(reference)) {
      img = as(img, Class = "array")
      aimg = as.antsImage(img)
      aimg = antsCopyImageInfo(
        target = aimg, 
        reference = reference)
      return(aimg)
    }
  }
  if (  is.nifti(img) | is.character(img) ) {
    fname = checkimg(img)
    stopifnot(file.exists(fname))
    img = antsImageRead(fname)
    return(img)
  }
  if ( is.antsImage(img) ) {
    return(img)    
  }   
  stop("img not class nifti or antsImage")
  return(NULL)
}


# To clean registration files:
clean_reg_files <- function(outprefix){
    file.remove(paste0(outprefix, "1InverseWarp.nii.gz"))
    file.remove(paste0(outprefix, "1Warp.nii.gz"))
    file.remove(paste0(outprefix, "0GenericAffine.mat"))
}



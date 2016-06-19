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
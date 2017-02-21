# Original code version: Taki Shinohara - September 10, 2012
# Current code version: Jean-Philippe Fortin - July 29 2015

# Assuming images are registered and normalized beforehand
normalizeHM <- function(input.files, output.files=NULL, brain.mask=NULL, 
	type=c("T1", "T2", "FLAIR", "PD"), writeToDisk=FALSE, returnMatrix=TRUE, verbose=TRUE){
	
	type <- match.arg(type)
	# RAVEL correction procedure:
	if (!verbose) pboptions(type="none") 

	if (!is.null(brain.mask)){
		if (is(brain.mask, "character")){
			brain.mask <- readNIfTI(brain.mask, reorient=FALSE)
		}
		brain.indices <- brain.mask==1
	} else {
		stop("brain.mask must be provided.")
	}

	if (is.null(output.files)){
		output.files <- gsub(".nii.gz|.nii","_HM.nii.gz", input.files)
	}
	
    cat("[normalizeHM] Histogram matching is applied to each scan. \n")
	# Matrix of voxel intensities:
	V <- pblapply(input.files, function(x){
		brain  <- readNIfTI(x, reorient=FALSE)
        brain  <- .hm(brain, type=type)
        brain  <- as.vector(brain[brain.indices])
		brain
	})
	V <- do.call(cbind, V)

	if (writeToDisk){
		if (verbose) cat("[normalizeHM] Writing out the corrected images \n")
		pblapply(1:ncol(V), function(i){
			.write_brain(brain.norm = V[,i], output.file = output.files[i], brain.mask=brain.mask)
		})
	} 
	if (returnMatrix){
		return(V)
	}	
}



# Per subject histogram matching function:
# Assuming img is a nifti object:
.hm <- function(img, type=c("T1", "T2", "FLAIR", "PD")){
	type <- match.arg(type)
	#seed <- 123413
	#set.seed(seed)
	i.min   <- 0.01 # Trimming Quantiles
	i.max   <- 0.99
	i.s.min <- 0 # Output Range
	i.s.max <- 1
	h <- seq(0.1, 0.9, by=0.1) # Quantile Landmarks

	# This function does training histogram normalization as 
	# in Shah et al. (2011) and Nyul et al. (2000)
	#get.landmarks <- function(rawdata, i.min, i.max, i.s.min, i.s.max, h, mask) {
	#	scaled.data <- (rawdata-quantile(rawdata[mask>0],i.min)+i.s.min)/quantile(rawdata[mask>0],i.max)*i.s.max
	#	return(quantile(scaled.data[mask>0], probs=h))
	#}

	#This function does histogram normalization transformation as in Shah et al. (2011) and Nyul et al. (2000)
	do.hist.norm <- function(rawdata, i.min, i.max, i.s.min, i.s.max, h, m, mask) {

		m.obs <- quantile(rawdata[mask>0], probs = c(i.min, h, i.max))
		m.withends <- c(i.s.min, m, i.s.max)
		transformed.data <- rawdata

		transformed.data[transformed.data <= quantile(rawdata[mask>0], probs=i.min)] <- i.s.min
		transformed.data[transformed.data >= quantile(rawdata[mask>0], probs=i.max)] <- i.s.max

		for (hist.section.i in 1:(length(h)+1)) {
			which.data <- (rawdata[mask>0]<m.obs[hist.section.i+1])&(rawdata[mask>0]>=m.obs[hist.section.i])
			transformed.data[mask>0][which.data]<-(rawdata[mask>0][which.data]-m.obs[hist.section.i])/(m.obs[hist.section.i+1]-m.obs[hist.section.i])*(m.withends[hist.section.i+1]-m.withends[hist.section.i])+m.withends[hist.section.i]
		}
		return(transformed.data)
	}


	data(histogram_match_train_healthy)
	if (type=="T1"){
		land.m 	 <- apply(t1.landmarks, 2, mean)
	} else if (type=="T2"){
		land.m 	 <- apply(t2.landmarks, 2, mean)
	} else if (type=="FLAIR"){
		land.m 	 <- apply(flair.landmarks, 2, mean)
	} else if (type=="PD"){
		land.m 	 <- apply(pd.landmarks, 2, mean)
	}
	img.fg  <- 1*(img>mean(img))
	img <- do.hist.norm(img, i.min, i.max, i.s.min, i.s.max, h, land.m, img.fg)
	return(img)
}

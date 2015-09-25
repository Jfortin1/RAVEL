# Histogram matching function:
hm <- function(input, output, what=c("T1", "FLAIR")){
	what <- match.arg(what)
	img  <- readNIfTI(input)
	seed <- 123413
	set.seed(seed)
	i.min   <- 0.01 # Trimming Quantiles
	i.max   <- 0.99
	i.s.min <- 0 # Output Range
	i.s.max <- 1
	h <- seq(0.1, 0.9, by=0.1) # Quantile Landmarks

	# This function does training histogram normalization as 
	# in Shah et al. (2011) and Nyul et al. (2000)
	get.landmarks<-function(rawdata, i.min, i.max, i.s.min, i.s.max, h, mask) {
		scaled.data <- (rawdata-quantile(rawdata[mask>0],i.min)+i.s.min)/quantile(rawdata[mask>0],i.max)*i.s.max
		return(quantile(scaled.data[mask>0],probs=h))
	}

	#This function does histogram normalization transformation as in Shah et al. (2011) and Nyul et al. (2000)
	do.hist.norm<-function(rawdata, i.min, i.max, i.s.min, i.s.max, h, m, mask) {

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

	load("../../objects/histogram_match_train_healthy.rda")
	if (what=="T1"){
		img.m 	<- apply(t1.landmarks, 2, mean)
	} else {
		img.m 	<- apply(flair.landmarks, 2, mean)
	}
	
	img.fg   <- 1*(img>mean(img))
	img <- do.hist.norm(img, i.min, i.max, i.s.min, i.s.max, h, img.m, img.fg)
	writeNIfTI(img, output)
}





























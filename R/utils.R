
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





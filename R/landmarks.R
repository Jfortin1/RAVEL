#' @title Estimated Landmarks for Various Imaging Modalities
#'
#' @description A list of histogram landmarks for histogram matching.
#'
#' @format A list with 4 elements, which are \code{flair}, \code{pd}, \code{t1}, 
#' \code{t2}.  Each of these elements contains a matrix of 9 columns, each different
#' quantiles of the distribution of voxel intensites and 33 rows, corresponding to
#' healthy individuals
"landmarks"

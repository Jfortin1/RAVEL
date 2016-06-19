library(oro.nifti)
library(Biobase)
library(pbapply)
library(S4Vectors)
library(matrixStats)
library(fslr)

dir <- "/Users/Jean-Philippe/RAVELData/data"
template.file <- file.path(dir,"JHU_MNI_SS_T1.nii.gz")
brainMask.file <- file.path(dir,"JHU_MNI_SS_T1_Brain_Mask.nii.gz")
input.files <- file.path(dir,paste0("scan", 1:3, "_processed.nii.gz"))


brainMask <- readNIfTI(brainMask_path, reorient=FALSE)
brain.indices <- brainMask==1
template <- readNIfTI(template_path, reorient=FALSE)
template[brainMask==0] <- 0

# Matrix of voxel intensities:
V <- pblapply(input.files, function(x){
    brain <- readNIfTI(x, reorient=FALSE)
    brain <- as.vector(brain[brain.indices])
    brain
})
V <- do.call(cbind, V)

 

setClass("mriSet",
         contains = "eSet",
         slots  = c(template="nifti")
         )


object <- obj
setValidity("mriSet", function(object) {
    msg <- NULL
    n.rows <- dim(object)[1]
    n.cols <- dim(object)[2]
    if (sum(template!=0)!=n.rows){
        msg <- "template is not compatible with the intensities matrix."
    }
    msg <- validMsg(NULL, assayDataValidMembers(assayData(object),
                                                c("intensities")))

    if (is.null(msg)) TRUE else msg
})



basic <- function(input.files, template.file, brainMask.file){
    brainMask <- readNIfTI(brainMask.file, reorient=FALSE)
    brain.indices <- brainMask==1
    template <- readNIfTI(template.file, reorient=FALSE)
    template[brainMask==0] <- 0

    # Matrix of voxel intensities:
    V <- pblapply(input.files, function(x){
        brain <- readNIfTI(x, reorient=FALSE)
        brain <- as.vector(brain[brain.indices])
        brain
    })
    V <- do.call(cbind, V)
    feature.names <- paste0("Voxel_",1:nrow(V))
    sample.names <- gsub(".nii.gz","",basename(input.files))
    rownames(V) <- feature.names
    colnames(V) <- sample.names

    obj <- mriSet(intensities=V, template=template)
    
    # Phenotype data:
    sampleNames(obj) <- sample.names
    temp <- data.frame(files=input.files)
    rownames(temp) <- sample.names
    pData(obj) <- temp

    # Feature data
    featureNames(obj) <- feature.names
    temp <- data.frame(template_indices=which(brain.indices))
    rownames(temp) <- feature.names
    featureData(obj) <-  temp
    obj
}




a <- basic(input.files, template.file, brainMask.file)


mriSet <- function(intensities = new("matrix"), template=new("nifti"), ...){
    mriSet <- new("mriSet", intensities = intensities, template=template, ...)
    mriSet
}





.show.mriSet <- function(object) {
    cat(paste0("An object of class ", class(object), " with the following data: \n \n") )
    cat("Imaging data:", paste(dim(object)[[1]], "features (voxels),", dim(object)[[2]], "samples (scans)"), "\n")
    cat("  feature names:",
        paste(assayDataElementNames(object), collapse=", "), "\n")
    cat("Phenotype data: ")
    show(phenoData(object))
    cat("Feature data: ")
    show(featureData(object))
}


setMethod("show", "mriSet", function(object) {
    .show.mriSet(object)
    #.show.annotation(annotation(object))
})






setGeneric("getIntensities", function(object) standardGeneric("getIntensities"))
setGeneric("getTemplate", function(object) standardGeneric("getTemplate"))
#setGeneric("getBrainMask", function(object) standardGeneric("getBrainMask"))
setGeneric("scanNames", function(object) standardGeneric("scanNames"))
setGeneric("scanNames<-", function(object, value) standardGeneric("scanNames<-"))
setGeneric("intensities<-", function(object, value) standardGeneric("intensities<-"))
setGeneric("intensities", function(object) standardGeneric("intensities"))





setMethod("getIntensities", signature(object = "mriSet"),
          function(object) {
              assayDataElement(object, "intensities")
          })
setMethod("intensities", signature(object = "mriSet"),
          function(object) {
              assayDataElement(object, "intensities")
          })

setMethod("intensities<-", signature(object = "mriSet", value="matrix"),
          function(object, value) {
              assayDataElement(object, "intensities") <- value
              object
          })



setMethod("getTemplate", signature(object = "mriSet"),
          function(object) {
              object@template
          })
#setMethod("getBrainMask", signature(object = "mriSet"),
#          function(object) {
#              object@brainMask
#          })



setReplaceMethod("pData", signature(object = "mriSet", value = "DataFrame"),
                 function(object, value) {
                     df <- as.data.frame(value)
                     pData(object) <- df
                     object
                 })
setMethod("scanNames<-", signature(object = "mriSet", value = "vector"),
                 function(object, value) {
                     sampleNames(object) <- value
                     object
                 })

setMethod("scanNames<-", signature(object = "mriSet", value = "vector"),
                 function(object, value) {
                     sampleNames(object) <- value
                     object
                 })



setMethod("featureData<-", signature(object = "mriSet", value = "data.frame"),
                 function(object, value) {
                     featureData(object) <- new("AnnotatedDataFrame",value)
                     object
                 })

getMeanVar <- function(mriSet, ...){
    means <- rowMeans(getIntensities(mriSet), ...)
    vars <- rowVars(getIntensities(mriSet), ...)
    cbind(means=means, vars=vars)
}


plotMeanVar <- function(mriSet, na.rm=TRUE, how.many=10000, sqrt=TRUE, ...){
    temp <- getMeanVar(mriSet, na.rm=na.rm)
    n <- nrow(temp)
    if (n > how.many){
        indices <- sample(n, how.many)
    } else {
        indices <- 1:n
    }
    x <- temp[indices,1]
    y <- temp[indices,2]
    if (sqrt) y <- sqrt(y)
    plot(x,y, ...)
}

getFeatureData <- function(mriSet){
    featureData(mriSet)@data
}


addTemplateSeg <- function(mriSet, seg.alg = "FAST", ...){
    require(fslr)
    temp <- getTemplate(mriSet)
    seg  <- fast(temp, ...)
    seg  <- seg[template!=0]
    pp <- pData(mriSet)
    pData(mriSet) <- data.frame(tissue_FAST=seg, pData(mriSet))
}







# Testing:
seg <- fast(template)



plotMeanVar(obj)

obj <- mriSet(intensities=V, template=template)
#pData(obj) <- data.frame(sampleNames=c(1,2,4), age=c(1,5,2))
template <- getTemplate(obj)
brainMask <- getBrainMask(obj)








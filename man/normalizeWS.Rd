\name{normalizeWS}
\alias{normalizeWS}
\title{
WhiteStripe intensity normalization. 
}
\description{
 WhiteStripe intensity normalization. 
}
\usage{
normalizeWS(input.files, output.files=NULL, brain.mask=NULL, 
  WhiteStripe_Type=c("T1","T2","FLAIR"), writeToDisk=FALSE, returnMatrix=FALSE, verbose=TRUE)
}
\arguments{
  \item{input.files}{
    Vector of filenames for the input images. Must be NIfTI files. 
  }
  \item{output.files}{
    Optional vector of filenames for the output images. By default, will be the \code{input.files} with suffix "WS".
  }
  \item{brain.mask}{
    Filename for the brain binary mask specifying the template space brain. Must be a NIfTI file. 
  }
  \item{WhiteStripe_Type}{
    What modality is used for WhiteStripe? Should be one of T1, T2 or FLAIR.
  }
  \item{writeToDisk}{
    Should the normalized scans be saved to the disk? 
  }
  \item{returnMatrix}{
    Should the matrix of normalized intensities be returned?
  }
  \item{verbose}{
    Should messages be printed?
  }
}

\value{
if \code{returnMatrix} is \code{FALSE}, no value returned, but WhiteStripe-normalized images are saved. If \code{returnMatrix} is \code{TRUE}, WhiteStripe-normalized images are saved and a matrix of normalized intensities is returned. 
}
\author{
Jean-Philippe Fortin
}
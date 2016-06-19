\name{preprocessWS}
\alias{preprocessWS}
\title{
WhiteStripe intensity normalization. 
}
\description{
 WS stuff. 
}
\usage{
preprocessWS(input.files, output.files=NULL, brain.mask=NULL, 
  WhiteStripe_Type="T1", verbose=TRUE, writeToDisk=FALSE, returnMatrix=FALSE)
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
    What modality is used for WhiteStripe? For the moment, only "T1" is supported. 
  }
  \item{verbose}{
    Should messages be printed?
  }
  \item{writeToDisk}{
    Should the normalized scans be saved to the disk? 
  }
  \item{returnMatrix}{
    Should the matrix of normalized intensities be returned?
  }
}

\value{
if \code{returnMatrix} is \code{FALSE}, no value returned, but WhiteStripe-normalized images are saved. If \code{returnMatrix} is \code{TRUE}, WhiteStripe-normalized images are saved and a matrix of normalized intensities is returned. 
}
\author{
Jean-Philippe Fortin
}
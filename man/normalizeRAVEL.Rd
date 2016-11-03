\name{normalizeRAVEL}
\alias{normalizeRAVEL}
\title{
RAVEL intensity normalization
}
\description{
 RAVEl intensity normalization.  
}
\usage{
normalizeRAVEL(input.files, output.files=NULL, brain.mask=NULL, 
  control.mask=NULL, WhiteStripe=TRUE, 
  WhiteStripe_Type=c("T1","T2", "FLAIR"),  
  k=1, returnMatrix=TRUE, writeToDisk=FALSE, verbose=TRUE)
}
\arguments{
  \item{input.files}{
    Vector of filenames for the input images. Must be NIfTI files. 
  }
  \item{output.files}{
    Optional vector of filenames for the output images. 
    By default, will be the \code{input.files} with suffix "RAVEL".
  }
  \item{brain.mask}{
    Filename for the brain binary mask specifying the template space brain.
    Must be a NIfTI file. 
  }
  \item{control.mask}{
    Filename for the control region binary mask to be used in RAVEL. 
    Must be a NIfTI file.
  }
  \item{WhiteStripe}{
    Should White Stripe intensity normalization be performed prior to RAVEL?. 
  }
  \item{WhiteStripe_Type}{
    What modality is used for WhiteStripe? Should be one of T1, T2 or FLAIR.
  }
  \item{k}{
    Number of unwanted factors to be included in the RAVEL model. 
  }
  \item{returnMatrix}{
    Should the matrix of normalized intensities be returned?
  }
  \item{writeToDisk}{
    Something
  }
  \item{verbose}{
    Should messages be printed?
  }
}

\value{
if \code{returnMatrix} is \code{FALSE}, no value returned, but RAVEL-corrected 
images are saved. If \code{returnMatrix} is \code{TRUE}, RAVEL-corrected images 
are saved and a matrix of normalized intensities is returned. 
}
\author{
Jean-Philippe Fortin
}
\name{normalizeRaw}
\alias{normalizeRaw}
\title{
Create matrix of voxel intensities without normalization.
}
\description{
Create matrix of voxel intensities without normalization.
}
\usage{
normalizeRaw(input.files, output.files=NULL, brain.mask=NULL, verbose=TRUE, writeToDisk=FALSE, returnMatrix=TRUE)
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
  \item{verbose}{
    Should messages be printed?
  }
  \item{writeToDisk}{
    Should the scans be saved to the disk? 
  }
  \item{returnMatrix}{
    Should the matrix of intensities be returned?
  }
}

\value{
if \code{returnMatrix} is \code{FALSE}, no value returned, but images are saved. If \code{returnMatrix} is \code{TRUE}, images are saved and a matrix of intensities is returned. 
}
\author{
Jean-Philippe Fortin
}
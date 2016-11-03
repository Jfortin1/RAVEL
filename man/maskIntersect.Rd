\name{maskIntersect}
\alias{maskIntersect}
\title{
Intersection of a list of brain masks
}
\description{
Intersection of a list of brain masks. 
}
\usage{
maskIntersect(list, output.file = NULL, prob=1, reorient = FALSE, returnObject = TRUE,
	writeToDisk = TRUE, verbose = TRUE)
}
\arguments{
  \item{list}{
    List of the mask input files (nifti files) or mask nifti objects. 
  }
  \item{output.file}{
    Optional filename for the output image. 
  }
  \item{prob}{
    Probability of intersection. Ex: 0.9 will return a mask for voxels that intersect 90\% of the time. 
  }
  \item{reorient}{
    Should the input images be reoriented? FALSE by default. 
  }
  \item{returnObject}{
    Should a nifti object be returned?
  }
  \item{writeToDisk}{
    Should the intersection mask be saved to the disk? 
  }
  \item{verbose}{
    Should messages be printed?
  }
}
\value{
If \code{returnObject} is \code{TRUE}, a \code{nifti} object is returned, and if \code{writeToDisk} is \code{TRUE}, the intersection mask is saved to disk as a NIfTI file. 
}
\author{
Jean-Philippe Fortin
}
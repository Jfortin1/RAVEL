\name{oro2ants}
\alias{oro2ants}
\title{Convert between nifti and antsImage}
\description{
NIfTI data can be converted between \code{nifti} 
(from the oro.nifti package) and \code{antsImage} objects.
}
\usage{
oro2ants(img, reference=NULL)
}
\arguments{
\item{img}{Object of class \code{nifti}}
\item{reference}{Object of class \code{antsImage} to 
copy header information (origin, spacing, direction) (experimental)}
}
\value{
Object of class \code{antsImage}
}
\author{
John Muschelli
}
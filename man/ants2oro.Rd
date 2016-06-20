\name{ants2oro}
\alias{ants2oro}
\title{
Convert between antsImage and nifti
}
\description{
NIfTI data can be converted between \code{antsImage} 
(from the ANTsR package) and nifti S4 objects.
}
\usage{
ants2oro(img, reorient = FALSE)
}
\arguments{
\item{img}{Object of class \code{antsImage}}
\item{reorient}{Reorientation passed to \code{\link{readnii}}}
}

\value{
Object of class \code{nifti}
}
\author{
John Muschelli
}
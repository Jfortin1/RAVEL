\name{tempants}
\alias{tempants}
\title{
Create temporary nii.gz file from antsImage
}
\description{
Takes in an object of class \code{antsImage}, writes it
to a temp file, and appends .nii or .nii.gz
}
\usage{
tempants(x, gzipped=TRUE)
}
\arguments{
\item{x}{Object of class \code{antsImage}}
\item{gzipped}{Should tempfile be gzipped or not?}
}

\value{
Character value of filename
}
\author{
John Muschelli
}
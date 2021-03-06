% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotDNA.R
\name{indexToRange}
\alias{indexToRange}
\title{Find contiguous ranges within a vector of indices}
\usage{
indexToRange(index)
}
\arguments{
\item{index}{Indices to be condensed e.g. from \code{\link{which}} or \code{\link{grep}}}
}
\value{
A data frame with columns: 
     \describe{
       \item{start:}{Start of a contiguous region}
       \item{end:}{End of a contiguous region}
     }
}
\description{
Take a vector of indices and returns ranges of contiguous regions.
}
\examples{
indexToRange(c(1:10,11,14,16,17:20))
}

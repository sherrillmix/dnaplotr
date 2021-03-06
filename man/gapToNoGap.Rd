% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotDNA.R
\name{gapToNoGap}
\alias{gapToNoGap}
\title{Convenience function to convert gapped coordinates to what the coordinates would be without gaps}
\usage{
gapToNoGap(refSeq, coords, gapChars = c("*", ".", "-"))
}
\arguments{
\item{refSeq}{the sequence containing gaps}

\item{coords}{coordinates on the gapped refSeq to be converted into equivalent nongap coordinates}

\item{gapChars}{characters interpreted as gaps}
}
\value{
A vector of coordinates in the gap free reference sequence
}
\description{
Take a reference sequence with gaps and coordinates in that gapped reference
sequence and translate the coordinates to the corresponding positions in the
gap free sequence. Note that positions landing on a gap are converted to the
first 5' nongap position e.g. the fifth position of 'G-----G' is converted
to 1 and the second position of '--G--GG' is converted to 0
}
\examples{
gapToNoGap('AA--AA-A',c(1:8))
}

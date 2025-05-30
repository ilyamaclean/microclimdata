% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/climate.R
\name{haduk_process}
\alias{haduk_process}
\title{Process UK Met Office climate data}
\usage{
haduk_process(r, tme, pathin)
}
\arguments{
\item{r}{a SpatRaster object used as a template region for which to extract data
from downloaded files.}

\item{tme}{a POSIXlt object giving the time period for which to download data}

\item{pathin}{directory to which data are saved using `haduk_download`.}
}
\value{
a list of wrapped SpatRasters of the following:
\describe{
 \item{precipitation}{daily precipitation (mm)}
 \item{tasmax}{Daily maximum temperature (deg C)}
 \item{tasmin}{Daily minimum temperature (deg C)}
}
}
\description{
This function processes daily 1km HadUK-Grid historic temperature
}
\seealso{
[haduk_download()], [haduk_blend()]
}

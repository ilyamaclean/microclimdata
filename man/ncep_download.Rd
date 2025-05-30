% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/climate.R
\name{ncep_download}
\alias{ncep_download}
\title{Downloads NCEP-NCAR Reanalysis climate data}
\usage{
ncep_download(r, tme)
}
\arguments{
\item{r}{a SpatRaster giving a template area for which to download data}

\item{tme}{a POSIXlt object giving the time period for which to download data}
}
\value{
a list of the following:
\describe{
 \item{tme}{POSIXlt object of period covered by download - typically longer than `tme`}
 \item{clim data}{a list of wrapped multi=layer SpatRasters of each climate variable}
 \item{tmin.2m}{Minimum temperature (K)}
 \item{tmax.2m}{Maximum temperature (K)}
 \item{shum.2m}{Specific humidity kg/kg}
 \item{pres.sfc}{Sea-level pressure Pa}
 \item{dswrf.sfc}{Downward shortwave radiation (W/m^2)}
 \item{dlwrf.sfc}{Downward longwave radiation (W/m^2)}
 \item{uwnd.10m}{u vector of wind speed at 10m (m/s)}
 \item{vwnd.10m}{v vector of wind speed at 10m (m/s)}
 \item{prate.sfc}{precipitation rate kg/s}
}
}
\description{
This function downloads six-hourly historic climate data from the NCEP-NCAR reanalysis project
}
\seealso{
[ncep_process()]
}

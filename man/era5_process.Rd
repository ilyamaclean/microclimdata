% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/climate.R
\name{era5_process}
\alias{era5_process}
\title{Processes ERA5 climate data}
\usage{
era5_process(
  req,
  pathin,
  r,
  tme,
  out = "grid",
  lat = NA,
  long = NA,
  resampleout = FALSE
)
}
\arguments{
\item{req}{a list of the request details sent to climateData Store as returned
by `era5_download`. If set to NA, all files in directory `pathin` are read in.}

\item{pathin}{directory in which files downloaded by `era5_download` are saved.}

\item{r}{a SpatRaster object used as a template for cropping and reprojecting the
outputs and resampling and masking outputs if specified.}

\item{out}{one of `grid` or `point` specifying whether to return a list of
arrays of climate data or a data.frame for a single point location.}

\item{lat}{latitude (decimal degrees) - alternative to supplying a template
raster if data are required for a point location.}

\item{resampleout}{optional logical indicating whether to resample and mask outputs
to match `r`.}

\item{lon}{longitude (decimal degrees) - alternative to supplying a template
raster if data are required for a point location.}
}
\value{
if `out = "point"` then a a data.frame of hourly weather comprising
the following columns:
\describe{
 \item{obs_time}{POSIXlt object of dates and times}
 \item{temp}{Temperature (deg C)}
 \item{relhum}{Relative humidity (percentage)}
 \item{pres}{Sea-level pressure (kPa)}
 \item{swdown}{Downward shortwave radiation on the horizontal (W/m^2)}
 \item{difrad}{Downward diffuse radiation (W/m^2)}
 \item{lwdown}{Downward longwave radiation W/m^2}
 \item{windspeed}{wind speed at 2m (m/s)}
 \item{winddir}{wind direction (decimal degrees)}
 \item{precip}{precipitation mm}
}

If `out = "grid"`a list of wrapped SpatRasters of the above climate variables.
}
\description{
This function processes hourly ERA5 climate data for use with `microgrid`
or `micropoint`
}
\seealso{
[era5_download()]
}

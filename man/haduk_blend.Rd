% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/climate.R
\name{haduk_blend}
\alias{haduk_blend}
\title{Combine coarse-resolution and UK Met Office data}
\usage{
haduk_blend(coarsedata, hadukdata)
}
\arguments{
\item{coarsedata}{coarse-resolution hourly climate data as returned by `ncep_process`
or `era5_process`.}

\item{hadukdata}{1km daily precipitation and temperature data as returned by
`haduk_process`}
}
\value{
a list comprising `climarrayr` - wrapped SpatRasters of 1km hourly
climate data and `tme` - a POSIXlt object of dates and times corresponding to
data in `climarrayr`. The dataset `climarrayr` is itself a list of wrapped
SpatRasters of the following variables:
\describe{
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
}
\description{
This function blends coarse resolution and daily 1km HadUK-Grid
historic temperature  and precipitation data to provide hourly 1km climate data
}
\seealso{
[haduk_download()], [haduk_process()]
}

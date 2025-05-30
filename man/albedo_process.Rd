% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ground.R
\name{albedo_process}
\alias{albedo_process}
\title{Processes downloaded MODIS albedo data}
\usage{
albedo_process(r, pathin)
}
\arguments{
\item{r}{a SpatRaster with a coordinate reference system defined. Used
to the location for which to crop downloaded data too.}

\item{pathin}{directory to which data are downloaded. Same as `pathout` in
[albedo_download()]}
}
\value{
a SpatRaster of white-sky albedo
}
\description{
Processes 500m grid resolution MODIS albedo data downlaoded using
[albedo_download()] and returns a single terra::SpatRaster for the area
specified by `r`
}
\details{
Downloaded 500 m grid resolution albedo data derived from MODIS imagery
are processed and mosaiced to return white-sky albedo only.
}
\seealso{
[albedo_download()], [albedo_fromaerial()], [albedo_adjust()]
}

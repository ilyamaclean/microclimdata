#' @title Download Digital Elevation Data
#'
#' @description Downloads Digital Elevation Data from the Amazon Web Server
#' @param r a SpatRaster object used as a template for downloading data with an EPSG code defined.
#' @param msk optional logical indicating whether to apply `r` as a mask to downloaded
#' data (see details),
#' @param zeroasna optional logical indicating whether to set values wiht an elevation of zero
#' (e.g. sea) as NA. Generally not needed if `msk = TRUE` (see details).
#' @return a SpatRaster of digital elevation data
#' @details the data returned will match the resolution, extent and coordinate reference
#' system of `r`. The resolution of available data varies globally however, so in some instances
#' it may be, that after resampling to match `r` the resulting data may be derived by
#' interpolation. A map showing the resolution of available data is available from
#' https://github.com/tilezen/joerd/blob/master/docs/images/footprints-preview.png.
#' For use with or use with microclimate packages, the returned SpatRaster, and hence
#' `r` should have a Mercator type Coordinate Reference System such that X, y and Z
#' units are all in metres. The data available from AWS blend digital elevation and
#' bathymetry data, the latter generally only available at coarser resolution. This
#' can cause oddities in coastal regions, where it is sensible to use `r` as a land-sea
#' mask by setting `msk` to `TRUE`.Alternatively by setting `zeroasna` to `TRUE` anything
#' below sea-level is set to `NA`.
#' @import elevatr
#' @export
dem_download<-function(r, msk = TRUE, zeroasna = FALSE) {
  if (!curl::has_internet()) {
    message("Please connect to the internet and try again.")
    return(NULL)
  }
  reso<-min(res(r))
  ll<-.latlongfromrast(r)
  z<-ceiling(log((cos(ll$lat*pi/180)*2*pi*6378137)/(256*reso),2))
  z[z>14]<-14
  prj<-sf::st_crs(r)
  dtm<-elevatr::get_aws_terrain(r,z,prj)
  dtm<-resample(dtm,r)
  if (msk) dtm<-mask(dtm, r)
  if (zeroasna) dtm[dtm<=0]<-NA
  return(dtm)
}
#' @title Downloads soil data
#'
#' @description Downloads data on soil properties from the https://soilgrids.org/
#' @param r a SpatRaster object used as a template for downloading data (see details).
#' @param pathdi directory to which data are temporarily or permenently downloaded (default working directory).
#' @param deletefiles optional logical indicating whether to delete the original downloaded files for each depth
#' @return a multi-layer SpatRaster of depth-weighted soil physical-properties:
#' #' \describe{
#'  \item{bdod}{Bulk density in cg/cm^3}
#'  \item{clay}{Mass fraction of clay particles g/kg}
#'  \item{sand}{Mass fraction of sand particles g/kg}
#'  \item{silt}{Mass fraction of silt particles g/kg}
#'  }
#' @seealso [soildata_downscale()]
#' @import httr
#' @import terra
#' @export
soildata_download<-function(r, pathdir = getwd(), deletefiles = TRUE) {
  # get bounding box
  dir.create(pathdir,showWarnings=FALSE)
  e<-ext(r)
  subsetx<-paste0("&SUBSET=X(",e$xmin,",",e$xmax,")")
  subsety<-paste0("&SUBSET=Y(",e$ymin,",",e$ymax,")")
  # get subsetting crs
  crs_info <- crs(r)
  epsg_code <- sub(".*ID\\[\"EPSG\",(\\d+)\\].*", "\\1", crs_info)
  subsetcrs<-paste0("&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/",epsg_code)
  outcrs<-paste0("&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/",epsg_code)
  # datasets
  # classes
  vars<-c("bdod","clay","sand","silt")
  dps_mn<-c(0,5,15,30,60,100)
  dps_mx<-c(5,15,30,60,100,200)
  ro<-r
  for (ii in 1:4) {
    for (jj in 1:6) {
      base_url<-paste0("https://maps.isric.org/mapserv?map=/map/",vars[ii],".map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage")
      idform<-paste0("&COVERAGEID=",vars[ii],"_",dps_mn[jj],"-",dps_mx[jj],"cm_Q0.5&FORMAT=GEOTIFF_INT16")
      request_url<-paste0(base_url,idform,subsetx,subsety,subsetcrs,outcrs)
      response <- GET(request_url)
      response
      if (response$status_code == 200) {
        fo<-paste0(pathdir,vars[ii],"_",dps_mn[jj],"_",dps_mx[jj],"cm.tif")
        writeBin(content(response, "raw"), fo)
      } else {
        stop(paste0("Bad query request for ",vars[ii]," ",dps_mn[jj],"-",dps_mx[jj],"cm"))
      }
    }
    dps<-c(5,10,15,30,40,100)
    wgts<-c(1,0.5,0.25,0.125,0.0625,0.03125)*dps
    wgts<-wgts/sum(wgts)
    fi<-paste0(pathdir,vars[ii],"_",dps_mn[1],"_",dps_mx[1],"cm.tif")
    ri<-rast(fi)*wgts[1]
    if(deletefiles) unlink(fi)
    for (jj in 2:6) {
      fi<-paste0(pathdir,vars[ii],"_",dps_mn[jj],"_",dps_mx[jj],"cm.tif")
      ri<-rast(fi)*wgts[jj]+ri
      if(deletefiles) unlink(fi)
    }
    if(ii == 1) ro<-resample(ro,ri)
    ro<-c(ro,ri)
  }
  ro<-ro[[-1]]
  names(ro)<-vars
  return(ro)
}
#' @title Downscales soil data using landcover data
#'
#' @description Fills NAs and Downscales soil data derived using `soildata_download` using
#' landcover data.
#' @param soildata a multi-layer SpatRaster of soil physical properties as returned by `soildata_download`.
#' @param landcover a high-resolution SpatRaster of landcover as returned by e.g. `lcover_download`
#' @param water value in land cover corresponding to water. Used to clean soil data
#' @return a multi-layer SpatRaster of downscaled depth-weighted soil physical-properties
#' @seealso [soildata_download()]
#' @import terra
#' @export
soildata_downscale <- function(soildata, landcover, water = 80) {
  # get coarse mask
  landcover[landcover == water] <- NA
  msk <- resample(landcover, soildata[[1]])
  bden <- soildata[[1]]
  bden[bden == 0] <- NA
  bden <- .fillna(bden, msk)
  # clay mass fraction
  clay <- soildata[[2]]
  clay[clay == 0] <- NA
  clay <- .fillna(clay, msk)
  # sand mass fraction
  sand <- soildata[[3]]
  sand[sand == 0] <- NA
  sand <- .fillna(sand, msk)
  # silt mass fraction
  silt <- soildata[[4]]
  silt[silt == 0] <- NA
  silt <- .fillna(silt, msk)
  soildata <- c(bden,clay,sand,silt)
  soildataf <- .sharpen(soildata[[1]], landcover)
  for (i in 2:4)  soildataf <- c(soildataf, .sharpen(soildata[[i]], landcover))
  return(soildataf)
}
#' @title Derives a soil type from soil data
#'
#' @description Derives a numerical code corresponding to likely soil type from bulk density and the
#' sand, clay and silt mass fraction.
#' @param soildata a multi-layer SpatRaster of soil physical properties as returned by `soildata_download`.
#' or `soildata_downscale`.
#' @return a SpatRaster of soil types one of:
#' \describe{
#'  \item{1}{Sand}
#'  \item{2}{Loamy sand}
#'  \item{3}{Sandy loam}
#'  \item{4}{Loam}
#'  \item{5}{Silt loam}
#'  \item{6}{Sandy clay loam}
#'  \item{7}{Clay loam}
#'  \item{8}{Silty clay loam}
#'  \item{9}{Sandy clay}
#'  \item{10}{Silty clay}
#'  \item{11}{Clay}
#' }
soildata_gettype <- function(soildata) {
  bden <- soildata[[1]] / 100
  clay <- soildata[[2]] / 1000
  sand <- soildata[[3]] / 1000
  silt <- soildata[[4]] / 1000
  # soil type
  soiltype <- .rast(getsoiltypecpp(as.matrix(bden, wide = TRUE), as.matrix(clay, wide = TRUE),
                                   as.matrix(sand, wide = TRUE), as.matrix(silt, wide = TRUE)), bden)
  return(soiltype)
}

#' @title Download MODIS albedo data
#'
#' @description Downloads 500m grid resolution MODIS albedo data from LPDAAC_ECS.
#' @param r a SpatRaster with a coordinate reference system defined. Used
#' to the location for which to download data.
#' @param tme POSIXlt object covering time period for which data are required
#' (see details).
#' @param pathout directory to which data will be downloaded.
#' @param credentials a data.frame of user credentials (see vignette for
#' setting this up)
#' @param pathtopython directory in which Python is saved (see
#' vignette instructions for installing and setting up Python)
#' @seealso [albedo_process()], [albedo_fromaerial()], [albedo_adjust()]
#' @details 500 m grid resolution albedo data derived from MODIS imagery are downloaded
#' from LPDAAC_ECS using the `luna` package. Data are available globally from 16th
#' February 2000 onward. Multiple files are downloaded and stored in the directory
#' given by `pathout`. Use [albedo_process()] to tile these together into one file.
#' See vignette for examples.
#' @import luna
#' @import terra
#' @export
#' @rdname albedo_download
albedo_download<-function(r, tme, pathout, credentials) {
  tsed<-tme[length(tme)]
  tmod<-as.POSIXlt(0,origin="2000-02-18",tz="UTC")
  if (tsed < tmod) stop("No data available prior to 2000-02-18")
  # Download MODIS DATA
  e<-ext(r)
  r2<-rast(e)
  crs(r2)<-crs(r)
  r2<-project(r2, "EPSG:4326")
  e2<-ext(r2)
  st<-substr(as.character(tme[1]),1,10)
  ed<-substr(as.character(tme[length(tme)]),1,10)
  mf<-luna::getNASA("MCD43A3", st, ed, aoi = e2, version = "061", download=FALSE)
  username <- credentials$username[1]
  password <- credentials$password[1]
  if (length(mf) > 0) {
    luna::getNASA("MCD43A3", st, ed, aoi = e2, version = "061", download=TRUE,
                  path=pathout,username=username,password=password,server="LPDAAC_ECS")
  } else {
    stop("No data for specified location or time period")
  }
}
#' @title Processes downloaded MODIS albedo data
#'
#' @description Processes 500m grid resolution MODIS albedo data downlaoded using
#' [albedo_download()] and returns a single terra::SpatRaster for the area
#' specified by `r`
#' @param r a SpatRaster with a coordinate reference system defined. Used
#' to the location for which to crop downloaded data too.
#' @param pathin directory to which data are downloaded. Same as `pathout` in
#' [albedo_download()]
#' @seealso [albedo_download()], [albedo_fromaerial()], [albedo_adjust()]
#' @details Downloaded 500 m grid resolution albedo data derived from MODIS imagery
#' are processed and mosaiced to return white-sky albedo only.
#' @returns a SpatRaster of white-sky albedo
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimdata, .registration = TRUE
#' @export
#' @rdname albedo_process
albedo_process<-function(r, pathin)  {
  lst <- list.files(pathin)
  if (length(lst) > 0) {
    pb <- utils::txtProgressBar(min = 0, max = length(lst) + 1, style = 3)
    fi <- paste0(pathin, lst[1])
    modr <- rast(fi)[[30]]
    # get extent of r
    bbx <- rast(ext(r))
    crs(bbx) <- crs(r)
    bbx <- project(bbx, crs(modr))
    e <- ext(bbx)
    e$xmin <- e$xmin - 1000
    e$xmax <- e$xmax + 1000
    e$ymin <- e$ymin - 1000
    e$ymax <- e$ymax + 1000
    modr <- extend(modr, e)
    modr <- crop(modr, e)
    utils::setTxtProgressBar(pb, 1)
    for (i in 2:length(lst)) {
      fi <- paste0(pathin, lst[i])
      mr <- rast(fi)[[30]]
      mr <- extend(mr, e)
      mr <- crop(mr, e)
      modr <- c(modr, mr)
      utils::setTxtProgressBar(pb, i)
    }
  } else stop("No files to process!")
  m <- apply3D(as.array(modr))
  albedo <- .rast(m, modr)
  albedo <- project(albedo, crs(r))
  albedo <-crop(albedo, ext(r))
  utils::setTxtProgressBar(pb, i + 1)
  return(albedo)
}
#' @title Derives albedo from aerial photographs
#'
#' @description Derives albedo from red-green-blue and colour-infrared aerial photographs
#' @param RGBimage a 3-layer red-green-blue SpatRaster with a coordinate reference system
#' defined and at least partially over-lapping with CIRimage.
#' @param CIRimage a 3-layer colour-infra-red SpatRaster with a coordinate reference system
#' defined and at least partially over-lapping with RGBimage.
#' @param RGBbandmins minimum wavelengths of bands in RGBimage (nm).
#' Default assumes image ordered Red-Green-Blue
#' @param RGBbandmaxs maximum wavelengths of bands in RGBimage (nm).
#' Default assumes image ordered Red-Green-Blue
#' @param CIRbandmins minimum wavelengths of bands in CIRimage (nm).
#' Default assumes image ordered Near-infrared, Red, Green
#' @param CIRbandmaxs maximum wavelengths of bands in CIRimage (nm).
#' Default assumes image ordered Near-infrared, Red, Green
#' @seealso [albedo_adjust()]
#' @details derives albedo by weighting each image band by the band width and spectral
#' density of light emitted from the sun in that band range. Note that technically this
#' function does not derive true albedo as the brightness and contrast of the aerial
#' image may have been artificially altered during processing. use [albedo_adjust()] to
#' deal with this  problem
#' @returns a SpatRaster of albedo
#' @import terra
#' @export
#' @rdname albedo_fromaerial
albedo_fromaerial <- function(RGBimage, CIRimage, RGBbandmins = c(620, 495, 450),
                              RGBbandmaxs = c(750, 570, 495),
                              CIRbandmins = c(750, 620, 495),
                              CIRbandmaxs = c(900, 750, 570)) {
  # Check resolutions and resample to coarser resolutions if not matching
  res1 <- res(RGBimage)[1]
  res2 <- res(CIRimage)[1]
  if (res1 > res2) CIRimage <- resample(CIRimage, RGBimage)
  if (res2 > res1) RGBimage <- resample(RGBimage, CIRimage)
  # Check extents and intersect
  e1 <- ext(RGBimage)
  e2 <- ext(CIRimage)
  e <- ext(max(e1$xmin, e2$xmin), min(e1$xmax, e2$xmax),
           max(e1$ymin, e2$ymin), min(e1$ymax, e2$ymax))
  RGBimage <- crop(RGBimage, e)
  CIRimage <- crop(CIRimage, e)
  # Create weights
  wgt1 <- sum(.Planck(seq(RGBbandmins[1], RGBbandmaxs[1], by = 1))/(10^13))
  wgt2 <- sum(.Planck(seq(RGBbandmins[2], RGBbandmaxs[2], by = 1))/(10^13))
  wgt3 <- sum(.Planck(seq(RGBbandmins[3], RGBbandmaxs[3], by = 1))/(10^13))
  wgt4 <- sum(.Planck(seq(CIRbandmins[1], CIRbandmaxs[1], by = 1))/(10^13))
  wgt5 <- sum(.Planck(seq(CIRbandmins[2], CIRbandmaxs[2], by = 1))/(10^13))
  wgt6 <- sum(.Planck(seq(CIRbandmins[3], CIRbandmaxs[3], by = 1))/(10^13))
  # calculate albedo
  albedo <- (RGBimage[[1]] * wgt1 + RGBimage[[2]] * wgt2 + RGBimage[[3]] * wgt3 +
               CIRimage[[1]] * wgt4 + CIRimage[[2]] * wgt5 + CIRimage[[3]] * wgt6)  /
    (wgt1 + wgt2 + wgt3 + wgt4 + wgt5 + wgt6)
  albedo <- albedo / 255
  albedo <- mask(albedo, RGBimage[[1]])
  return(albedo)
}
#' @title Corrects albedo derived from aerial photographs
#'
#' @description Using MODIS, derived albedo, corrects albedo derived using function
#' [albedo_fromaerial()] to correct for brightness and contrast in the aerial image
#' @param photoalbedo a SpatRaster of albedo derived using [albedo_fromaerial()].
#' @param modisalbedo a SpatRaster of albedo derived from MODIS data using [albedo_download()]
#' and [albedo_process()] usually of much coarser resoltuion than `photoalbedo`.
#' @returns a SpatRaster of adjusted albedo
#' @import terra
#' @export
#' @rdname albedo_adjust
albedo_adjust<-function(photoalbedo, modisalbedo) {
  # crop modis to match extent of aerial
  if (crs(photoalbedo) != crs(modisalbedo)) modisalbedo<-project(modisalbedo, crs(photoalbedo))
  v1 <- as.vector(crop(modisalbedo, ext(photoalbedo)))
  v1 <-v1[is.na(v1) == FALSE]
  v2 <- as.vector(resample(photoalbedo, modisalbedo))
  v2 <-v2[is.na(v2) == FALSE]
  # logit transform
  lv1 <- log(v1/ (1 - v1))
  lv2 <- log(v2/ (1 - v2))
  # if length of either is 1
  n <- min(length(v1), length(v2))
  if (n == 1) {
    mu <- lv1 - lv2
    lphoto <- log(photoalbedo / (1 - photoalbedo)) + mu
    albedo <- 1 / (1 + exp(-lphoto))
  } else {
    mum <- mean(lv1) - mean(lv2)
    mus <- sd(lv1) / sd(lv2)
    lphoto <- log(photoalbedo / (1 - photoalbedo))
    me <- mean(as.vector(lphoto), na.rm = TRUE)
    lphoto <- ((lphoto - me) * mus) + mum + me
    albedo <- 1 / (1 + exp(-lphoto))
  }
  return(albedo)
}
#' @title Derives leaf and ground reflectance
#'
#' @description Derives leaf and ground reflectance from albedo.
#'
#' @param alb a SpatRaster of albedo values.
#' @param lai a SpatRaster of leaf afrea index values.
#' @param x a SpatRaster of leaf inclination coefficient values as returned by [x_calc()].
#' @param plotprogress optional logical indicating whether or not to keep track ofprogress by plotting interim results.
#' @param maxiter maximum number of iterations
#' @param tol convergence tolerance
#' @param bwgt value between zero and one indicating by how much to back-weight results during iteration. HIgher values
#' ensure more stable convergence, but may require more iteration to achieve convergence.
#' @returns a list comprising SpatRasters of leaf and ground reflectance
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimdata, .registration = TRUE
#' @import terra
#' @export
#' @rdname reflectance_calc
reflectance_calc <- function(alb, lai, x, plotprogress = TRUE, maxiter = 50, tol = 0.001, bwgt = 0.5) {
  # calculate intersecting extents
  e1 <- intersect(ext(lai), ext(alb))
  e <- intersect(e1, ext(x))
  lai <- crop(lai, e)
  alb <- crop(alb, e)
  x <- crop(x, e)
  # check dims
  all_same <- compareGeom(lai, alb, x)
  if (all_same) {
    tst <- exp(-mean(as.vector(lai), na.rm = T))
    lref <- (x * 0 + 0.5) * (1 - wgt) + wgt * alb
    gref <- x * 0 + 0.15
    mxdif <- tol * 10
    paim <- as.matrix(lai, wide = TRUE)
    xm <- as.matrix(x, wide = TRUE)
    albm <- as.matrix(alb, wide = TRUE)
    itr <- 1
    while (mxdif > tol) {
      if (tst < 0.5) {
        lref2 <- .rast(find_lref(paim, as.matrix(gref, wide = TRUE), xm, albm), x)
        lref2 <- .fillna(lref2, x, zerotoNA = FALSE)
        gref2 <- .rast(find_gref(as.matrix(lref2, wide = TRUE), paim, xm, albm), x)
        gref2 <- .fillna(gref2, x, zerotoNA = FALSE)
      } else {
        gref2 <- .rast(find_gref(as.matrix(lref, wide = TRUE), paim, xm, albm), x)
        gref2 <- .fillna(gref2, x, zerotoNA = FALSE)
        lref2 <- .rast(find_lref(paim, as.matrix(gref, wide = TRUE), xm, albm), x)
        lref2 <- .fillna(lref2, x, zerotoNA=FALSE)
      }
      gref <- bwgt * gref + (1 - bwgt) * gref2
      lref <- bwgt * lref + (1 - bwgt) * lref2
      mxdif1 <- mean(abs(as.vector(gref) - as.vector(gref2)), na.rm = TRUE)
      mxdif2 <- mean(abs(as.vector(lref) - as.vector(lref2)), na.rm = TRUE)
      mxdif <- max(mxdif1, mxdif2)
      if (plotprogress & itr%%3 == 0) {
        tp1 <- paste0("Ground difference from previous: ", round(mxdif1, 4))
        tp2 <- paste0("Leaf difference from previous: ", round(mxdif2, 4))
        par(mfrow = c(1, 2))
        plot(gref, main = tp1, cex.main = 1)
        plot(lref, main = tp2, cex.main = 1)
      }
      itr <- itr+1
      if (itr > maxiter) mxdif <- 0
    }
  } else (stop("Geometries of input rasters do not match"))
  return(list(gref = gref, lref = lref))
}
#' @title Create vegetation inputs for grid model
#'
#' @description Creates an object of class vegparams used as an input to `microclimf`
#' or `microgrid`.
#'
#' @param soildata a multilayer SpatRaster of soil data as returned by [soildata_download()].
#' @param refldata a list of ground and leaf reflectances as returned by [reflectance_calc()].
#' @param landcover a SpatRaster of landcover class as returned by [lcover_download()].
#' @param water integer value in `landcover` corresponding to water (e.g. 80 for ESA dataset)
#' @returns an object of class `soilcharac`, namely a list of the following objects:
#' \describe{
#'   \item{soiltype}{a PackedSpatRaster of soil types specified using integer values}
#'   \item{groundr}{a PackedSpatRaster object of ground reflectances}
#'   \item{rho}{a PackedSpatRaster object of bulk densities (Mg / m^3)}
#'   \item{Vm}{a PackedSpatRaster object of the volumetric mineral fraction of soil}
#'   \item{Vq}{a PackedSpatRaster object of the volumetric quartz fraction of soil}
#'   \item{Mc}{a PackedSpatRaster object of the mass fraction of clay}
#' }
#' @import terra
#' @import httr
#' @export
#' @rdname create_soilgrid
create_soilgrid <- function(soildata, refldata, landcover, water = 80) {
  # downscale soil data using landcover
  soildata <- soildata_downscale(soildata, landcover, water)
  # get soil type
  soiltype <- soildata_gettype(soildata)
  names(soiltype) <- "soiltype"
  names(refldata$gref) <- "gorund reflectance"
  organic <- 1 - (soildata$clay + soildata$sand + soildata$silt) /1000
  soilp <- list(soiltype = wrap(soiltype),
                groundr = wrap(refldata$gref),
                rho = wrap(soildata$bdod / 100),
                Vm =  wrap((soildata$bdod /265) * (1- organic)),
                Vq = wrap((soildata$bdod / 265) * ((soildata$sand / 1000) * 0.9
                            + (soildata$silt / 1000) * 0.5)),
                Mc = wrap(soildata$clay / 1000))
  class(soilp) <- "soilcharac"
  return(soilp)
}
#' @title Create Soil inputs for grid model
#'
#' @description Creates an object of class soilparams used as an input to `micropoint`
#'
#' @param soildata a multilayer SpatRaster of soil data as returned by [soildata_download()].
#' @param refldata a list of ground and leaf reflectances as returned by [reflectance_calc()].
#' @param dtm a SpatRaster of elevations as returned by [dem_download()].
#' @param landcover a SpatRaster of landcover class as returned by [lcover_download()].
#' @param water integer value in `landcover` corresponding to water (e.g. 80 for ESA dataset)
#' @param lat optionally, latitude (decimal degrees) - see details.
#' @param long optionally, longitude (decimal degrees) - see details.
#' @details If `lat` and `long` are unspecified, point data are selelected for
#' the location and the centre of supplied `landcover`.
#' @returns an object of class `groundparams`, namely a list of the following objects:
#' \describe{
#'   \item{gref}{ground reflectance}
#'   \item{slope}{slope of ground surface (deg)}
#'   \item{aspect}{aspect of ground surface (deg from north)}
#'   \item{em}{emissivity of ground surface}
#'   \item{rho}{Soil bulk density (Mg / m^3)}
#'   \item{Vm}{Volumetric mineral fraction of soil}
#'   \item{Vq}{Volumetric quartz fraction of soil}
#'   \item{Mc}{Mass fraction of clay}
#'   \item{b}{Shape parameter for Campbell soil moisture model}
#'   \item{Psie}{Matric potential (J / m^3)}
#'   \item{Smax}{Volumetric water content at saturation}
#'   \item{Smin}{Residual water content}
#'   \item{alpha}{Shape parameter of the van Genuchten model (cm^-1)}
#'   \item{n}{Pore size distribution parameter}
#'   \item{Ksat}{Saturated hydraulic conductivity (cm / day)}
#' }
#' @import terra
#' @import httr
#' @export
#' @rdname create_soilpoint
create_soilpoint <- function(soildata, refldata, dtm, landcover, lat = NA, long = NA, water = 80) {
  if (class(lat) == "logical") {
    e<-ext(landcover)
    xy <- data.frame(x = (e$xmin + e$xmax) / 2,
                     y = (e$ymin + e$ymax) / 2)
  } else {
    # convert latlong to local crs
    xy <- data.frame(x = long, y = lat)
    xy <- sf::st_as_sf(xy, coords = c("x", "y"),
                       crs = 4326)
    xy <- sf::st_transform(xy, crs(landcover))
    xy <- data.frame(x = sf::st_coordinates(xy)[1], y = sf::st_coordinates(xy)[2])
  }
  soilcg <- create_soilgrid(soildata, refldata, landcover, water)
  soiltype <- rast(soilcg$soiltype)
  groundr <- rast(soilcg$groundr)
  rho <- rast(soilcg$rho)
  Vm <- rast(soilcg$Vm)
  Vq <- rast(soilcg$Vq)
  Mc <- rast(soilcg$Mc)
  # calaculate slope
  slope <- terrain(dtm, v = "slope", unit = "degrees")
  aspect <- terrain(dtm, v = "aspect", unit = "degrees")

  # Extract 1
  soilt <- extract(soiltype, xy)[,2]
  groundp <- list(gref = extract(groundr, xy)[,2],
                  slope = extract(slope, xy)[,2],
                  aspect = extract(aspect, xy)[,2],
                  em = 0.97,
                  rho = extract(rho, xy)[,2],
                  Vm = extract(Vm, xy)[,2],
                  Vq = extract(Vq, xy)[,2],
                  Mc = extract(Mc, xy)[,2])
  # get additional variables
  groundp$b <- soiltable$b[soilt]
  groundp$Psie <- -soiltable$psi_e[soilt]
  groundp$Smax <- soiltable$Smax[soilt]
  groundp$Smin <- soiltable$Smin[soilt]
  groundp$alpha <- soiltable$alpha[soilt]
  groundp$n <- soiltable$n[soilt]
  groundp$Ksat <- soiltable$Ksat[soilt]
  class(groundp) <- "groundparams"
  return(groundp)
}


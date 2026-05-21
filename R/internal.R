#' Check if input is a SpatRaster or PackedSpatRaster and convert to matrix if it is
#' @import terra
.is <- function(r) {
  if (class(r)[1] == "PackedSpatRaster") r<-rast(r)
  if (class(r)[1] != "matrix") {
    if (dim(r)[3] > 1) {
      y<-as.array(r)
    } else y<-as.matrix(r,wide=TRUE)
  } else y<-r
  y
}
#' Create SpatRaster object using a template
#' @import terra
.rast <- function(m,tem) {
  r<-rast(m)
  ext(r)<-ext(tem)
  crs(r)<-crs(tem)
  r
}
#' Create bounding box from raster
.bbox<-function(r) {
  e<-ext(r)
  xy<-data.frame(x = c(e$xmin,e$xmin,e$xmax,e$xmax),
                 y = c(e$ymin,e$ymax,e$ymin,e$ymax))
  xy <- sf::st_as_sf(xy, coords = c("x", "y"),
                     crs = crs(r))
  ll <- sf::st_transform(xy, 4326)
  ll <- data.frame(sf::st_coordinates(ll))
  out <- c(min(ll$X),min(ll$Y),max(ll$X),max(ll$Y))
  return(out)
}
#' Checks whether file has been downloaded completely
#' @import terra
.chcktif <- function(file_path) {
  result <- tryCatch({
    # Attempt to read the file
    r <- rast(file_path)
    # Perform a basic check to ensure the raster is valid
    sample_values <- values(r, cells = seq_len(min(100, ncell(r))))
    # If successful, return 1
    return(1)
  }, error = function(e) {
    # Handle the error and log a message
    message(paste("Error processing file:", file_path, ":", e$message))
    # Return NA to indicate the file was not processed successfully
    return(NA)
  })
  return(result)
}
#' Applies quality mask to lai data
#' @import terra
.cleanlai<-function(lai,qual) {
  e1<-ext(lai)
  e2<-ext(qual)
  e<-e1
  e$xmin<-min(e1$xmin,e2$xmin)
  e$ymin<-min(e1$ymin,e2$ymin)
  e$xmax<-max(e1$xmax,e2$xmax)
  e$ymax<-max(e1$ymax,e2$ymax)
  lai<-extend(lai,e)
  qual<-extend(qual,e)
  qual[qual>1]<-NA
  lai<-mask(lai,qual)
  return(lai)
}
#' Finds majority class in aggregated SpaRaster grid cell
.majority_class <- function(x) {
  # Remove NAs
  x <- na.omit(x)
  if (length(x) == 0) return(NA)
  freq_table <- table(x)
  majority <- names(freq_table)[which.max(freq_table)]
  # Check if the majority class makes up more than half of the values
  if (max(freq_table) > (length(x)*0.75)) {
    return(as.numeric(majority))
  } else {
    return(NA)
  }
}
#' Sorts out NAs in a raster dataset due to coarsening
.rastna<-function(r,af,msk) {
  r<-mask(r,msk)
  rc<-resample(aggregate(r,af,"mean",na.rm=TRUE),r)
  m<-.is(r)
  mc<-.is(rc)
  s<-which(is.na(m))
  m[s]<-mc[s]
  ro<-.rast(m,r)
  return(ro)
}
# fills NA values in a SpatRaster
.fillna<-function(ri,msk,zerotoNA=TRUE) {
  if (zerotoNA) ri[ri==0]<-NA
  mx<-min(dim(ri)[1:2])
  if (mx > 2) {
    r2<-.rastna(ri,2,msk)
    ri[is.na(ri)]<-r2[is.na(ri)]
  }
  if (mx > 4) {
    r2<-.rastna(ri,4,msk)
    ri[is.na(ri)]<-r2[is.na(ri)]
  }
  if (mx > 8) {
    r2<-.rastna(ri,8,msk)
    ri[is.na(ri)]<-r2[is.na(ri)]
  }
  if (mx > 16) {
    r2<-.rastna(ri,16,msk)
    ri[is.na(ri)]<-r2[is.na(ri)]
  }
  mm<-mean(as.vector(ri),na.rm=TRUE)
  ri[is.na(ri)]<-mm
  ri<-mask(ri,msk)
  return(ri)
}
#' Used to monitor downloads using rgee
#' @import rgee
.monitor_task <- function(task_id, sleep_time = 30) {
  repeat {
    # Get task status
    task_info <- ee_monitoring(task_id)

    # Print task status
    cat("Task status: ", task_info$state, "\n")

    # Check if the task is completed
    if (task_info$state %in% c('COMPLETED', 'FAILED', 'CANCELLED')) {
      cat("Task finished with status: ", task_info$state, "\n")
      break
    }
    # Sleep for the specified time before checking again
    Sys.sleep(sleep_time)
  }
}
#' performs circular NA approx
.circular_approx <- function(x) {
  # Create a wrapped version of the vector by appending the first element at the end
  wrapped_x <- c(x, x[1])
  # Find the positions of the non-NA values (ignoring the last element added for wrapping)
  non_na_pos <- which(!is.na(wrapped_x[-length(wrapped_x)]))

  # Perform linear interpolation using approx, specifying that x should be treated as circular
  approx_result <- approx(x = c(non_na_pos, length(x) + 1), # Add the last element for circular wrap
                          y = wrapped_x[c(non_na_pos, 1)],   # Wrap the last element with the first
                          xout = seq_along(x),
                          rule = 2)$y
  return(approx_result)
}
#' Get latitude and longitude of centre of r
.latlongfromrast<-function (r) {
  e <- ext(r)
  xy <- data.frame(x = (e$xmin + e$xmax)/2, y = (e$ymin + e$ymax)/2)
  xy <- sf::st_as_sf(xy, coords = c("x", "y"),
                     crs = crs(r))
  ll <- sf::st_transform(xy, 4326)
  ll <- data.frame(lat = sf::st_coordinates(ll)[2], long = sf::st_coordinates(ll)[1])
  return(ll)
}
#' Latitudes from SpatRaster object
.latsfromr <- function(r) {
  e <- ext(r)
  lts <- rep(seq(e$ymax - res(r)[2] / 2, e$ymin + res(r)[2] / 2, length.out = dim(r)[1]), dim(r)[2])
  lts <- array(lts, dim = dim(r)[1:2])
  lts
}
#' Longitudes from SpatRaster object
.lonsfromr <- function(r) {
  e <- ext(r)
  lns <- rep(seq(e$xmin + res(r)[1] / 2, e$xmax - res(r)[1] / 2, length.out = dim(r)[2]), dim(r)[1])
  lns <- lns[order(lns)]
  lns <- array(lns, dim = dim(r)[1:2])
  lns
}
#' lats and longs from SpatRaster object, including reprojection
#' @import terra
.latslonsfromr <- function(r) {
  lats<-.latsfromr(r)
  lons<-.lonsfromr(r)
  xy<-data.frame(x=as.vector(lons),y=as.vector(lats))
  xy <- sf::st_as_sf(xy, coords = c('x', 'y'), crs = crs(r))
  ll <- sf::st_transform(xy, 4326)
  ll <- data.frame(lat = sf::st_coordinates(ll)[,2],
                   long = sf::st_coordinates(ll)[,1])
  lons<-array(ll$long,dim=dim(lons))
  lats<-array(ll$lat,dim=dim(lats))
  return(list(lats=lats,lons=lons))
}
#' Download on climate variable using the RNCEP package
#' @import RNCEP
#' @import terra
.get_ncepone<-function(r, tme, climvar) {
  # buffer time sequence by a day
  st<-as.numeric(tme[1])-24*3600
  ed<-as.numeric(tme[length(tme)])+24*3600
  tme2<-as.POSIXlt(seq(st,ed,by=3600*24),origin="1970-01-01 00:00", tz="UTC")
  # get lats and longs range
  e<-ext(r)
  xy<-data.frame(x=c(e$xmin,e$xmin,e$xmax,e$xmax),
                 y=c(e$ymin,e$ymax,e$ymin,e$ymax))
  xy <- sf::st_as_sf(xy, coords = c("x", "y"),
                     crs = crs(r))
  ll <- sf::st_transform(xy, 4326)
  ll <- data.frame(lat = sf::st_coordinates(ll)[,2], long = sf::st_coordinates(ll)[,1])
  # Get years
  yrs<-unique(tme2$year)+1900
  ro<-list()
  # Download data for each year and store as list
  tta<-as.POSIXlt(0,origin="1900-01-01 00", tz = "UTC")
  for (yr in 1:length(yrs)) {
    s<-which(tme2$year+1900==yrs[yr])
    mths<-unique(tme2$mon[s]+1)
    a<-NCEP.gather(climvar, level = 'gaussian',
                   years.minmax = c(yrs[yr],yrs[yr]),
                   months.minmax = c(min(mths):max(mths)),
                   lat.southnorth = c(min(ll$lat),max(ll$lat)),
                   lon.westeast = c(min(ll$lon),max(ll$lon)),
                   reanalysis2 = FALSE, return.units = FALSE, status.bar = FALSE)
    tt<-unlist(dimnames(a)[3])
    tt<-strptime(tt, format = "%Y_%m_%d_%H", tz = "UTC")
    tta<-c(tta,tt)
    x<-as.numeric(colnames(a))
    y<-as.numeric(row.names(a))
    eo<-ext(min(x)-1.875/2,max(x)+1.875/2,
            min(y)-1.9047/2,max(y)+1.9047/2)
    ri<-rast(a)
    ext(ri)<-eo
    crs(ri)<-"EPSG:4326"
    ro[[yr]]<-ri
  }
  tta<-tta[-1]
  # Convert to SpatRast
  roo<-ro[[1]]
  if (length(yrs) > 1) {
    for (yr in 2:length(yrs)) roo<-c(roo,ro[[yr]])
  }
  return(list(roo=roo,tme=tta))
}
#' Extracts a NCEP variable from a SpatRaster returtened by get_ncep
#' @import terra
.extractvar<-function(cr,ll) {
  cr<-rast(cr)
  e<-ext(cr)
  if (e$xmax > 180) cr<-terra::shift(cr,dx=-360)
  xy<-data.frame(x=ll$lon,y=ll$lat)
  v<-as.numeric(extract(cr,xy)[1,])
  v<-v[-1]
  return(v)
}
#' Calculates diurnal temperature range adjustment for NCEP coastal grid cells
#' @import terra
.calcmu<-function(lat,lon) {
  if (lon > 180) lon<-lon-360
  lsea<-rast(landsea)
  lsf<-rast(landseafrac)
  xy<-data.frame(x=lon,y=lat)
  bin<-as.numeric(extract(lsea,xy)[1,])
  frac<-as.numeric(extract(lsf,xy)[1,])
  bin<-bin[-1]
  frac<-frac[-1]
  mu<-1
  if (bin == 0 && frac !=0) {
    lfrac<-log(frac/(1-frac))
    ldtr<-(0.770381-0.056130 *lfrac)
    mu<-exp(ldtr)/1.841305
    mu[mu<1]<-1
  }
  return(mu)
}
#' Converts NCEP tmax and tmin data to hourly
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimdata, .registration = TRUE
.hourlytemp<-function(tmn,tmx,tmeb,lat,lon,mu=1) {
  tmet<-as.POSIXlt(tmeb-3600*6)
  st<-which(tmet$hour==0)[1]
  ed<-which(tmet$hour==18)
  ed<-ed[length(ed)]
  s<-c(st:ed)
  tmnd<-matrix(tmn[s],ncol=4,byrow=TRUE)
  tmxd<-matrix(tmx[s],ncol=4,byrow=TRUE)
  tmnd<-apply(tmnd,1,min)
  tmxd<-apply(tmxd,1,max)
  dtr<-(tmxd-tmnd)*mu
  te<-(tmxd+tmnd)/2
  tmxd<-te+0.5*dtr
  tmnd<-te-0.5*dtr
  tmed <- matrix(as.numeric(tmet[s]), ncol = 4, byrow = TRUE)
  tmed <- apply(tmed, 1, mean)
  tmed <- as.POSIXlt(tmed, origin = "1970-01-01 00:00", tz = "UTC")
  year<-tmed$year+1900
  mth<-tmed$mon+1
  day<-tmed$mday
  if (lon > 180) lon<-lon-360
  tc<-hourlytempv(tmnd-273.15,tmxd-273.15,year,mth,day,lat,lon)
  tmed<-tmet[s]
  st<-tmed[1]
  ed<-tmed[length(tmed)]
  tmeh<-as.POSIXlt(seq(st,ed,by=3600),tz="UTC")
  return(list(tc=tc,tmeh=tmeh))
}
#' Converts all NCEP data to hourly for a single location
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimdata, .registration = TRUE
.hourlyall<-function(tmn,tmx,pk6,sp6,sw6,lw6,u6,v6,pr6,tmeb,tme,lat,lon) {
  mu<-.calcmu(lat,lon)
  tcc<-.hourlytemp(tmn,tmx,tmeb,lat,lon,mu)
  mn<-min(tme)
  mx<-max(tme)
  s<-which(tcc$tmeh>=mn & tcc$tmeh <=mx)
  tc<-tcc$tc[s]
  # (2) Pressure
  n<-(length(tmeb)-1)*6+1
  tme1<-tmeb+3600*6
  pkht<-spline(pk6~tme1,n=n)
  s<-which(pkht$x>=mn & pkht$x<=mx)
  pk<-pkht$y[s]/1000
  # (3) Relative humidity
  sph<-spline(sp6~tme1,n=n)$y[s]
  rh<-converthumidityCpp(sph, tc, pk)
  # (4) Shortwave radiation
  csr<-clearskyradv(tmeb$year+1900,tmeb$mon+1,tmeb$mday,tmeb$hour,lat,lon%%360,0,rep(1,length(tmeb$year)))
  csf<-sw6/csr
  csf[csf>1]<-1
  csf[is.infinite(csf)]<-NA
  if (is.na(csf[1])) csf[1]<-0.5
  if (is.na(csf[length(csf)])) csf[length(csf)]<-0.5
  csf<-zoo::na.approx(csf)
  tme2<-tmeb+3600*3
  csfht<-spline(csf~tme2,n=n)
  s2<-which(csfht$x>=mn & csfht$x<=mx)
  csfh<-csfht$y[s2]
  csfh[csfh<0]<-0
  csfh[csfh>1]<-1
  csrh<-clearskyradhourly(tme$year+1900,tme$mon+1,tme$mday,tme$hour,lat,lon%%360,0,rep(1,length(tme)))
  swr<-csfh*csrh
  # (5) Diffuse radiation
  swd<-difpropvCpp(swr,tme$year+1900,tme$mon+1,tme$mday,tme$hour,lat,lon%%360)*swr
  # (6) Longwave radiation
  tc6<-tempsix(tcc$tc)
  skyem<-lw6[1:length(tc6)]/(5.67*10^-8*(tc6+273.15)^4)
  skyem[skyem>1]<-1
  skyem[skyem<0]<-0
  n2<-(length(tc6)-1)*6+1
  skt<-spline(skyem~tme2[1:length(tc6)],n=n2)
  s3<-which(skt$x>=mn & skt$x<=mx)
  lwr<-skt$y[s3]*5.67*10^-8*(tc+273.15)^4
  # (7) Wind speed and direction
  u<-spline(u6~tme1,n=n)$y
  v<-spline(v6~tme1,n=n)$y
  ws<-sqrt(u[s]^2+v[s]^2)*0.7477849
  wd<-(atan2(u[s],v[s])*180/pi+180)%%360
  # (8) Precipitation
  ph<-prectohour(pr6)[s]
  return(list(temp=tc,relhum=rh,pres=pk,swdown=swr,difrad=swd,lwdown=lwr,windspeed=ws,winddir=wd,precip=ph))
}

#' fill in NAs and expand daily sst data to hourly
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimdata, .registration = TRUE
.expand_sst<-function(sst) {
  r<-sst[[1]]*0+1
  r[is.na(r)]<-1
  sst[is.nan(sst)] <- NA
  sst<-.rast(na_fill_array(.is(sst),.is(r),2),r)
  sst<-.rast(expandssttohour(.is(sst)),r)
  return(sst)
}
# ' equivelent of mcera5::extract_clima but corrected
.extract_clima <- function(nc, r, resampleout = FALSE) {
  # get time
  nc_file <- ncdf4::nc_open(nc)
  tme <- as.POSIXlt(ncdf4::ncvar_get(nc_file, "valid_time"),
                    origin="1970-01-01 00:00", tz = "UTC")
  ncdf4::nc_close(nc_file)
  # extract variables
  varn <- c("t2m", "d2m", "sp", "u10" , "v10",  "tp", "msdwlwrf", "fdir", "ssrd", "lsm")
  rlst <- list()
  for (i in 1:9) rlst[[i]]<-rast(nc, subds=varn[i])
  rlst[[10]]<-rast(nc, subds=varn[10])[[1]]
  lsm <-  rlst[[10]]
  tc <- rlst[[1]]-273.15
  # coastal correction
  if (any(terra::values(lsm) < 1)) {
    # Calculate daily average
    # Indices to associate each layer with its yday
    ind <- rep(1:(dim(tc)[3]/24), each = 24)
    # Average across days
    tmean <- terra::tapp(tc, ind, fun = mean, na.rm = T)
    # Repeat the stack 24 times to expand back out to original timeseries
    tmean <- rep(tmean, 24)
    # Sort according to names so that the stack is now in correct order: each
    # daily mean, repeated 24 times
    # your pasted command properly sorts the names, X1 to X365 (or X366)
    tmean <- tmean[[paste0("X", sort(rep(seq(1:(dim(tc)[3]/24)), 24)))]]
    m <- (1 - lsm) * 1.285 + 1
    tdif <- (tc - tmean) * m
    tc<- tmean + tdif
  }
  rlst[[1]]<-tc
  # crop variables
  e<-ext(r)
  rr<-rast(e)
  crs(rr)<-crs(r)
  rll<-project(rr,"EPSG:4326")
  ell<-ext(rll)
  for (i in 1:10) rlst[[i]]<-crop(rlst[[i]],ell,snap='out')
  if (resampleout) {
    for (i in 1:10) {
      if (crs(r) != crs(rll)) {
        rlst[[i]]<-project(rlst[[i]],r)
      } else {
        rlst[[i]]<-resample(rlst[[i]],r)
      }
      rlst[[i]]<-mask(rlst[[i]],r)
    }
  }
  rlsto<-list()
  rte<-rlst[[1]]
  rte<-rte[[1]]
  # convert variables
  rlsto[[1]]<-rlst[[1]]
  rh <- .rast((satvapCpp(.is(rlst[[2]])-273.15) /  satvapCpp(.is(rlst[[1]]))) * 100, rte)
  rh[rh > 100]<-100
  rlsto[[2]]<-rh
  rlsto[[3]] <- rlst[[3]]/1000
  rlsto[[4]] <- rlst[[9]]/3600
  dni <- rlst[[8]]/3600
  ll<-.latslonsfromr(rte)
  si <- .rast(solarindexarray(tme$year+1900, tme$mon+1, tme$mday, tme$hour, ll$lats, ll$lons), rte)
  rlsto[[5]]  <- rlsto[[4]]  - (si * dni)
  rlsto[[6]]  <- rlst[[7]]
  rlsto[[7]]<-sqrt(rlst[[4]]^2+rlst[[5]]^2)*0.7477849 # Wind speed (m/s)
  rlsto[[8]]<-(atan2(rlst[[4]],rlst[[5]])*180/pi+180)%%360
  rlsto[[9]]<- rlst[[6]] * 1000
  names(rlsto)<-c("temp","relhum","pres","swdown","difrad","lwdown","windspeed","winddir","precip")
  for (i in 1:9)  time(rlsto[[i]])<-as.POSIXct(tme)
  return(rlsto)
}
#' get the right UKCP decade
.find_ukcp_decade<-function(collection=c('land-gcm','land-rcm'),startdate,enddate){
  collection<-match.arg(collection)
  if(class(startdate)[1]!="POSIXlt" | class(enddate)[1]!="POSIXlt") stop("Date parameters NOT POSIXlt class!!")
  rcm_decades<-c('19801201-19901130','19901201-20001130','20001201-20101130','20101201-20201130',
                 '20201201-20301130','20301201-20401130','20401201-20501130','20501201-20601130',
                 '20601201-20701130','20701201-20801130')
  gcm_decades<-c('19791201-19891130','19891201-19991130','19991201-20091130','20091201-20191130',
                 '20191201-20291130','20291201-20391130','20391201-20491130','20491201-20591130',
                 '20591201-20691130','20691201-20791130')
  if(collection=='land-gcm') ukcp_decades<-gcm_decades else ukcp_decades<-rcm_decades
  decade_start<-lubridate::ymd(sapply(strsplit(ukcp_decades,"-"), `[`, 1))
  decade_end<-lubridate::ymd(sapply(strsplit(ukcp_decades,"-"), `[`, 2))
  decades<-ukcp_decades[which(decade_end>startdate & decade_start<enddate)]
  return(decades)
}
#' extract dates from downloaded ukcp18 nc file
.get_ukcp18_dates<-function(ncfile){
  netcdf_data <-ncdf4::nc_open(ncfile)
  time_hours <- ncdf4::ncvar_get(netcdf_data,"time")
  ncdf4::nc_close(netcdf_data)
  years<-floor(time_hours/(360*24))+1970
  months<-ceiling((time_hours-(years-1970)*(360*24)) / (30*24) )
  days<-ceiling( (time_hours-((years-1970)*(360*24)) - ((months-1)*30*24) )/ 24 )
  hours<-time_hours - ((years-1970)*(360*24)) - ((months-1)*30*24) - ((days-1)*24)
  ukcp_dates<-paste(years,sprintf("%02d",months),sprintf("%02d",days),sep="-")
  return(ukcp_dates)
}
#' correct dates form downloaded ukcp18 nc file
.correct_ukcp_dates<-function(ukcp_dates){
  years<-as.numeric(sapply(strsplit(ukcp_dates,"-"), getElement, 1))
  months<-as.numeric(sapply(strsplit(ukcp_dates,"-"), getElement, 2))
  days<-as.numeric(sapply(strsplit(ukcp_dates,"-"), getElement, 3))

  # Shift March day by plus one
  sel<-which(months==3)
  days[sel] <- days[sel] + 1

  # Shift Feb day by minus 1
  sel <- which(months==2)
  days[sel] <-  days[sel] - 1

  # 1st of Feb as 31 Jan
  sel <- which(days == 0)
  days[sel] <- 31
  months[sel] <- 1

  # 29 of Feb as 1 Mar
  sel <- which(months == 2 & days == 29)
  months[sel] <- 3
  days[sel] <- 1

  # Construct new date variable
  real_dates<-as.POSIXlt(ISOdate(years, months, days))
  real_dates<-trunc(real_dates,"day")
  return(real_dates)
}
#' convert downloaded ukcp18 nc from 360 day to 365 day
.fill_calendar_data<-function(ukcp_r, real_dates, testplot=FALSE, plotdays=c(89:91,242:244)){
  # Assign real dates to layer names and time values
  terra::time(ukcp_r)<-real_dates
  names(ukcp_r)<-real_dates
  # Insert and fill missing date layers of spatrast - could use zoo::na.approx or na.spline
  ukcp_r<-terra::fillTime(ukcp_r)
  if(testplot)  plot(ukcp_r[[plotdays]],main=paste(time(ukcp_r)[plotdays], 'before filling'))
  ukcp_r<-terra::approximate(ukcp_r,method="linear")
  if(testplot)  plot(ukcp_r[[plotdays]],main=paste(time(ukcp_r)[plotdays], 'after filling'))
  return(ukcp_r)
}
#' apply bias correction to UKCP18 data
.biascorrect <- function(hist_obs, hist_mod, fut_mod, rangelims = 1.05, silent = FALSE) {
  if (!inherits(hist_obs, "SpatRaster")) stop("hist_obs must be a SpatRaster")
  if (!inherits(hist_mod, "SpatRaster")) stop("hist_mod must be a SpatRaster")
  if (!inherits(fut_mod, "SpatRaster")) stop("hist_mod must be a SpatRaster")
  # reproject and crop if necessary
  if (crs(hist_mod) != crs(hist_obs)) hist_mod<-project(hist_mod,hist_obs)
  if (crs(fut_mod) != crs(hist_obs)) fut_mod<-project(fut_mod,hist_obs)
  hist_mod<-resample(hist_mod,hist_obs)
  fut_mod<-resample(fut_mod,hist_obs)
  # Convert to arrays
  a1<-.is(hist_obs)
  a2<-.is(hist_mod)
  a3<-.is(fut_mod)
  # Mask out any cells that are missing
  msk1<-apply(a1,c(1,2),mean,na.rm=T)
  msk2<-apply(a2,c(1,2),mean,na.rm=T)
  msk3<-apply(a3,c(1,2),mean,na.rm=T)
  msk<-msk1*msk2*msk3
  # Check whether dataset has more than 1000 entries per time-series and
  # subset if it does
  n<-dim(a1)[3]
  if (n > 1000) s<-sample(0:n,998,replace = FALSE)  # 998 so min and max can be tagged on
  ao<-array(NA,dim=dim(a3))
  # Create array for storing data
  counter<-0
  nn<-dim(a1)[1]*dim(a1)[2]
  if (silent == FALSE) pb<-utils::txtProgressBar(min = 0, max = nn, style = 3)
  for (i in 1:dim(a1)[1]) {
    for (j in 1:dim(a1)[2]) {
      if (silent == FALSE) utils::setTxtProgressBar(pb,counter)
      if (is.na(msk[i,j]) == FALSE) {
        v1<-a1[i,j,]
        v2<-a2[i,j,]
        v1 <- v1[order(v1)]
        v2 <- v2[order(v2)]
        if (n > 1000) {
          v1<-c(v1[1],v1[s],v1[length(v1)]) # tags on min and max value
          v2<-c(v2[1],v2[s],v2[length(v2)]) # tags on min and max value
        }
        # Apply gam
        m1 <- gam(v1~s(v2))
        v3<-a3[i,j,]
        xx <- as.numeric(predict.gam(m1, newdata = data.frame(v2 = v3)))
        if (is.na(rangelims) == FALSE) xx<-rangelimapply(v1, v2, v3, xx)
        ao[i,j,]<-xx
      }
      counter<-counter+1
    }
  }
  ao<-.rast(ao,hist_obs)
  return(ao)
}
#' @title Calculate moving average
#' @noRd
.mav <- function(x, n = 5) {
  y <- stats::filter(x, rep(1 / n, n), circular = TRUE, sides = 1)
  y
}
#' Correct precipitation
.precipcorrect <- function(hist_obs, hist_mod, fut_mod, rangelim = 1.05) {
  if (!inherits(hist_obs, "SpatRaster")) stop("hist_obs must be a SpatRaster")
  if (!inherits(hist_mod, "SpatRaster")) stop("hist_mod must be a SpatRaster")
  if (!inherits(fut_mod, "SpatRaster")) stop("hist_mod must be a SpatRaster")
  # reproject and crop if necessary
  if (crs(hist_mod) != crs(hist_obs)) hist_mod<-project(hist_mod,hist_obs)
  if (crs(fut_mod) != crs(hist_obs)) fut_mod<-project(fut_mod,hist_obs)
  hist_mod<-resample(hist_mod,hist_obs)
  fut_mod<-resample(fut_mod,hist_obs)
  # mask data
  hist_mod<-mask(hist_mod,hist_obs)
  # Calculate observed rainfall total and rain day frac
  rcount<-hist_obs
  rcount[rcount > 0] <-1
  rtot1<-apply(.is(hist_obs),c(1,2),sum)
  tfrac1<-apply(.is(rcount),c(1,2),sum)/dim(rcount)[3]
  # Calculate observed rainfall total and rain day frac
  rcount<-hist_mod
  rcount[rcount > 0] <-1
  rtot2<-apply(.is(hist_mod),c(1,2),sum)
  tfrac2<-apply(.is(rcount),c(1,2),sum)/dim(rcount)[3]
  # Calculate ratios
  mu_tot<-rtot1/rtot2
  mu_frac<-tfrac1/tfrac2
  # Apply correction to modelled future data
  rcount<-fut_mod
  rcount[rcount > 0] <-1
  rtot<-apply(.is(fut_mod),c(1,2),sum)*mu_tot
  tfrac<-(apply(.is(rcount),c(1,2),sum)/dim(rcount)[3])*mu_frac
  # Calculate regional precipitation
  rrain<-apply(.is(fut_mod),3,sum,na.rm=TRUE)
  rr2<-as.numeric(.mav(rrain,10))
  s<-which(rrain==0)
  rrain[s]<-rr2[s]
  s<-which(rrain==0)
  rrain[s]<-0.1
  rrain<-(rrain/max(rrain))*10
  # Adjust rainfall
  fut_mod<-mask(fut_mod,hist_obs)
  a<-as.array(fut_mod)
  mm<-matrix(as.vector(a),nrow=dim(a)[1]*dim(a)[2],ncol=dim(a)[3])
  rtot<-.rast(rtot,hist_obs)
  tfrac<-.rast(tfrac,hist_obs)
  rtot<-as.vector(t(rtot))
  rfrac<-as.vector(t(tfrac))
  mm<-rainadjustm(mm,rrain,rfrac,rtot)
  a2<-array(mm,dim=dim(a))
  a2<-prangelimapply(.is(hist_obs), .is(hist_mod), .is(fut_mod), a2, rangelims)
  out<-.rast(a2,fut_mod)
  return(out)
}
#' Used by albedo calculations - calculates wavelength specific spectral density
.Planck <- function(wavelength, temperature =  5504.85) {
  d <- wavelength * 1e-09
  h <- 6.6256e-34
  cc <- 299792458
  tt <- temperature + 273.15
  k <- 1.38054e-23
  b <- (2 * pi * h * cc^2) / (d^5 * (exp((h * cc) / (k * d * tt)) - 1))
  b
}
#' sharpens a coarse resolution raster using a fine-resolution categorical raster
.sharpen <- function(coarse, fine, msk = TRUE) {
  if (msk) {
    mskc <- resample(fine, coarse)
    coarse <- mask(coarse, mskc)
  }
  # Produce table of categorical fine scale variable
  ufine <- unique(as.vector(fine))
  ufine <- ufine[is.na(ufine) == FALSE]
  # Calculate the mean value of coarse for each unique value of fine
  mcoarse <- 0
  for (i in 1:length(ufine)) {
    # Calculate fraction of each fine type in each coarse grid cell
    ftype <- fine
    ftype[ftype != ufine[i]] <-0
    ftype[ftype == ufine[i]] <-1
    cfrac <- resample(ftype, coarse)
    # multiply coarse by fraction and sum
    mucoarse <- cfrac * coarse
    fsum <- sum(as.vector(mucoarse), na.rm=TRUE)
    # Sum the fractions across cells
    csum <- sum(as.vector(cfrac), na.rm=TRUE)
    # Calculate the mean coarse for each fine type
    mcoarse[i] <- fsum / csum
  }
  mx <- (max(mcoarse) / mean(mcoarse)) * 3
  # Calculate a fine-scale dataset of expected values of coarse
  finem <- as.matrix(fine, wide = TRUE)
  evals <- finem * 0
  for (i in 1:length(ufine)) {
    s <- which(finem == ufine[i])
    evals[s] <- mcoarse[i]
  }
  evals <- .rast(evals, fine)
  # Calculate multiplier by dividing actual by expected (coarse scale)
  mu <- coarse / resample(evals, coarse)
  mskc <- mskc * 0 +1
  mskc[is.na(mskc)] <- 1
  mu <- .fillna(mu, mskc)
  mu <- resample(mu, fine)
  mu[mu > mx] <- mx
  sharpened <- mu * evals
  return(sharpened)
}

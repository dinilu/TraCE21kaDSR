
.elementwise.all.equal <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})

#' Get vertical levels in a grid.
#'
#' @param grid Grid object (see loadeR package) from which to get the levels
#' @param level Level to be found within the levels in the grid object
#'
#' @return List object with two vectors inside: 
#' # level: the level searched.
#' # zRange: the position for the level.
#' @export
#'
#' @examples
#' TBW
getVerticalLevelPars <- function (grid, level) {
  gcs <- grid$getCoordinateSystem()
  if (gcs$hasVerticalAxis()) {
    levels <- loadeR::scanVarDimensions(grid)$level$Values
    if (is.null(level)) {
      if (length(levels) == 1) {
        level <- levels
        if (gcs$getVerticalAxis()$findCoordElement(level) < 
            0) {
          levelInd <- gcs$getVerticalAxis()$findCoordElement(0)
        }
      } else {
        stop("Variable with vertical levels: '@level' following the variable name is required\nPossible values: ", 
             paste(levels, collapse = ", "))
      }
    } else {
      if (any(.elementwise.all.equal(level, levels))) {
        levelInd <- gcs$getVerticalAxis()$findCoordElement(level)
      } else {
        stop("Vertical level not found\nPossible values: ", 
             paste(levels, collapse = ", "), call. = FALSE)
      }
    }
    zRange <- .jnew("ucar/ma2/Range", levelInd, levelInd)
  } else {
    if (!is.null(level)) {
      level <- level
    }
    zRange <- .jnull()
  }
  return(list(level = level, zRange = zRange))
}

utils::assignInNamespace("getVerticalLevelPars", getVerticalLevelPars, ns="loadeR")


#' Copy coordinates info from one grid to another
#'
#' @param x A grid object (see loadeR package) to be modified.
#' @param y A grid object from which to extract the coordinates info.
#'
#' @return The x grid object with the coordinates info from the y grid object.
#' @export
#'
#' @examples
copyXYCoords <- function(x, y){
  x$xyCoords <- y$xyCoords
  return(x)
} 


#' Compute wind speed from their horizontal (u) and vertical (v) components.
#'
#' @param u A grid object (see loadeR package) 
#' @param v 
#'
#' @return
#' @export
#'
#' @examples
compute_wind_speed <- function(u, v) { 
  message("[", Sys.time(), "] Computing wind speed (wss) from its horizontal (u) and vertical (v) components")
  u <- transformeR::gridArithmetics(u, u, operator="*")
  v <- transformeR::gridArithmetics(v, v, operator="*")
  
  ws <- transformeR::gridArithmetics(u, v, operator="+")
  
  ws$Data <- round(sqrt(ws$Data), 1)
  ws$Variable$varName <- "wss"
  ws$Variable$level <- NA
  attr(ws$Variable, "use_dictionary") <- TRUE
  attr(ws$Variable, "description") <- "Wind speed at surface level"
  attr(ws$Variable, "units") <- "m/s"
  attr(ws$Variable, "longname") <- "wind speed"
  
  return(ws)
}

# Load custom functions
DateSeq <- function(st, en, freq, frac=0) {
  require(zoo)
  st <- as.Date(as.yearmon(st)) 
  en <- as.Date(as.yearmon(en)) 
  result <- as.Date(as.yearmon(seq(st, en, by = paste(as.character(12/freq), "months"))), frac = frac)
  return(result)
}

loadTraceData <- function(file, var=NULL, lonLim=trace.lon, latLim=trace.lat, start_date="1551-01-01", end_date="1990-12-31", years=1961:1990, dictionary="../../Data/Trace21ka/dictionary.dic") {
  require(loadeR)
  require(transformeR)
  if(is.null(var)){
    stop("Argument var not defined. Please define a variable to be loaded.")
  }
  data <- loadeR::loadGridData(dataset = file,
                               var = var,
                               lonLim = lonLim,
                               latLim = latLim, 
                               dictionary = dictionary)
  data$Dates$start <- paste(DateSeq(start_date, end_date, 12, 0), "00:00:00 GMT", sep = " ")
  data$Dates$end <- paste(DateSeq(start_date, end_date, 12, 1), "00:00:00 GMT", sep = " ")
  
  data <- transformeR::subsetGrid(data, years=years)  
  
  return(data)
}

loadTrace <- function(file_list, var_list, lonLim=trace.lon, latLim=trace.lat, start_date="1551-01-01", end_date="1990-12-31", years=1961:1990, dictionary="../../Data/Trace21ka/dictionary.dic", var_selection=trace.model.var.names, compute_wss=TRUE){
  
  data <- mapply(loadTraceData, file=file_list, var=var_list, MoreArgs=list(lonLim= lonLim, latLim = latLim, years=years, dictionary=dictionary), SIMPLIFY=FALSE)
  
  names(data) <- var_list
  
  if(compute_wss == TRUE){
    data$wss <- compute_wind_speed(data$'u@992.5561', data$'v@992.5561')
  } 
  
  data <- transformeR::makeMultiGrid(data)
  data <- transformeR::subsetGrid(data, var=var_selection)
  
  return(data)
} 


loadCMIP <- function(var_list, var_new_list, indir = "../../Data/CMIP5/", rcp, mod, lonLim=cmip5.lon, latLim=cmip5.lat, years=NULL, dictionary="../../Data/CMIP5/dictionary.dic"){
  # var_list <- cmip5.vars
  # var_new_list <- cmip5.new.vars
  # indir <- "../../Data/CMIP5/"
  # rcp <- "rcp2.6"
  # mod <- cmip5.mods[[1]]
  # lonLim <- cmip5.lon
  # latLim <- cmip5.lat
  # years <- 2006:2100
  # dictionary <- "../../Data/CMIP5/dictionary.dic"
  
  if(rcp != "historical"){
    rcp_folder <- paste0("RCP_", substr(rcp, 4, 4), "_", substr(rcp, 6, 6))
    rcp_file <- paste0("rcp", substr(rcp, 4, 4), substr(rcp, 6, 6))
    dates <- "200601-210012"
  } else {
    rcp_folder <- rcp
    rcp_file <- rcp
    dates <- "185001-200512"
  } 
  
  file_list <- paste0(indir, rcp, "/", mod, "_", rcp_folder, "/", var_list, "_Amon_", mod, "_", rcp_file, "_r1i1p1_", dates, ".nc")
  
  data <- mapply(loadeR::loadGridData, dataset = file_list, var = var_new_list, MoreArgs=list(lonLim = lonLim, latLim = latLim, dictionary = dictionary), SIMPLIFY = FALSE)

  start_date <- paste0(years[[1]], "-01-01")
  end_date <- paste0(years[[length(years)]], "-12-31")
  
  data <- lapply(data, FUN = modifyDates, start_date, end_date)
  
  names(data) <- var_new_list
  
  data <- transformeR::makeMultiGrid(data)
  
  if(!is.null(years)){
    data <- transformeR::subsetGrid(data, years = years)
  } 
  
  return(data)
} 


modifyDates <- function(x, start_date="1961-01-01", end_date="1990-12-31") {
  
  x$Dates$start <- paste(DateSeq(start_date, end_date, 12, 0), "00:00:00 GMT", sep = " ")
  
  x$Dates$end <- paste(DateSeq(start_date, end_date, 12, 1), "00:00:00 GMT", sep = " ")
  
  return(x)
}

loadUerra <- function(file, var, lonLim = uerra.lon, latLim = uerra.lat, dictionary = "../../Data/UERRA/UERRA_dictionary.dic", new.coordinates){
  
  data <- loadeR::loadGridData(file, var = var, lonLim = lonLim, latLim = latLim, dictionary = dictionary)
  
  if(var == "pr"){
    data <- transformeR::upscaleGrid(data, times=2, aggr.fun=list(FUN=mean))
    
    data <- transformeR::interpGrid(data, new.coordinates = new.coordinates, method="bilinear")
  }
  
  data <- modifyDates(data)
  
  return(data)
}



loadManualTraceData <- function(file, var, trace.y1, trace.y2, lonLim = trace.lon, latLim = trace.lat, dictionary="../../Data/Trace21ka/dictionary.dic"){ 
  
  var <- loadeR::findVerticalLevel(var)
  
  dic <- loadeR:::dictionaryLookup(dictionary, var$var, "none")
  vocabulary <- climate4R.UDG::C4R.vocabulary()

  trace.nc <- ncdf4::nc_open(file)

  long_name <- ncdf4::ncatt_get(trace.nc, dic$short_name)$long_name
  
  if(!is.null(dictionary)){
    units <- as.character(vocabulary[grep(paste0("^", var$var, "$"), vocabulary$identifier), 3])
  }else{
    units <- ncdf4::ncatt_get(trace.nc, dic$short_name)$units
  } 
  
  trace.c4r <- list()
  trace.c4r[[1]] <- list()
  if(is.null(var$level)){
    trace.c4r[[1]][1]  <- var$var
  }else{
    trace.c4r[[1]][1]  <- paste0(var$var, "@", var$level) 
  } 
  trace.c4r[[1]][2]   <- list(NULL)
  names(trace.c4r[[1]]) <- c("varName", "level")
  attr(trace.c4r[[1]], "use_dictionary") <- TRUE
  attr(trace.c4r[[1]], "description") <- long_name
  attr(trace.c4r[[1]], "units") <- units
  attr(trace.c4r[[1]], "longname") <- var$var
  attr(trace.c4r[[1]], "daily_agg_cellfun") <- "none"
  attr(trace.c4r[[1]], "monthly_agg_cellfun") <- "none"
  attr(trace.c4r[[1]], "verification_time") <- "none"
  
  if(!is.null(var$level)){
    lev <- ncdf4::ncvar_get(trace.nc, "lev")
    lev <- which(.elementwise.all.equal(lev, var$level))
    trace.var <- ncdf4::ncvar_get(trace.nc, dic$short_name, c(1,1,lev,1))
  }else{
    trace.var <- ncdf4::ncvar_get(trace.nc, dic$short_name)
  }  

  trace.c4r[[2]] <- aperm(trace.var, c(3,2,1)) * dic$scale + dic$offset
  
  trace.c4r[[3]] <- list()
  # trace.c4r[[3]]$x <- ncvar_get(trace.nc, "lon")
  lon <- ncdf4::ncvar_get(trace.nc, "lon")
  lon[which(lon>180)] <- lon[which(lon>180)] - 360 
  trace.c4r[[3]]$x <- lon
  trace.c4r[[3]]$y <- ncdf4::ncvar_get(trace.nc, "lat")
  attr(trace.c4r[[3]], "projection") <- "LatLonProjection"
  attr(trace.c4r[[3]], "resX") <- (max(trace.c4r[[3]]$x) - min(trace.c4r[[3]]$x)) / (length(trace.c4r[[3]]$x) - 1)
  attr(trace.c4r[[3]], "resY") <- (max(trace.c4r[[3]]$y) - min(trace.c4r[[3]]$y)) / (length(trace.c4r[[3]]$y) - 1)
  
  trace.c4r[[2]] <- trace.c4r[[2]][,,order(trace.c4r[[3]]$x)]
  trace.c4r[[3]]$x <- sort(trace.c4r[[3]]$x)
  attr(trace.c4r[[2]], "dimensions") <- c("time", "lat", "lon")
  
  if(0 %in% trace.y1:trace.y2){
    n.years <- length(trace.y1:trace.y2) - 1    
  }else{
    n.years <- length(trace.y1:trace.y2)
  }  
  trace.c4r[[4]] <- list() 
  trace.c4r[[4]]$start <- paste(DateSeq("4000-01-01", paste0(3999+n.years, "-12-31"), 12, 0), "00:00:00 GMT", sep = " ")
  trace.c4r[[4]]$end <- paste(DateSeq("4000-01-01", paste0(3999+n.years, "-12-31"), 12, 1), "00:00:00 GMT", sep = " ")
  attr(trace.c4r[[4]], "subset") <- "subsetYears"
  attr(trace.c4r[[4]], "season") <- 1:12 

  names(trace.c4r) <- c("Variable", "Data", "xyCoords", "Dates")
  
  attr(trace.c4r, "dataset") <- file
  attr(trace.c4r, "R_package_desc") <- paste0("loadeR-v", packageVersion("loadeR"))
  attr(trace.c4r, "R_package_URL") <- "https://github.com/SantanderMetGroup/loadeR"
  attr(trace.c4r, "R_package_ref") <- "https://doi.org/10.1016/j.envsoft.2018.09.009"
  
  if(!is.null(lonLim) & !is.null(latLim)){
    trace.c4r <- transformeR::subsetGrid(trace.c4r, lonLim = lonLim, latLim = latLim)
  } 
  trace.c4r <- recalcGridResolution(trace.c4r)
  
  return(trace.c4r)
}

recalcGridResolution <- function(grid){
  attr(grid$xyCoords, "resX") <- (max(grid$xyCoords$x) - min(grid$xyCoords$x)) / (length(grid$xyCoords$x) - 1)
  attr(grid$xyCoords, "resY") <- (max(grid$xyCoords$y) - min(grid$xyCoords$y)) / (length(grid$xyCoords$y) - 1)
  return(grid)  
} 

loadManualTrace <- function(file_list, var_list, trace.y1, trace.y2, lonLim, latLim, dictionary="../../Data/Trace21ka/dictionary.dic", var_selection=trace.model.var.names, compute_wss=TRUE){
  
  data <- mapply(loadManualTraceData, file=file_list, var=var_list, MoreArgs=list(trace.y1, trace.y2, lonLim, latLim, dictionary), SIMPLIFY=FALSE)
  
  names(data) <- var_list
  
  if(compute_wss == TRUE){
    data$wss <- compute_wind_speed(data$'u@992.5561', data$'v@992.5561')
  } 
  
  data <- transformeR::makeMultiGrid(data)
  data <- transformeR::subsetGrid(data, var=var_selection)
  
  return(data)
} 

nc2sp_df <- function(grid, output.dir){
  sp <- transformeR::grid2sp(grid)
  df <- as.data.frame(sp)
  df <- df[,c(13,14,1:12)] 
  colnames(df) <- c("x", "y", 1:12)
  df
}


downscaleTrace <- function(i, new.data.list, var.names, y1.list, y2.list, lonLim, latLim, hist.trace, data, model, local.var, trace.model.var.names, global.nc.attributes){

  # i <- 36
  # new.data.list <- new.trace.file.names
  # var.names <- trace.var.names
  # y1.list <- years.y1
  # y2.list <- years.y2
  # lonLim <- trace.lon
  # latLim <- trace.lat
  # hist.trace <- hist.trace
  # data <- data
  # model <- model
  # local.var <- local.var
  # trace.model.var.names <- trace.model.var.names
  # global.nc.attributes <- global.nc.attributes

  if(!dir.exists(paste0("../../Output/Trace21ka/", local.var, "/dat"))){
    dir.create(paste0("../../Output/Trace21ka/", local.var, "/dat"), recursive = TRUE)
  } 
  
  new.data <- new.data.list[[i]] 
  y1 <- y1.list[[i]] 
  y2 <- y2.list[[i]] 
  
  new.trace <- loadManualTrace(new.data, var.names, y1, y2, lonLim, latLim, var_selection = trace.model.var.names)

  new.trace.xy <- copyXYCoords(new.trace, hist.trace)
  
  if(y1 == 400){
    real.years <- c(-y1:-1,1:-y2)
  }else{
    real.years <- -y1:-y2
  }  

  fake.years <- seq(4000, length=length(real.years))
  
  ydiff <- -y1 - 3999  
  
  for(j in 1:length(fake.years)){
  
    message("Calculating year: ", real.years[j], "...")
    
    new.trace.sub <- transformeR::subsetGrid(new.trace.xy, years = fake.years[j])

    new.data <- downscaleR::prepareNewData(new.trace.sub, data)
    
    pred <- downscaleR::downscalePredict(new.data, model)
    
    pred$Data <- round(pred$Data, 2)

    loadeR.2nc::grid2nc(pred, NetCDFOutFile = paste0("../../Output/Trace21ka/", local.var, "/", local.var, real.years[j], "_tmp.nc"), missval = -9999, globalAttributes = global.nc.attributes)

    infile <- paste0("../../Output/Trace21ka/", local.var, "/", local.var, real.years[j],  "_tmp.nc")
    outfile <- paste0("../../Output/Trace21ka/", local.var, "/", local.var, real.years[j], ".nc")

    message("New year: ", fake.years[j] + ydiff)
    
    # system(paste0("cdo -r setreftime,1950-01-01,00:00:00,hours -shifttime,", ydiff, "y ", infile, " ", outfile)) #not working  
    system2("cdo", c("-r", "setreftime,1950-01-01,00:00:00,hours", paste0("-shifttime,", ydiff, "y"), infile, outfile))
    file.remove(infile)
    
    pred.df <- nc2sp_df(pred)
    
    write.table(pred.df, file=paste0("../../Output/Trace21ka/", local.var, "/dat/", local.var, real.years[j], ".dat"), sep="\t")

    message("   ...done")
  } 
  return("Done")
} 


downscaleTraceBimodel <- function (i, new.data.list, var.names, y1.list, y2.list, lonLim, latLim, hist.trace, data, model, data.bin, model.bin, local.var, trace.model.var.names, global.nc.attributes){
  
  # i <- 36
  # new.data.list <- new.trace.file.names
  # var.names <- trace.var.names
  # y1.list <- years.y1
  # y2.list <- years.y2
  # lonLim <- trace.lon
  # latLim <- trace.lat
  # hist.trace <- hist.trace
  # data <- data
  # model <- model
  # data.bin <- data_bin
  # model.bin <- model_bin
  # local.var <- local.var
  # trace.model.var.names <- trace.model.var.names
  # global.nc.attributes <- global.nc.attributes
  
  if(!dir.exists(paste0("../../Output/Trace21ka/", local.var, "/dat"))){
    dir.create(paste0("../../Output/Trace21ka/", local.var, "/dat"))
  } 
  
  new.data <- new.data.list[[i]] 
  y1 <- y1.list[[i]] 
  y2 <- y2.list[[i]] 
  
  new.trace <- loadManualTrace(new.data, var.names, y1, y2, lonLim, latLim, var_selection = trace.model.var.names)
  
  new.trace.xy <- copyXYCoords(new.trace, hist.trace)
  
  if(y1 == 400){
    real.years <- c(-y1:-1,1:-y2)
  }else{
    real.years <- -y1:-y2
  }  
  
  fake.years <- seq(4000, length=length(real.years))
  
  for(j in 1:length(fake.years)){
    
    if(real.years[j] > 0){
      ydiff <- -y1 - 4000  
    }else{
      ydiff <- -y1 - 3999
    } 
    
    message("Calculating year: ", real.years[j], "...")
    
    new.trace.sub <- transformeR::subsetGrid(new.trace.xy, years = fake.years[j])
    
    new.data <- downscaleR::prepareNewData(new.trace.sub, data)
    
    pred.bin <- downscaleR::downscalePredict(new.data, model.bin)
    pred.cont <- downscaleR::downscalePredict(new.data, model)
    
    pred <- transformeR::gridArithmetics(pred.bin, pred.cont, operator = "*")
    
    pred$Data <- round(pred$Data, 2)
    
    loadeR.2nc::grid2nc(pred, NetCDFOutFile = paste0("../../Output/Trace21ka/", local.var, "/", local.var, real.years[j], "_tmp.nc"), missval = -9999, globalAttributes = global.nc.attributes)
    
    infile <- paste0("../../Output/Trace21ka/", local.var, "/", local.var, real.years[j],  "_tmp.nc")
    outfile <- paste0("../../Output/Trace21ka/", local.var, "/", local.var, real.years[j], ".nc")
    
    message("New year: ", fake.years[[j]] + ydiff)
    
    # system(paste0("cdo -r setreftime,1950-01-01,00:00:00,hours -shifttime,", ydiff, "y ", infile, " ", outfile))
    system2("cdo", c("-r", "setreftime,1950-01-01,00:00:00,hours", paste0("-shifttime,", ydiff, "y"), infile, outfile))
    file.remove(infile)
    
    pred.df <- nc2sp_df(pred)
    
    write.table(pred.df, file=paste0("../../Output/Trace21ka/", local.var, "/dat/", local.var, real.years[j], ".dat"), sep="\t")
    
    message("   ...done")
  } 
  return("Done")
} 


downscaleCMIP5 <- function(uerra, var_list, var_new_list, rcp, mod, lonLim = cmip5.lon, latLim = cmip5.lat, dictionary = "../../Data/CMIP5/dictionary.dic", indir = "../../Data/CMIP5/", outdir = "../../Output/CMIP5/", local.var, spatial.pars = cmip5.spatial.pars, method = "GLM", family.link = family.link, global.nc.attributes = global.nc.attributes){

  # uerra <- uerra
  # var_list <- cmip5.vars
  # var_new_list <- cmip5.new.vars
  # rcp <- "rcp6.0"
  # mod <- cmip5.mods[[3]]
  # lonLim <- cmip5.lon
  # latLim <- cmip5.lat
  # dictionary <- "../../Data/CMIP5/dictionary.dic"
  # indir <- "../../Data/CMIP5/"
  # outdir <- "../../Output/CMIP5/"
  # local.var <- "cld"
  # spatial.pars = cmip5.spatial.pars
  # family.link <- "gaussian"
  # global.nc.attributes <- global.nc.attributes

  if(!dir.exists(paste0(outdir, rcp, "/", mod, "/", local.var, "/dat"))){
    dir.create(paste0(outdir, rcp, "/", mod, "/", local.var, "/dat"), recursive = TRUE)
  } 
  
  hist.cmip5 <- loadCMIP(var_list = var_list, var_new_list = var_new_list, indir = indir, rcp = "historical", mod = mod, years=1961:1990)

  if(local.var == "pr"){
    uerra.bin <- transformeR::binaryGrid(uerra, condition = "GE", threshold = 1)
    
    data.bin <- downscaleR::prepareData(hist.cmip5, uerra.bin, spatial.predictors = cmip5.spatial.pars)
    data <- downscaleR::prepareData(hist.cmip5, uerra, spatial.predictors = cmip5.spatial.pars)
    
    model.bin <- downscaleR::downscaleTrain(data.bin, method = "GLM", family = binomial(link="logit"), predict = TRUE)
    model <- downscaleR::downscaleTrain(data, method = "GLM", family = family.link, predict = TRUE, condition = "GE", threshold = 1)
  } else {
    data <- downscaleR::prepareData(hist.cmip5, uerra, spatial.predictors = cmip5.spatial.pars)
    model <- downscaleR::downscaleTrain(data, method = method, family = family.link, predict = TRUE)
  } 

  rcp.cmip5.1 <- loadCMIP(var_list = var_list, var_new_list = var_new_list, indir = indir, rcp = "historical", mod = mod, years=1991:2005)
  rcp.cmip5.2 <- loadCMIP(var_list = var_list, var_new_list = var_new_list, indir = indir, rcp = rcp, mod = mod, years = 2006:2100)
  
  rcp.cmip5 <- transformeR::bindGrid(rcp.cmip5.1, rcp.cmip5.2, dimension = "time")

  new.data <- downscaleR::prepareNewData(rcp.cmip5, data)
  
  if( local.var == "pr"){
    pred.bin <- downscaleR::downscalePredict(new.data, model.bin)
    pred.cont <- downscaleR::downscalePredict(new.data, model)
    pred <- transformeR::gridArithmetics(pred.bin, pred.cont, operator = "*")
  } else {
    pred <- downscaleR::downscalePredict(new.data, model)
  } 

  pred$Data <- round(pred$Data, 2)
  
  loadeR.2nc::grid2nc(pred, NetCDFOutFile = paste0(outdir, rcp, "/", mod, "/", local.var, "/", local.var, "1991-2100.nc"), missval = -9999, globalAttributes = global.nc.attributes)
    
  for(i in 1:110){
    y <- c(1991:2100)[i] 
    yBP <- c(41:151)[i] 
    pred.i <- transformeR::subsetGrid(pred, years=y)
    
    pred.df <- nc2sp_df(pred.i)
    
    file.name <- paste0(outdir, rcp, "/", mod, "/", local.var, "/dat/", local.var, yBP, ".dat")
    write.table(pred.df, file = file.name, sep = "\t")
  } 
  return("Done")
} 

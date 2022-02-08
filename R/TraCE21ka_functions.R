
#' Title
#'
#' @param indir TBW
#' @param vars TBW
#' @param file.number TBW 
#'
#' @return TBW
#' @export
#'
#' @examples  #TBW
traceFileNames <- function(indir, vars = NULL, file.number = 1){
  # indir <- "../../Data/TraCE21ka/"
  if(is.null(vars)){
    vars <- trace.vars
  }
  years.n <- trace.years.n[file.number] 
  years.bp <- trace.years.bp[file.number]
  years.nums <- trace.years.nums[file.number]
  paste0(indir, vars, "/trace.", years.n, ".", years.bp, ".cam2.h0.", vars, ".", years.nums, ".nc")
} 

#' Title
#'
#' @param file TBW
#' @param var  TBW
#' @param lonLim  TBW
#' @param latLim  TBW
#' @param start_date  TBW
#' @param end_date  TBW
#' @param years  TBW
#' @param dictionary TBW 
#'
#' @return TBW
#' @export
#'
#' @examples #TBW
loadTraceData <- function(file, var=NULL, lonLim=trace.lon, latLim=trace.lat, start_date="1551-01-01", end_date="1990-12-31", years=1961:1990, dictionary=system.file("extdata", "TraCE21ka_dictionary.csv", package = "TraCE21kaDSR")) {
  if(is.null(var)){
    stop("Argument var not defined. Please define a variable to be loaded.")
  }
  data <- loadeR::loadGridData(dataset = file,
                               var = var,
                               lonLim = lonLim,
                               latLim = latLim, 
                               dictionary = dictionary)
  data$Dates$start <- paste(.DateSeq(start_date, end_date, 12, 0), "00:00:00 GMT", sep = " ")
  data$Dates$end <- paste(.DateSeq(start_date, end_date, 12, 1), "00:00:00 GMT", sep = " ")
  
  data <- transformeR::subsetGrid(data, years=years)  
  
  return(data)
}


#' Title
#'
#' @param file_list TBW
#' @param var_list TBW
#' @param lonLim TBW
#' @param latLim TBW
#' @param start_date TBW
#' @param end_date TBW
#' @param years TBW
#' @param dictionary TBW
#' @param var_selection TBW
#' @param compute_wss TBW
#'
#' @return TBW
#' @export
#'
#' @examples #TBW
loadTrace <- function(file_list, var_list, lonLim=trace.lon, latLim=trace.lat, start_date="1551-01-01", end_date="1990-12-31", years=1961:1990, dictionary=system.file("extdata", "TraCE21ka_dictionary.csv", package = "TraCE21kaDSR"), var_selection='trace.model.var.names', compute_wss=TRUE){
  
  data <- mapply(loadTraceData, file=file_list, var=var_list, MoreArgs=list(lonLim= lonLim, latLim = latLim, years=years, dictionary=dictionary), SIMPLIFY=FALSE)
  
  names(data) <- var_list
  
  if(compute_wss == TRUE){
    data$wss <- compute_wind_speed(data$'u@992.5561', data$'v@992.5561')
  } 
  
  data <- transformeR::makeMultiGrid(data)
  data <- transformeR::subsetGrid(data, var=var_selection)
  
  return(data)
} 



#' Title
#'
#' @param file TBW
#' @param var TBW
#' @param trace.y1 TBW
#' @param trace.y2 TBW
#' @param lonLim TBW
#' @param latLim TBW
#' @param dictionary TBW
#'
#' @return TBW
#' @export
#'
#' @examples #TBW
loadManualTraceData <- function(file, var, trace.y1, trace.y2, lonLim = trace.lon, latLim = trace.lat, dictionary=system.file("extdata", "TraCE21ka_dictionary.csv", package = "TraCE21kaDSR")){ 
  
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
  trace.c4r[[4]]$start <- paste(.DateSeq("4000-01-01", paste0(3999+n.years, "-12-31"), 12, 0), "00:00:00 GMT", sep = " ")
  trace.c4r[[4]]$end <- paste(.DateSeq("4000-01-01", paste0(3999+n.years, "-12-31"), 12, 1), "00:00:00 GMT", sep = " ")
  attr(trace.c4r[[4]], "subset") <- "subsetYears"
  attr(trace.c4r[[4]], "season") <- 1:12 
  
  names(trace.c4r) <- c("Variable", "Data", "xyCoords", "Dates")
  
  attr(trace.c4r, "dataset") <- file
  attr(trace.c4r, "R_package_desc") <- paste0("loadeR-v", utils::packageVersion("loadeR"))
  attr(trace.c4r, "R_package_URL") <- "https://github.com/SantanderMetGroup/loadeR"
  attr(trace.c4r, "R_package_ref") <- "https://doi.org/10.1016/j.envsoft.2018.09.009"
  
  if(!is.null(lonLim) & !is.null(latLim)){
    trace.c4r <- transformeR::subsetGrid(trace.c4r, lonLim = lonLim, latLim = latLim)
  } 
  trace.c4r <- recalcGridResolution(trace.c4r)
  
  return(trace.c4r)
}


#' Title
#'
#' @param file_list TBW
#' @param var_list TBW
#' @param trace.y1 TBW
#' @param trace.y2 TBW
#' @param lonLim TBW
#' @param latLim TBW
#' @param dictionary TBW
#' @param var_selection TBW
#' @param compute_wss TBW
#'
#' @return TBW
#' @export
#'
#' @examples #TBW
loadManualTrace <- function(file_list, var_list, trace.y1, trace.y2, lonLim, latLim, dictionary=system.file("extdata", "TraCE21ka_dictionary.csv", package = "TraCE21kaDSR"), var_selection='trace.model.var.names', compute_wss=TRUE){
  
  data <- mapply(loadManualTraceData, file=file_list, var=var_list, MoreArgs=list(trace.y1, trace.y2, lonLim, latLim, dictionary), SIMPLIFY=FALSE)
  
  names(data) <- var_list
  
  if(compute_wss == TRUE){
    data$wss <- compute_wind_speed(data$'u@992.5561', data$'v@992.5561')
  } 
  
  data <- transformeR::makeMultiGrid(data)
  data <- transformeR::subsetGrid(data, var=var_selection)
  
  return(data)
} 


#' Title
#'
#' @param i TBW
#' @param new.data.list TBW
#' @param var.names TBW
#' @param y1.list TBW
#' @param y2.list TBW
#' @param lonLim TBW
#' @param latLim TBW
#' @param hist.trace TBW
#' @param data TBW
#' @param model TBW
#' @param local.var TBW
#' @param trace.model.var.names TBW
#' @param global.nc.attributes TBW
#'
#' @return TBW
#' @export
#'
#' @examples #TBW
downscaleTrace <- function(i, new.data.list, var.names, y1.list, y2.list, lonLim, latLim, hist.trace, data, model, local.var, trace.model.var.names, global.nc.attributes){
  
  # i <- 36
  # new.data.list <- new.trace.file.names
  # var.names <- trace.var.names
  # y1.list <- trace.years.y1
  # y2.list <- trace.years.y2
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
    
    utils::write.table(pred.df, file=paste0("../../Output/Trace21ka/", local.var, "/dat/", local.var, real.years[j], ".dat"), sep="\t")
    
    message("   ...done")
  } 
  return("Done")
} 


#' Title
#'
#' @param i TBW
#' @param new.data.list TBW
#' @param var.names TBW
#' @param y1.list TBW
#' @param y2.list TBW
#' @param lonLim TBW
#' @param latLim TBW
#' @param hist.trace TBW
#' @param data TBW
#' @param model TBW
#' @param data.bin TBW
#' @param model.bin TBW
#' @param local.var TBW
#' @param trace.model.var.names TBW
#' @param global.nc.attributes TBW
#'
#' @return TBW
#' @export
#'
#' @examples #TBW
downscaleTraceBimodel <- function (i, new.data.list, var.names, y1.list, y2.list, lonLim, latLim, hist.trace, data, model, data.bin, model.bin, local.var, trace.model.var.names, global.nc.attributes){
  
  # i <- 36
  # new.data.list <- new.trace.file.names
  # var.names <- trace.var.names
  # y1.list <- trace.years.y1
  # y2.list <- trace.years.y2
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
    
    utils::write.table(pred.df, file=paste0("../../Output/Trace21ka/", local.var, "/dat/", local.var, real.years[j], ".dat"), sep="\t")
    
    message("   ...done")
  } 
  return("Done")
} 



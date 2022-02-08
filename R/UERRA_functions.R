#' Load UERRA dataset in climate4R grid format
#'
#' @param file TBW
#' @param var TBW
#' @param lonLim TBW
#' @param latLim TBW
#' @param dictionary TBW
#' @param new.coordinates TBW
#'
#' @return TBW
#' @export
#'
#' @examples # TBW
loadUerra <- function(file, var, lonLim = uerra.lon, latLim = uerra.lat, dictionary = system.file("extdata", "UERRA_dictionary.csv", package = "TraCE21kaDSR"), new.coordinates){
  
  data <- loadeR::loadGridData(file, var = var, lonLim = lonLim, latLim = latLim, dictionary = dictionary)
  
  if(var == "pr"){
    data <- transformeR::upscaleGrid(data, times=2, aggr.fun=list(FUN=mean))
    
    data <- transformeR::interpGrid(data, new.coordinates = new.coordinates, method="bilinear")
  }
  
  data <- modifyDates(data)
  
  return(data)
}


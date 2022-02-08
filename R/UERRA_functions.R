#' Title
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
#' @examples TBW
loadUerra <- function(file, var, lonLim = uerra.lon, latLim = uerra.lat, dictionary = "../../Data/UERRA/UERRA_dictionary.dic", new.coordinates){
  
  data <- loadeR::loadGridData(file, var = var, lonLim = lonLim, latLim = latLim, dictionary = dictionary)
  
  if(var == "pr"){
    data <- transformeR::upscaleGrid(data, times=2, aggr.fun=list(FUN=mean))
    
    data <- transformeR::interpGrid(data, new.coordinates = new.coordinates, method="bilinear")
  }
  
  data <- modifyDates(data)
  
  return(data)
}


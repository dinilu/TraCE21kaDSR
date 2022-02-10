# Define parameters and argument objects
# cmip5.lon <- c(-11.25, 12.50)
# cmip5.lat <- c(27, 45)
# 
# new.trace.file.names <- mapply(FUN=function(x, y, z, w){paste0("../../Data/Trace21ka/", w, "/trace.", x, ".", y, ".cam2.h0.", w, ".", z, ".nc")}, 'trace.years.n', 'trace.years.bp', 'trace.years.nums', MoreArgs = list(w='trace.vars'), SIMPLIFY = FALSE)
# 
# cmip5.spatial.pars <- list(which.combine = 'cmip5.new.vars',
#                            v.exp = .7,
#                            rot = FALSE)
# 
# local.pars.M21 <- list(n = 1, vars = local.var)
# 
# local.pars.M24 <- list(n = 4, vars = local.var)
# 
# local.pars.M31 <- list(n = 1, vars = 'trace.model.var.names')
# 
# local.pars.M34 <- list(n = 4, vars = 'trace.model.var.names')
# 
# global.nc.attributes <- list("author" = "Diego Nieto Lugilde & Daniel Romera Romera", "institution" = "Universidad de Cordoba", "email" = "bv2nilud@uco.es")
# 
# 
# 
# 
# family.link <- "gaussian"
# 
# spatial.pars <- list(which.combine = TraCE21kaDSR:::trace.final.var.names,
#                      v.exp = .7,
#                      rot = FALSE)
# 
# hist.trace = loadHistoricalTraceGrids(traceFileNames("../Data/TraCE21ka/"))
# 
# uerra = loadUerra("../Data/UERRA/UERRA-HARMONIE/2m_temperature/latlon/1961-90_2m_temperature.nc", "tas")
# 
# require(downscaleR)
# data = prepareData(hist.trace, uerra, spatial.predictors = spatial.pars)
# 
# model = downscaleTrain(data, method = "GLM", family = family.link, predict = TRUE)
# 
# global.nc.attributes = list("author" = "Diego Nieto Lugilde & Daniel Romera Romera", "institution" = "Universidad de Cordoba", "email" = "bv2nilud@uco.es")
# 
# downscaleTrace(1, "../Output/Trace21ka/", "../Data/TraCE21ka/", hist.trace = hist.trace, mod_data = data, model = model, global.nc.attributes = global.nc.attributes)

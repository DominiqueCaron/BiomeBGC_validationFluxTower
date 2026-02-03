## Everything in this file and any files in the R directory are sourced during `simInit()`;
## all functions and objects are put into the `simList`.
## To use objects, use `sim$xxx` (they are globally available to all modules).
## Functions can be used inside any function that was sourced in this module;
## they are namespaced to the module, just like functions in R packages.
## If exact location is required, functions will be: `sim$.mods$<moduleName>$FunctionName`.
defineModule(sim, list(
  name = "BiomeBGC_validationFluxTower",
  description = "",
  keywords = "",
  authors = c(
    person("Dominique", "Caron", email = "dominique.caron@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Céline", "Boisvenue", email = "celine.boisvenue@nrcan-rncan.gc.ca", role = "ctb")
  ),
  childModules = character(0),
  version = list(BiomeBGC_validationFluxTower = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "BiomeBGC_validationFluxTower.Rmd"),
  reqdPkgs = list("SpaDES.core (>= 3.0.4)", "ggplot2"),
  parameters = bindrows(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter("resolution", "numeric", 250, NA, NA,
                    ""),
    defineParameter("targetCRS", "character", "+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs", NA, NA,
                    ""),
    defineParameter(".plots", "character", "screen", NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used - e.g., a hash of the study",
                    "area obtained using `reproducible::studyAreaName()`"),
    ## .seed is optional: `list('init' = 123)` will `set.seed(123)` for the `init` event only.
    defineParameter(".seed", "list", list(), NA, NA,
                    "Named list of seeds to use for each event (names)."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?")
  ),
  inputObjects = bindrows(
    expectsInput(objectName = "studyArea", objectClass = NA, desc = NA, sourceURL = NA),
    expectsInput(objectName = "towerMetaData", objectClass = NA, desc = NA, sourceURL = NA),
    expectsInput(objectName = "towerFluxDaily", objectClass = NA, desc = NA, sourceURL = NA),
    expectsInput(objectName = "towerFluxMonthly", objectClass = NA, desc = NA, sourceURL = NA),
    expectsInput(objectName = "towerFluxAnnual", objectClass = NA, desc = NA, sourceURL = NA),
    expectsInput(objectName = "rasterToMatch", objectClass = NA, desc = NA, sourceURL = NA),
    expectsInput(objectName = "dailyOutput", objectClass = NA, desc = NA, sourceURL = NA),
    expectsInput(objectName = "annualSummary", objectClass = NA, desc = NA, sourceURL = NA)
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "validationSummary", objectClass = NA, desc = NA)
  )
))

doEvent.BiomeBGC_validationFluxTower = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      
      # schedule future event(s)
      sim <- scheduleEvent(sim, end(sim), "BiomeBGC_validationFluxTower", "compareGPP")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "BiomeBGC_validationFluxTower", "save")
    },
    plot = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event
      
      plotFun(sim) # example of a plotting function
      # schedule future event(s)
      
      # e.g.,
      #sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "BiomeBGC_validationFluxTower", "plot")
      
      # ! ----- STOP EDITING ----- ! #
    },
    save = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event
      
      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function
      
      # schedule future event(s)
      
      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "BiomeBGC_validationFluxTower", "save")
      
      # ! ----- STOP EDITING ----- ! #
    },
    compareGPP = {
      
      # Evaluate daily predictions
      # format the flux tower data
      colToKeep <- c("TIMESTAMP", "GPP_DT_VUT_REF")
      towerData <- sim$towerFluxDaily[ , colToKeep]
      
      # switch -9999 to NA
      towerData[towerData == -9999] <- NA
      
      # format dates
      dates <- as.Date(as.character(towerData[,"TIMESTAMP"]), format = "%Y%m%d")
      # remove Feb 29th (not in BiomeBGC)
      feb29 <- format(dates, "%d") == "29" & format(dates, "%m") == "02"
      towerData <- towerData[!feb29,]
      
      years <- format(dates[!feb29], "%Y")
      nyear <- length(unique(years))
      towerData <- data.table(
        year = as.integer(years),
        day = rep(c(1:365), nyear),
        totGPPmean = towerData[,"GPP_DT_VUT_REF"]
      ) |> na.omit()
      
      dtOut <- merge(towerData, sim$dailyOutput[,.(year, timestep, day, daily_gpp)])
      # put in the same units and rename columns
      dtOut <- dtOut[, .(timestep, year, day, obsGPP = totGPPmean, predGPP = daily_gpp*1000)]
      
      # summarize the fit
      resid <- dtOut$predGPP - dtOut$obsGPP
      
      MAE <- mean(abs(resid))
      RMSE <- sqrt(mean(resid^2))
      R2 <- cor(dtOut$predGPP, dtOut$obsGPP) ^ 2
      Bias <- mean(resid)
      
      sim$validationSummary <- data.frame(
        MAE_daily = MAE,
        RMSE_daily = RMSE,
        R2_daily = R2,
        Bias_daily = Bias
      )
      
      # Evaluate monthly-level predictions
      towerData <- sim$towerFluxMonthly[ , colToKeep]
      
      # switch -9999 to NA
      towerData[towerData == -9999] <- NA
      
      # format dates
      dates <- as.Date(paste0(as.character(towerData[,"TIMESTAMP"]), "01"), format = "%Y%m%d")
      
      years <- format(dates, "%Y")
      nyear <- length(unique(years))
      towerData <- data.table(
        year = as.integer(years),
        month = rep(c(1:12), nyear),
        totGPPmean = towerData[,"GPP_DT_VUT_REF"]
      ) |> na.omit()
      
      dtOut <- merge(towerData, sim$monthlyAverages[,.(year, month, daily_gpp)])
      # put in the same units and rename columns
      dtOut <- dtOut[, .(year, month, obsGPP = totGPPmean, predGPP = daily_gpp*1000)]
      
      # summarize the fit
      resid <- dtOut$predGPP - dtOut$obsGPP
      
      sim$validationSummary$MAE_monthly <- mean(abs(resid))
      sim$validationSummary$RMSE_monthly <- sqrt(mean(resid^2))
      sim$validationSummary$R2_monthly <- cor(dtOut$predGPP, dtOut$obsGPP) ^ 2
      sim$validationSummary$Bias_monthly <- mean(resid)
      
      # Evaluate year-level predictions
      towerData <- sim$towerFluxAnnual[ , colToKeep]
      
      # switch -9999 to NA
      towerData[towerData == -9999] <- NA
      
      years <- as.numeric(towerData[,"TIMESTAMP"])
      towerData <- data.table(
        year = as.integer(years),
        totGPPmean = towerData[,"GPP_DT_VUT_REF"]
      ) |> na.omit()
      
      dtOut <- merge(towerData, sim$annualAverages[,.(year, daily_gpp)])
      # put in the same units and rename columns
      dtOut <- dtOut[, .(year, obsGPP = totGPPmean, predGPP = daily_gpp*1000*365)]
      
      # summarize the fit
      resid <- dtOut$predGPP - dtOut$obsGPP
      
      sim$validationSummary$MAE_annual <- mean(abs(resid))
      sim$validationSummary$RMSE_annual <- sqrt(mean(resid^2))
      sim$validationSummary$R2_annual <- cor(dtOut$predGPP, dtOut$obsGPP) ^ 2
      sim$validationSummary$Bias_annual <- mean(resid)
      
    },
    event2 = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event
      
      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function
      
      # schedule future event(s)
      
      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + increment, "BiomeBGC_validationFluxTower", "templateEvent")
      
      # ! ----- STOP EDITING ----- ! #
    },
    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}


### template for save events
Save <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for plot events
plotFun <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sampleData <- data.frame("TheSample" = sample(1:10, replace = TRUE))
  Plots(sampleData, fn = ggplotFn) # needs ggplot2
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event1
Event1 <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event1Test1 <- " this is test for event 1. " # for dummy unit test
  # sim$event1Test2 <- 999 # for dummy unit test
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event2
Event2 <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event2Test1 <- " this is test for event 2. " # for dummy unit test
  # sim$event2Test2 <- 777  # for dummy unit test
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  #cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")
  
  if (!suppliedElsewhere('towerMetaData', sim)) {
    
    stop("User needs to provide the metadata and ancillary data of the EC site.")
    
  }
  
  if (!suppliedElsewhere('towerMetaData', sim)) {
    
    stop("User needs to provide the flux data data of the EC site.")
    
  }
  
  if (!suppliedElsewhere('studyArea', sim)) {
    
    lat <- sim$towerMetaData[which(sim$towerMetaData$VARIABLE == "LOCATION_LAT"), "DATAVALUE"] |> as.numeric()
    lon <- sim$towerMetaData[which(sim$towerMetaData$VARIABLE == "LOCATION_LONG"), "DATAVALUE"] |> as.numeric()
    sim$studyArea <- vect(data.frame(lon = lon, lat = lat), geom=c("lon", "lat"), crs="EPSG:4326")
    sim$studyArea <- postProcessTo(sim$studyArea, projectTo =  P(sim)$targetCRS)
    
  }
  
  if (!suppliedElsewhere('rasterToMatch', sim)) {
    
    rtm <- terra::rast(buffer(sim$studyArea, P(sim)$resolution), 
                       res = c(P(sim)$resolution, P(sim)$resolution))
    terra::crs(rtm) <-  P(sim)$targetCRS
    rtm[] <- 1
    
    sim$rastertoMatch <- rtm
  }
  
  return(invisible(sim))
}

ggplotFn <- function(data, ...) {
  ggplot2::ggplot(data, ggplot2::aes(TheSample)) +
    ggplot2::geom_histogram(...)
}


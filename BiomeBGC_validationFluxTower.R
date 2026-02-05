## Everything in this file and any files in the R directory are sourced during `simInit()`;
## all functions and objects are put into the `simList`.
## To use objects, use `sim$xxx` (they are globally available to all modules).
## Functions can be used inside any function that was sourced in this module;
## they are namespaced to the module, just like functions in R packages.
## If exact location is required, functions will be: `sim$.mods$<moduleName>$FunctionName`.
defineModule(sim, list(
  name = "BiomeBGC_validationFluxTower",
  description = "Compares Biome-BGC gross primary productivity predictions to Eddy covariance flux tower estimates.",
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
  reqdPkgs = list("SpaDES.core (>= 3.0.4)", "ggplot2", "terra", "data.table", "reproducible"),
  parameters = bindrows(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter("resolution", "numeric", 250, NA, NA,
                    "Defines the resolution for the raster created for the study site."),
    defineParameter("targetCRS", "character", "+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs", NA, NA,
                    "Defines the resolution for the raster created for the study site."),
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
    expectsInput("towerCoordinates", "vector",
                 desc = paste("Must be provided by the user.",
                              "A named vector with the longitude and latitude of the EC tower.", 
                              "The names of the variables must be 'lon' and 'lat'.")
    ), 
    expectsInput("towerDailyFlux", "data.frame", 
                 desc = paste("Must be provided by the user.",
                              "The daily data from the Eddy covariance flux tower.",
                              "A csv named 'XXX_XX-XXX_FLUXNET_FLUXMET_DD_XXXX-XXXX_XXX.csv'",
                              "is usually available when downloading flux tower data.")
    ),
    expectsInput("towerMonthlyFlux", "data.frame", 
                 desc = paste("Must be provided by the user.",
                              "The monthly data from the Eddy covariance flux tower.",
                              "A csv named 'XXX_XX-XXX_FLUXNET_FLUXMET_MM_XXXX-XXXX_XXX.csv'",
                              "is usually available when downloading flux tower data.")
    ),
    expectsInput("towerAnnualFlux", "data.frame", 
                 desc =paste("Must be provided by the user.",
                             "The annual data from the Eddy covariance flux tower.",
                             "A csv named 'XXX_XX-XXX_FLUXNET_FLUXMET_YY_XXXX-XXXX_XXX.csv'",
                             "is usually available when downloading flux tower data.")
    ),
    expectsInput("dailyOutput", "data.frame", 
                 desc = paste("The daily GPP predicted by Biome-BGC.",
                              "This is an output of the module BiomeBGC_core.")
    ),
    expectsInput("monthlyAverages", "data.frame", 
                 desc = paste("The GPP predicted by Biome-BGC averaged at the monthly resolution.",
                              "This is an output of the module BiomeBGC_core.")
    ),
    expectsInput("annualAverages", "data.frame", 
                 desc = paste("The GPP predicted by Biome-BGC averaged at the monthly resolution.",
                              "This is an output of the module BiomeBGC_core.")
    )
  ),
  outputObjects = bindrows(
    createsOutput("studyArea", "SpatVector", 
                  desc = "A point vector locating the Eddy covariance flux tower."),
    createsOutput("rasterToMatch", "SpatVector", 
                  desc = "A 1 pixel raster."),
    createsOutput("validationSummary", "data.frame", 
                  desc = "A data frame with validation metrics.")
  )
))

doEvent.BiomeBGC_validationFluxTower = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      
      # schedule future event(s)
      sim <- scheduleEvent(sim, end(sim), "BiomeBGC_validationFluxTower", "compareGPP")
      # schedule plotting
      if (anyPlotting(P(sim)$.plots)) sim <- scheduleEvent(sim, end(sim), "BiomeBGC_validationFluxTower", "plot", eventPriority = 12)
    },
    plot = {
      figPath <- file.path(outputPath(sim), "BiomeBGC_validationFluxTower")
      #1.  Plot comparing daily GPP
      plot_dt <- mergeGPPdata(
        towerData = sim$towerDailyFlux,
        BiomeBGCData = sim$dailyOutput,
        timescale = "day"
      )
      
      Plots(
        plot_dt,
        fn = GPPplot,
        filename = "dailyGPP",
        types = "png",
        path = figPath,
        ggsaveArgs = list(width = 7, height = 7, units = "in", dpi = 300)
      )
      Plots(
        plot_dt,
        fn = dailyGPPtimeseries,
        filename = "dailyGPPtimeseries",
        types = "png",
        path = figPath,
        ggsaveArgs = list(width = 12, height = 7, units = "in", dpi = 300)
      )
      
      #2.  Plot comparing monthly GPP
      plot_dt <- mergeGPPdata(
        towerData = sim$towerMonthlyFlux,
        BiomeBGCData = sim$monthlyAverages,
        timescale = "month"
      )
      
      Plots(
        plot_dt,
        fn = GPPplot,
        filename = "monthlyGPP",
        types = "png",
        path = figPath,
        ggsaveArgs = list(width = 7, height = 7, units = "in", dpi = 300)
      )
      Plots(
        plot_dt,
        fn = monthlyGPPtimeseries,
        filename = "monthlyGPPtimeseries",
        types = "png",
        path = figPath,
        ggsaveArgs = list(width = 12, height = 7, units = "in", dpi = 300)
      )
      
      # 3. Annual GPP plot
      plot_dt <- mergeGPPdata(
        towerData = sim$towerAnnualFlux,
        BiomeBGCData = sim$annualAverages,
        timescale = "year",
        confInt = 95
      )
      Plots(
        plot_dt,
        fn = annualGPPplot,
        filename = "annualGPP",
        types = "png",
        path = figPath,
        ggsaveArgs = list(width = 12, height = 7, units = "in", dpi = 300)
      )
      
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
      
      #1. Evaluate daily predictions
      dayComparison <- mergeGPPdata(
        towerData = sim$towerDailyFlux,
        BiomeBGCData = sim$dailyOutput,
        timescale = "day"
      ) |> na.omit()
      
      # summarize the fit
      resid <- dayComparison$BGC_GPP - dayComparison$tower_GPP
      
      MAE <- mean(abs(resid))
      RMSE <- sqrt(mean(resid^2))
      R2 <- cor(dayComparison$BGC_GPP, dayComparison$tower_GPP) ^ 2
      Bias <- mean(resid)
      
      sim$validationSummary <- data.frame(
        MAE_daily = MAE,
        RMSE_daily = RMSE,
        R2_daily = R2,
        Bias_daily = Bias
      )
      
      #2. Evaluate monthly-level predictions
      monthComparison <- mergeGPPdata(
        towerData = sim$towerMonthlyFlux,
        BiomeBGCData = sim$monthlyAverages,
        timescale = "month"
      ) |> na.omit()
      
      # summarize the fit
      resid <- monthComparison$BGC_GPP - monthComparison$tower_GPP
      
      sim$validationSummary$MAE_monthly <- mean(abs(resid))
      sim$validationSummary$RMSE_monthly <- sqrt(mean(resid^2))
      sim$validationSummary$R2_monthly <- cor(monthComparison$BGC_GPP, monthComparison$tower_GPP) ^ 2
      sim$validationSummary$Bias_monthly <- mean(resid)
      
      #3. Evaluate year-level predictions
      yearComparison <- mergeGPPdata(
        towerData = sim$towerAnnualFlux,
        BiomeBGCData = sim$annualAverages,
        timescale = "year"
      ) |> na.omit()
      
      # summarize the fit
      resid <- yearComparison$BGC_GPP - yearComparison$tower_GPP
      
      sim$validationSummary$MAE_annual <- mean(abs(resid))
      sim$validationSummary$RMSE_annual <- sqrt(mean(resid^2))
      sim$validationSummary$R2_annual <- cor(yearComparison$BGC_GPP, yearComparison$tower_GPP) ^ 2
      sim$validationSummary$Bias_annual <- mean(resid)
      
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


dailyGPPtimeseries <- function(gpp){
  dt <- copy(gpp)
  dt[, date := as.Date(paste0(year, day), format = "%Y%j")]
  ggplot(data = dt) +
    geom_line(aes(x = date, y = BGC_GPP, col = "Biome-BGC")) +
    geom_line(aes(x = date, y = tower_GPP, col = "EC tower")) +
    scale_color_manual(name = NULL, values = c("Biome-BGC" = "darkblue", "EC tower" = "red")) +
    labs(x = "Time", y = "GPP (gC/m^2/day)") +
    scale_x_date(date_breaks = "year", date_labels = "%Y") +
    theme_classic()
}

monthlyGPPtimeseries <- function(gpp){
  dt <- copy(gpp)
  dt[, date := as.Date(paste0(year, formatC(month, digits = 1, flag = "0"), "01"), format = "%Y%m%d")]
  ggplot(data = dt) +
    geom_line(aes(x = date, y = BGC_GPP, col = "Biome-BGC")) +
    geom_line(aes(x = date, y = tower_GPP, col = "EC tower")) +
    scale_color_manual(name = NULL, values = c("Biome-BGC" = "darkblue", "EC tower" = "red")) +
    labs(x = "Time", y = "GPP (gC/m^2/day)") +
    scale_x_date(date_breaks = "year", date_labels = "%Y") +
    theme_classic()
}

GPPplot <- function(gpp){
  GPPlims <- range(c(gpp$BGC_GPP, gpp$tower_GPP), na.rm = T)
  GPPlims <- c(floor(GPPlims[1]), ceiling(GPPlims[2]))
  ggplot(data = gpp) +
    geom_point(aes(x = tower_GPP, y = BGC_GPP), alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_x_continuous(expand = c(0, 0), limits = GPPlims) + 
    scale_y_continuous(expand = c(0, 0), limits = GPPlims)  +
    labs(x = "EC tower GPP (gC/m^2/day)", y = "Biome-BGC (gC/m^2/day)") +
    theme_classic()
}

annualGPPplot <- function(gpp){
  ggplot(gpp) +
    geom_pointrange(aes(x = as.factor(year), y = tower_GPP, ymin = tower_GPP_min, ymax = tower_GPP_max, col = "EC tower"), position = position_nudge(x = -0.1)) +
    geom_point(aes(x = as.factor(year), y = BGC_GPP, col = "Biome-BGC"), position = position_nudge(x = 0.1)) +
    scale_color_manual(name = NULL, values = c("Biome-BGC" = "darkblue", "EC tower" = "red")) +
    labs(x = "Year", y = "GPP (gC/m^2/yr)") +
    theme_classic()
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
  
  if (!suppliedElsewhere('towerCoordinates', sim)) {
    
    stop("User needs to provide the metadata and ancillary data of the EC site.")
    
  }
  
  if (any(!suppliedElsewhere(c('towerDailyFlux', 'towerMonthlyFlux', 'towerAnnualFlux'), sim))) {
    
    stop("User needs to provide the flux data of the EC site.")
    
  }
  
  if (!suppliedElsewhere('studyArea', sim)) {
    
    lat <- sim$towerCoordinates["lat"] |> as.numeric()
    lon <- sim$towerCoordinates["lon"] |> as.numeric()
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


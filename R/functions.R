
# Check packages are up to date -------------------------------------------
renv::restore(lockfile = here::here('renv.lock'))

# Set options -------------------------------------------------------------

# default time zone
Sys.setenv(TZ='UTC')

# download timeout
# set to 300 sec to allow get_OISST_data() to request larger files
options(timeout = max(3000, getOption("timeout")), scipen = 1000)

# plotting options
library(ggplot2)
ggplot2::theme_set(ggplot2::theme_light())

# Create directories and key objects/values -------------------------------

# directories
if (dir.exists(here::here('input_data')) == FALSE) dir.create(here::here('input_data'))
if (dir.exists(here::here('input_data/envir')) == FALSE) dir.create(here::here('input_data/envir'))
if (dir.exists(here::here('output_data')) == FALSE) dir.create(here::here('output_data'))
if (dir.exists(here::here('output_data/plots_man')) == FALSE) dir.create(here::here('output_data/plots_man'))
if (dir.exists(here::here('output_data/envir_extract')) == FALSE) dir.create(here::here('output_data/envir_extract'))

# paths
envir <- here::here("input_data/envir")
# colonyCoord <- here::here("input_data/colony.csv")

# values
# loc.date <- format(Sys.Date(), "%d%b%Y") 
coast <- rnaturalearth::ne_countries(scale = 10, returnclass = 'sf')
envir_res <- 0.05 
min_pop <- 10000
colord <- c("Penguin", "Baccalieu", "Gull", "Great", "Corbin", "Middle Lawn", "Green", "Grand Colombier", "Scatarie", "Country", "Little White", "Long", "Bon Portage", "Kent", "Little Duck", "Great Duck")
mp <- '+proj=aea +lon_0=-60 +lat_1=30 +lat_2=60 +lat_0=45 +datum=WGS84 +units=m +no_defs'
 
# Functions - Java --------------------------------------------------------

# load java location from config.yml, ignored file specific to each collaborator
java_loc <- config::get("java_loc")

getJava <- function(out_folder = here::here(dirname(java_loc))) {

  jurl <- 'https://download.java.net/java/GA/jdk23/3c5b90190c68498b986a97f276efd28a/37/GPL/openjdk-23_windows-x64_bin.zip'

  if (dir.exists(out_folder) == FALSE) {
    dir.create(java_loc, recursive = TRUE)
    download.file(url = jurl,
                  destfile = paste0(out_folder,'/openjdk.zip'),
                  method = "auto",
                  quiet = FALSE,
                  mode="wb",
                  cacheOK = F)
    unzip(paste0(out_folder,'/openjdk.zip'), exdir = out_folder)
    unlink(paste0(out_folder,'/openjdk.zip'))
  } else print(paste('Java found in', java_loc))
}


# Settings - CMS ----------------------------------------------------------
# Set up copernicus python environment
# If using this code you should folllow help files in PelagicMagic::make_cms_envir 
# and create a config.yml file pointing to the location of your copernicus python environment. 
# A dummy version of the config file is included in the project folder. Change the file path to point to 
# your saved environment.

cm_env <- here::here(config::get("cms_env"))
PelagicMagic::make_cms_envir(my_env = cm_env)

# Functions --------------------------

# * Standardize raster resolution ----

rast_res <- function(r, res, fun = 'mean', na.rm = TRUE, focalNA, w = 3) {

  rn <- names(r)
  res_adj <- round(res/terra::res(r)[1])
  if (res_adj > 1) {
    r <- terra::aggregate(r, fact=res_adj, fun=fun, na.rm = na.rm)
  }

  if (terra::res(r)[1]) {
    suppressWarnings(
      r <- terra::disagg(r, fact=round(terra::res(r)[1]/res), method = 'bilinear')
    )
  }

  if (focalNA == TRUE) r <- terra::focal(r, fun = mean, na.policy = 'only', w = w)
  names(r) <- rn
  r
}

# * Calculate raster gradient ----
# Function to calculate new rasters that apply a focal function to an existing raster set.
# Calculations are done for every date in a list of dates, if a file with that date already
# exists in the destination folder, the date is skipped

calc_gradient_rast <- function(layers = fn, dates, destination, newlabel, w = 3, fun = 'sd', na.rm = T) {

  if (dir.exists(destination) == FALSE) dir.create(destination)
  #r <- fix_extent(fn)
  #if (terra::has.time(r) == FALSE) stop('Raster layers must have times')
  b <- terra::rast(layers)
  for (d in as.character(dates)) {
    out_file <- paste0(destination, '/',newlabel,'-', gsub("-","",d),'.nc')
    if (file.exists(out_file) == TRUE) {message(paste('Skipping', out_file, 'already exists'))
    } else {
      ##rn <- terra::rast(fn[grep(gsub('-','',d), fn)])
      rn <- b[[which(as.Date(terra::time(b)) == d)]]
      tt <- terra::time(rn)
      if (terra::ext(rn)[2] > 180) rn <- terra::rotate(rn)
      rn <- terra::focal(rn, w = w, fun = fun, na.rm = na.rm, na.policy = 'omit')
      terra::time(rn) <- tt
      message(paste('Creating file', out_file))
      terra::writeCDF(rn, out_file, varname = newlabel, missval = NA)
    }
  }
}

# Competition layers (inter- and intra-colony) ----
# Calculate the null density of birds from neighbouring colonies

calc_competition_rast <- function(colonies,
                                  colonydist,
                                  logvalues = TRUE,
                                  destination = envir,
                                  res = envir_res,
                                  focalNA = FALSE) {

  if (dir.exists(paste0(destination,'/intercolony/')) == FALSE) dir.create(paste0(destination,'/intercolony/'))
  if (dir.exists(paste0(destination,'/intracolony/')) == FALSE) dir.create(paste0(destination,'/intracolony/'))
  fn <- list.files(colonydist, full.names = T)
  cn <- sub('.nc', '', basename(fn))

  for (i in 1:length(fn)) {

    #col_idx <- grep(cn[i], colonies$colony)
    col_idx <- which(cn[i]==colonies$colony)

    pop <- colonies$population[col_idx]

    r <- terra::rast(fn[i])
    r <- rast_res(r = r, res = res, fun = 'mean', na.rm = TRUE, focalNA = focalNA, w = 3)
    names(r) <- 'colonydist'
    val <- data.frame(
      idx = 1:terra::ncell(r),
      terra::crds(r, na.rm=FALSE),
      terra::values(r, na.rm=FALSE),
      terra::cellSize(r, unit = 'km')
    )
    val$colonydist <- ifelse(is.nan(val$colonydist), NA, val$colonydist)

    temp <- val |>
      dplyr::arrange(colonydist) |>
      dplyr::mutate(
        cumarea = cumsum(area),
        comp = pop/cumarea
      )

    temp$comp <- ifelse(is.na(temp$colonydist), NA, temp$comp)

    if (i == 1) rr <- terra::rast(temp[,c('x','y','comp')], type="xyz", crs = terra::crs(r))
    if (i > 1) {
      temp <- terra::rast(temp[,c('x','y','comp')], type="xyz", crs = terra::crs(r))
      rr <- c(rr, temp)
    }
    #terra::plot(rr[[i]], main = cn[i])
  }
  names(rr) <- cn

  val <- terra::values(rr, na.rm=FALSE)

  for (i in 1:length(names(rr))) {
    if (logvalues == TRUE) terra::values(rr[[i]]) <- log(terra::values(rr[[i]]))
    #terra::plot(rr[[i]], main = names(rr)[i])
    terra::writeCDF(rr[[i]], paste0(destination,'/intracolony/',names(rr)[i],'.nc'),
                    longname = 'Intracolony competition (birds/km^2)',
                    varname = 'Intracolony',
                    unit = 'birds/km^2',
                    overwrite = T)
  }

  out <- rr
  for (i in 1:length(names(rr))) {
    if (logvalues == TRUE) terra::values(out[[i]]) <- log(rowSums(val[,-i]))
    if (logvalues == FALSE) terra::values(out[[i]]) <- (rowSums(val[,-i]))
    #terra::plot(out[[i]], main = names(out)[i])
    terra::writeCDF(out[[i]], paste0(destination,'/intercolony/',names(out)[i],'.nc'),
                    longname = 'Intercolony competition (birds/km^2)',
                    varname = 'density',
                    unit = 'birds/km^2',
                    overwrite = T)
  }

  fn <- list.files(paste0(destination,'/intracolony'), full.names = T)
  temp <- terra::rast(fn)
  if (logvalues == TRUE) temp <- exp(temp)
  r <- temp[[1]]
  terra::values(r) <- rowSums(terra::values(temp))
  if (logvalues == TRUE) r <- log(r)
  terra::writeCDF(r, paste0(destination,'/competition.nc'),
                  longname = 'Intraspecific competition (birds/km^2)',
                  varname = if (logvalues == TRUE) 'log density' else 'density',
                  unit = if (logvalues == TRUE) 'log birds/km^2' else 'birds/km^2',
                  overwrite = T)

  out
}


# Make rasters combining dynamic predictors by date -----------------------
# This function combines all dynamic predictors into a single raster by date. It ensures that
# all variables are using the same grid and NAs are filled using terra::focal

make_dyn_rasters <- function(dates, in_dir, out_dir, base_rast, vars, w = 5) {
  if (dir.exists(out_dir) == FALSE) dir.create(out_dir)

  for (d in as.character(dates)) {
    out <- terra::rast()
    out_file <- paste0(out_dir, '/', gsub('-','',d), '.nc')
    if (file.exists(out_file) == FALSE) {
      for (v in vars) {
        fn <- list.files(paste0(in_dir,'/',v), pattern = v, full.names = T, recursive = T)
        r <- terra::rast(fn)
        un <- terra::units(r)[1]
        vn <- terra::varnames(r)[1]
        r <- r[[as.Date(terra::time(r)) == as.Date(d)]]
        r <- terra::ifel(is.na(r), NA, r)
        r <- terra::focal(r, w=w, fun='mean',na.rm =T, na.policy='only')
        r <- terra::resample(r, base_rast, method = 'bilinear')
        names(r) <- v
        terra::units(r) <- un
        terra::varnames(r) <- vn
        out <- c(out,r)
      }
      terra::writeCDF(out, out_file, missval = NA, split=TRUE,
                      overwrite = T)
      print(paste('Finished:', d))
    }
  }
}

# simplified extraction function ------------------------------------------
# Simplified extraction function that uses combined standardized rasters for
# obtaining environmental predictors

get_preds <- function(data, stat_rast, col_rast, dyn_dir,
                      daynight = TRUE, verbose = TRUE,
                      convar = TRUE, static = TRUE, coldist = TRUE, dyn = TRUE) {

  if (convar) {
    # logpop
    if ('population' %in% names(data)) {
      data$logpop <- log(data$population)
    } else {stop('population is not a variable in data')}

    # region
    data$region <- ifelse(data$collat > 46, 'NL', 'SS')

    # daynight
    if (daynight == TRUE) {
      r <- suncalc::getSunlightPosition(
        data = data.frame(
          date = data$time,
          lon = sf::st_coordinates(data)[,1],
          lat = sf::st_coordinates(data)[,2]),
        keep = c('altitude')
      )
      data$daynight <- ifelse(as.numeric(r$altitude) > 0, 'day', 'night')
    }
  }

  # extract static variables
  if (static) {
    stat_rast <- terra::ifel(is.na(stat_rast), NA, stat_rast)
    data <- terra::extract(stat_rast, data, ID = FALSE, bind = T)
  }

  # extract colony distance variables
  if (coldist) {
    col_rast <- terra::ifel(is.na(col_rast), NA, col_rast)
    data$colonydist <- terra::extract(col_rast, data, ID = FALSE, layer = data$colony)$value
  }

  # extract dynamic variables
  if (dyn) {
    fn <- list.files(dyn_dir, full.names = T)
    date_idx <- sort(as.character(unique(data$date)))
    all_data <- data.frame()
    for (d in date_idx){
      dyn_rast <- terra::rast(fn[grep(gsub('-','',d), fn)])
      temp <- terra::extract(dyn_rast, data[as.character(data$date) == d,], ID = FALSE, bind = T)
      all_data <- rbind(all_data, sf::st_as_sf(temp))
      if (verbose == TRUE) print(paste('Finished:', d))
    }

    # insert 0 values for areas on the shelf which cannot have migrant micronekton
    all_data$mnkc_hmlmeso[is.na(all_data$mnkc_hmlmeso)] <- 0
    all_data$mnkc_mumeso[is.na(all_data$mnkc_mumeso)] <- 0
    data <- all_data
  }

  data

}

# # Functions - Pseudo absences ---------------------------------------------
get_buffer_pa <- function(data = dat,
                          landmask = r, # raster with water = 1, land = NA
                          nsamp = 1, # number of psuedo-absences per point
                          max_dist = max_step, # maximum distance from observed point
                          #k = 100, # number of particles to generate
                          dist_weight = F
) {

  names(landmask) <- 'lsmask'

  data$loc_idx <- 1:nrow(data)

  tt <- data |>
    dplyr::ungroup() |>
    dplyr::group_by(colony, collon, collat) |>
    dplyr::summarise(
      n = dplyr::n()
    )

  #for (i in 1:nrow(tt)) {
  out <- lapply(1:nrow(tt), FUN = function(i) {

    # sample locations from landmask raster masked to buffer distance
    mb <- sf::st_buffer(tt[i,], dist = max_dist * 1000)
    m <- terra::mask(landmask, mb)
    rp <- terra::spatSample(m, size = nsamp * tt$n[i], method = 'random', replace = T,
                            as.points = T, xy = TRUE, na.rm = T, cells = F)
    temp <- sf::st_as_sf(rp)
    temp$loc_idx <- data$loc_idx[data$colony == tt$colony[i] ]

    # add metadata about obs points
    nn <- names(data)[!(names(data) %in% names(temp))]
    temp <- dplyr::inner_join(temp, sf::st_drop_geometry(data[,c('loc_idx',nn)]), by = dplyr::join_by(loc_idx))

    temp <- temp[,names(data)]
    temp
  }
  )

  out <- do.call(rbind, out)

  out$obs <- 0
  data$obs <- 1

  rbind(data, out)

}

# Functions - Digest SDM predictions --------------------------------------

# * Bird density rasters -----
# Main function 'sdm_density_rast'
# Requires three nested functions
# Takes nBirds prediction with posterior distribution from SDM
# and makes single mean density raster with 90% CI layer for plotting

# INPUTS for 'sdm_density_rast'
# colony = ONE colony name OR "all" (for rangewide output)
# date = three character month OR list of months for multiple outputs
#        Also allows for yyyymm01 date format representing the month
# brStatus = "incubation", "chick-rearing"       [if list both, get sum]
# behaviour = "foraging", "transit"              [if list both, get sum]
# tod = "day", "night"                           [if list both, get weighted sum]
# sunLoc = c(lon, lat)
#          ONE location to determine proportion of 24hr period that is light
# source = path to predicted nBirds rasters, each with posterior distribution as layers
# destination = path to save output density rasters

## Nested Function 1
# Read and sum layers across multiple rasters (all inputs 250 layers, output 250 layers)
# Allow for single raster request (e.g. Kent, June, incubation, foraging, night)
sum_layers_across_rasts <- function(in_files) {

  if (length(in_files) == 0) stop('No files met selection criteria')

  if (length(in_files) == 1) {

    r <- terra::rast(in_files)

  } else {

    r <- terra::sprc(in_files)
    r <- terra::mosaic(r, fun = 'sum')

  }

  return(r)

}


## Nested function 2
# if(todLab == "all"), adjust nBirds for day and night length
sum_layers_across_rasts_wNightDay <- function(source = source,
                                              pattern = p,
                                              date = dd,
                                              sunLoc = sunLoc) {

  in_files <- list.files(path = source, pattern = pattern, full.names = T, recursive = T)

  if (length(in_files) == 0) stop('No files met selection criteria')

  # make day and night rasters with 250 layers representing sum of all statuses and behaviours requested
  if(length(in_files[grepl('day', in_files)]) == 1) {

    rDay <- terra::rast(in_files[grepl('day', in_files)])
    rNight <- terra::rast(in_files[grepl('night', in_files)])

  } else {

    in_files_day <- in_files[grepl('day', in_files)]
    rDay <- sum_layers_across_rasts(in_files = in_files_day)

    in_files_night <- in_files[grepl('night', in_files)]
    rNight <- sum_layers_across_rasts(in_files = in_files_night)

  }

  # Adjust nBirds in each raster for proportion day v. night on 15th of month
  pSun <- suncalc::getSunlightTimes(date = date, lon = sunLoc[1], lat = sunLoc[2], keep = c('sunrise','sunset')) |>
    dplyr::mutate(daylength = as.numeric(difftime(sunset, sunrise, units = 'hours')),
                  daylength = ifelse(is.na(daylength), 0, daylength),
                  pSun = daylength/24)
  pSun <- pSun$pSun

  rDay <- rDay*pSun
  rNight <- rNight*(1-pSun)

  # Sum adjusted rasters
  r <- rDay + rNight

  return(r)

}

## Nested function 3
# Convert raster to density, get median and 90% CI
birds_to_densityMedCI <- function(rast,
                                  units = 'km',
                                  probs = c(0.5, 0.05, 0.95)) {

  # Convert to density using actual cell size (default birds/km2)
  d <- rast/terra::cellSize(rast, unit = units)

  # Default gets median density and 90% CI across the layers
  d <- terra::quantile(d, probs=probs, na.rm = T)

}


## Main function
sdm_density_rast <- function(colony = NULL,
                             date,
                             brStatus = c('incubation','chick-rearing'),
                             behaviour = c('foraging','transit'),
                             tod = c('night','day'),
                             sunLoc = c(-52.7988, 48.1165),    # defaults to Baccalieu
                             source,
                             destination) {

  if (dir.exists(destination) == FALSE) dir.create(destination)
  if (length(colony) > 1) stop('colony must be length 1')

  # Deal with variables vs file naming labels
  if (colony == "all") colLab <- "AllColonies" else colLab <- colony
  if (length(brStatus)>1) brStatusLab <- "combined" else brStatusLab <- brStatus
  if (length(behaviour)>1) behaviourLab <- "both" else behaviourLab <- behaviour
  if (length(tod)>1) todLab <- "daily" else todLab <- tod

  # After making colLab, load list of all colonies to 'colony' if requested
  if (colony == "all") colony <- list.dirs(source, full.names = F, recursive = F)

  # Find all combinations of variables to pattern match when reading rasters
  pGrid <- expand.grid(colony = colony, status = brStatus, behaviour = behaviour, tod = tod)

  # Loop to make separate output for each date/month
  for (i in 1:length(date)){

    ## This should work for single colony or all colonies (add EVERYTHING together unless night-day involved) ##

    # make a consistent date object when you begin looping on date [better to convert all input months before loop?]
    if (date[i] %in% month.abb) dd <- as.Date(paste0('2025',date[i],15), format = '%Y%b%d') else {
      dd <- as.Date(date[i], format = '%Y-%m-%d')
    }

    # Make out_file name, check if product already exists
    out_file <- paste0(destination,"/", colLab, "_", date[i], "_", brStatusLab, "_", behaviourLab, "_", todLab, "_density.tif")

    if (file.exists(out_file) == TRUE) {message(paste('Skipping', out_file, 'already exists'))
    } else {

      p <- paste(paste0(pGrid$colony, "_", date[i], "_", pGrid$status, "_", pGrid$behaviour, "_", pGrid$tod), collapse = '|')

      # If todLab = daily, nBirds need to be adjusted by proportion light and dark before summing rasters
      # If not, nBirds rasters can simply be summed across colonies, breeding status(es), behaviour(s)
      if(todLab == "daily"){

        r <- sum_layers_across_rasts_wNightDay(source = source,
                                               pattern = p,
                                               date = dd,
                                               sunLoc = sunLoc)

      } else {

        in_files <- list.files(source, pattern = p, full.names = T, recursive = T)

        r <- sum_layers_across_rasts(in_files = in_files)

      }

      # Convert resulting nBirds raster to density, get quantiles of interest
      d <- birds_to_densityMedCI(rast = r, units = 'km', probs = c(0.5, 0.05, 0.95))

      # Save
      terra::writeRaster(d,
                         out_file,
                         overwrite=TRUE,
                         gdal=c("of=COG",
                                "COPY_SRC_OVERVIEWS=YES",
                                "COMPRESS=DEFLATE"),
                         filetype='COG'
      )

    }

  }

}

### Function for extracting densities from summary rasters.

extract_dens <- function(source = od, extract_areas = pa, modlab = 'Breeding') {
  fn <- list.files(source, full.names = T)

  out <- lapply(fn, function(f) {
    r <- terra::rast(f)
    ds <- terra::extract(r, extract_areas, weights = T, exact = T, fun = 'mean', na.rm = T)
    data.frame(
      ProdArea = extract_areas$ProdArea,
      area_km2 = extract_areas$area_km2,
      Province = extract_areas$Province,
      license_type = extract_areas$license_type,
      colony = strsplit(basename(f), split = '_')[[1]][1],
      month = strsplit(sub('.tif','',basename(f)), split = '_')[[1]][2],
      med = ds$q0.5,
      ll = ds$q0.05,
      ul = ds$q0.95)
  })
  out <- do.call(rbind, out)

  out$colony <- factor(out$colony, levels = colord)
  out$month <- factor(out$month, levels = month.abb, labels = month.name)
  out$model <- modlab
  out
}


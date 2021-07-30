#' Convert every grid-cell of a gridded pracipitation product into a station for SWAT+
#'
#' @param product Precipitation product as 'SpatRast' that includes all the days of the study period as layers.
#' @param dem Digital elevation model to be used to extract the elevation of the stations.
#' @param output Folder to write the station files (data and metadata).
#' @param hru.shape Spatial object of type 'SpatVector' that contains the HRUs as obtained from QSWAT+
#' @param start Character object that contains the starting date of the model in the format '%Y%m%d' (e.g., '19810101').
#'
#' @return
#' @export
#'
#' @examples
prec2stations <- function(product, dem, output, hru.shape, start){

  if(terra::crs(product) != terra::crs(hru.shape))
    hru.shape.pj <- terra::project(hru.shape, product)

  # Crop the product to the catchment extent
  product <- terra::crop(product, hru.shape.pj, snap = "out")

  # Raster to points
  product.points <- terra::as.points(product[[1]])

  # Defining variable name
  var <- "pcp"

  # Extracting elevation
  if(terra::crs(dem) != terra::crs(product.points))
    dem <- terra::project(dem, terra::crs(product.points))

  elev <- terra::extract(dem, product.points)[,2]

  # Defining which points don't have elevation
  tmp  <- which(is.nan(elev))
  
  # Neglecting those points in the 'product.points'
  if(length(tmp) > 0){
    elev <- elev[-tmp]
    product.points <- product.points[-tmp,]
  }


  # Setting the latitude and longitude based on the virtual stations
  lat   <- terra::geom(product.points)[,4]
  lon   <- terra::geom(product.points)[,3]

  # Extract timeseries
  product.ts <-terra::extract(product, product.points)

  # Genetating the 'pcp.txt' file
  id    <- 1:nrow(product.ts)
  nam   <- paste0("basin", var, id)

  metadata <- data.frame(ID = id, NAME = nam, LAT = lat,
                         LONG = lon, ELEVATION = elev)

  # Writting metadata file
  write.table(metadata, paste0(output, "/", var, ".txt"),
              row.names = FALSE, sep = ",", quote = FALSE)

  ## Writting the single files
  for(i in 1:length(id)){

    ts <- product.ts[i,]

    if(names(ts)[1] == "ID") {
      ts <- round(as.numeric(ts),2)[2:length(ts)]
    } else {
      ts <- round(as.numeric(ts),2)
    }

    ts <- data.frame(ts)
    names(ts) <- start

    write.table(ts, paste0(output, "/", nam[i], ".txt"),
                row.names = FALSE, sep = ",", quote = FALSE)

  }

}


#' Convert every grid-cell of a gridded maximum and minimum temperature product into a station for SWAT+
#'
#' @param tmax Maximum temperatures product as 'SpatRast' that includes all the days of the study period as layers.
#' @param tmin Minimum temperatures product as 'SpatRast' that includes all the days of the study period as layers.
#' @param dem Digital elevation model to be used to extract the elevation of the stations.
#' @param output Folder to write the station files (data and metadata).
#' @param hru.shape Spatial object of type 'SpatVector' that contains the HRUs as obtained from QSWAT+
#' @param start Character object that contains the starting date of the model in the format '%Y%m%d' e.g., '19810101'.
#'
#' @return
#' @export
#'
#' @examples
temp2stations <- function(tmax, tmin, dem, output, hru.shape, start){

  if(terra::crs(tmax) != terra::crs(hru.shape))
    hru.shape.pj <- terra::project(hru.shape, tmax)

  # Crop the product to the catchment extent
  tmax <- terra::crop(tmax, hru.shape.pj, snap = "out")
  tmin <- terra::crop(tmin, hru.shape.pj, snap = "out")

  # Raster to points
  product.points <- terra::as.points(tmax[[1]])

  # Defining variable name
  var <- "tmp"

  # Extracting elevation
  if(terra::crs(dem) != terra::crs(product.points))
    dem <- terra::project(dem, terra::crs(product.points))

  elev <- terra::extract(dem, product.points)[,2]

  # Defining which points don't have elevation
  tmp  <- which(is.nan(elev))
  
  # Neglecting those points in the 'product.points'
  if(length(tmp) > 0){
    elev <- elev[-tmp]
    product.points <- product.points[-tmp,]
  }

  # Setting the latitude and longitude based on the virtual stations
  lat   <- terra::geom(product.points)[,4]
  lon   <- terra::geom(product.points)[,3]

  # Extract timeseries
  tmax.ts <-terra::extract(tmax, product.points)
  tmin.ts <-terra::extract(tmin, product.points)

  # Genetating the 'pcp.txt' file
  id    <- 1:nrow(tmax.ts)
  nam   <- paste0("basin", var, id)

  metadata <- data.frame(ID = id, NAME = nam, LAT = lat,
                         LONG = lon, ELEVATION = elev)

  # Writting metadata file
  write.table(metadata, paste0(output, "/", var, ".txt"),
              row.names = FALSE, sep = ",", quote = FALSE)

  ## Writting the single files
  for(i in 1:length(id)){

    ts.max <- tmax.ts[i,]
    ts.min <- tmin.ts[i,]

    if(names(ts.max)[1] == "ID") {
      ts.max <- round(as.numeric(ts.max),2)[2:length(ts.max)]
    } else {
      ts.max <- round(as.numeric(ts.max),2)
    }

    if(names(ts.min)[1] == "ID") {
      ts.min <- round(as.numeric(ts.min),2)[2:length(ts.min)]
    } else {
      ts.min <- round(as.numeric(ts.min),2)
    }

    ts <- data.frame(ts.max, ts.min)

    write.table(start, paste0(output, "/", nam[i], ".txt"),
                row.names = FALSE, sep = ",", quote = FALSE, col.names = FALSE)

    write.table(ts, paste0(output, "/", nam[i], ".txt"),
                append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

  }

}


#' Generate a data.frame with the land cover of the HRUs for every land cover given as an input
#'
#' @author Oscar M. Baez-Villanueva
#' @param hru.shape Spatial object of type 'SpatVector' that contains the HRUs as obtained from QSWAT+.
#' @param lc Raster object of type 'SpatRaster' that contains the land cover for the years to be analysed.
#' @param lookup.table A three column data.frame. The first column stores the land use values from the 'lc' object. The
#'    second column stores the SWAT_CODE as in the lookup table that is given to SWAT+. Finally, the third column represents
#'    the land use (lum) identifiers from SWAT+ that can be obtained from the 'landuse.lum' file.
#' @param fun Function to be used to calculate the land cover over each HRU. The 'modal' function from the terra package
#'    is reccomended.
#' @param fact.diss Factor of disaggregation of the land cover maps. Sometimes, when the HRUs are very small,
#'    this function can return NAs. Therefore, the LU maps are disaggregated using the nearest neighbour. A value of 1 means
#'    no disaggregation. See terra::disaggregate for more details.
#'
#' @return This function returns a object type 'data.frame' for all HRUs with their respective land use values obtaied through the function 'fun'
#' @export
#'
#' @examples
extract_changes <- function(hru.shape, lc, lookup.table, fun = modal, fact.diss = NULL) {   # change the HRU area to KM2

  # Cropping the land cover stack to the catchment extent
  shape.proj <- terra::project(hru.shape, terra::crs(lc))

  # Cropping the land cover stack to the catchment extent
  lc <- terra::crop(lc, shape.proj, snap = "out")

  # Projecting the lc object
  lc <- terra::project(lc, terra::crs(hru.shape), method = "near")

  # Dissagregating the land cover to avoid NAs in small HRUs
  if(!is.null(fact.diss))
    lc <- terra::disaggregate(lc, fact.diss, method = "near")
  
  # Masking the lc object
  lc <- terra::mask(lc, hru.shape)

  # Extract the values of the hru's according to the defined function
  vals <- terra::extract(lc, hru.shape, fun = fun)
  vals <- data.frame(apply(vals, 2, round, digits = 0))

  # Setting IDs HRUs
  hrus.df <- terra::values(hru.shape)
  ID      <- hrus.df$HRUS

  vals$ID <- as.numeric(ID)
  vals    <- vals[order(vals$ID),]

  # Getting the SWAT+ codes from the vals data.frame
  swat.codes <- vals

  return.codes <- function(x, df, n){
    r <- 0

    for(i in 1:length(x)){
      r[i] <- df[x[i],n]
    }
    return(r)
  }


  for(i in 2:ncol(swat.codes)){
    swat.codes[,i] <- return.codes(swat.codes[,i], df = lookup.table, n = 2)
  }

  # Getting the SWAT+ land use management ids from the vals data.frame
  lum.ids <- vals

  for(i in 2:ncol(lum.ids)){
    lum.ids[,i] <- return.codes(lum.ids[,i], df = lookup.table, n = 3)
  }

  # Generating the list that will be returned
  result <- list(LU_values = vals,
                 SWAT_codes = swat.codes,
                 LUM_ids = lum.ids)

  return(result)
}

#' Plot the changes of the given land use layers as pie charts
#'
#' @author Oscar M. Baez-Villanueva
#' @param df.changes Object of type 'data.frame' obtained from
#'    applying the function 'extract_changes'.
#' @param ncolumns Number of columns of the resulting figure.
#' @param df.codes Data frame of two columns. The first column 'ID'
#'    must contain the codes of the land uses, while the second column
#'    'names' contain the respective land cover names.
#'
#' @return This function returns a 3D pie plot for every land cover layer
#'    given to the 'extract_changes' function.
#' @export
#'
#' @examples
plot_changes <- function(df.changes, ncolumns, df.codes){

  n     <- ncol(df.changes) - 1
  nrows <- ceiling(n/ncolumns)

  par(mfrow = c(nrows, ncolumns))

  for(i in 1:n){

    slices <- df.changes[,i+1]
    slices <- data.frame(ftable(slices))
    labels <- df.codes[which(df.codes[,1] %in% slices$slices),2]
    slices <- slices$Freq
    main   <- colnames(df.changes)[i+1]

    plotrix::pie3D(slices, labels = labels, explode = 0.1, main = main, labelcex = 0.55, height = 0.1)
  }

}



#' Format ground-based precipitation data according to `SWAT+` requirements
#'
#' @param data zoo object with the precipitation time series of all stations. The dates should have the format \%Y-\%m-\%d.
#' @param metadata data.frame object containing the metadata of the precipitation stations. This data.frame must have five fields in the following specific order:
#     ID: number of the stations, starting from 1 to the total number of stations; NAME: the code of each station; LAT: the latitude coordinates in WGS84;
#'   LONG: the longitude coordinates in WGS84; and ELEVATION: the elevation of the station in metres.
#' @param output character object with the path of the directory where the data will be stored.
#'
#' @return This function formats the ground-based precipitation stations data according to SWAT+ requirements.
#' @export
#'
#' @examples
formatPstations <- function(data, metadata, output){
  
  names_metadata <- c("ID", "NAME", "LAT", "LONG", "ELEVATION")
  
  if(!(identical(names(metadata), names_metadata)))
    stop("The metadata must have the following columns: ", paste(names_metadata, collapse = ", "))
  
  if(!zoo::is.zoo(data))
    stop("The data object must be in zoo format with the dates as '%Y-%m-%d'")

   # Writting metadata file
  write.table(metadata, paste0(output, "/pcp.txt"),
              row.names = FALSE, sep = ",", quote = FALSE)
  
  start <- zoo::index(data)[1]
  start <- gsub('[-]', '', start)
  
  ## Writting the single files
  for(i in 1:ncol(data)){
    
    code <- names(data)[i]
    ts   <- as.numeric(data[,i])
    ts   <- data.frame(ts)
    
    names(ts) <- start
    
    write.table(ts, paste0(output, "/", code, ".txt"),
                row.names = FALSE, sep = ",", quote = FALSE)
    
  }
  
}

#' Format ground-based temperature data according to `SWAT+` requirements
#'
#' @param dataTmax zoo object with the maximum temperature time series of all stations. The dates should have the format \%Y-\%m-\%d.
#' @param dataTmin zoo object with the minimum temperature time series of all stations. The dates should have the format \%Y-\%m-\%d.
#' @param metadata data.frame object containing the metadata of the precipitation stations. This data.frame must have five fields in the following specific order:
#     ID: number of the stations, starting from 1 to the total number of stations; NAME: the code of each station; LAT: the latitude coordinates in WGS84;
#'   LONG: the longitude coordinates in WGS84; and ELEVATION: the elevation of the station in metres.
#' @param output character object with the path of the directory where the data will be stored.
#'
#' @return This function formats the ground-based precipitation stations data according to SWAT+ requirements.
#' @export
#'
#' @examples
formatTstations <- function(dataTmax, dataTmin, metadata, output){
  
  names_metadata <- c("ID", "NAME", "LAT", "LONG", "ELEVATION")
  
  if(!(identical(names(metadata), names_metadata)))
    stop("The metadata must have the following columns: ", paste(names_metadata, collapse = ", "))
  
  if(!(nrow(dataTmax) == nrow(dataTmin)))
    stop("dataTmax and dataTmin have different number of rows")
  
  if(!(ncol(dataTmax) == ncol(dataTmin)))
    stop("dataTmax and dataTmin have different number of columns")
  
  if(!(identical(names(dataTmax), names(dataTmin))))
    stop("The order of the stations is not the same in 'dataTmax' and 'dataTmin'")
  
  if(!zoo::is.zoo(dataTmax) & !zoo::is.zoo(dataTmin))
    stop("The dataTmax and dataTmin objects must be in zoo format with the dates as '%Y-%m-%d'")
  
  # Writting metadata file
  write.table(metadata, paste0(output, "/tmp.txt"),
              row.names = FALSE, sep = ",", quote = FALSE)
  
  start <- zoo::index(dataTmax)[1]
  start <- gsub('[-]', '', start)
  
  ## Writting the single files
  for(i in 1:ncol(dataTmax)){
    
    code <- names(dataTmax)[i]

    ts.max <- dataTmax[,i]
    ts.min <- dataTmin[,i]
    
    ts <- data.frame(ts.max, ts.min)
    
    write.table(start, paste0(output, "/", code, ".txt"),
                row.names = FALSE, sep = ",", quote = FALSE, col.names = FALSE)
    
    write.table(ts, paste0(output, "/", code, ".txt"),
                append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
   }
  
}


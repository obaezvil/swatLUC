


#' Generate a new land use change scenario in SWAT+
#'
#' @param scen.name Object type 'character' with the name of the new scenario to be created.
#' @param textInOut.path Path to the default 'TxtInOut' folder. The new scenario will be created from this file.
#' @param hru.shape Spatial object of type 'SpatVector' that contains the HRUs that were created in the `QSWAT+` project.
#' @param lc Raster object of type 'SpatRaster' that contains the land cover for the years to be analysed.
#' @param lookup.table A three column 'data.frame' object. The first column stores the land use values from the 'lc' 
#'    object. The second column stores the SWAT_CODE as in the lookup table that is given to `SWAT+`. Finally, the 
#'    third column represents the land use (lum) identifiers from `SWAT+` as stored in the 'landuse.lum' file.
#' @param skip.copy Logical object. Set tu true only if the files from the default 'TxtInOut' folder have been 
#'    already copied to the directory of the new scenario.
#' @param change.yr Object type 'character' containing the year(s) when the land use changes take place (in chronological order).
#' @param jday Numerical value that indicates the julian day that the land use change(s) will take place (default jday = 1).
#' @param fact.diss Factor of disaggregation of the land cover maps. Sometimes, when the HRUs are very small,
#'    this function can return NAs. Therefore, the LU maps are disaggregated using the nearest neighbour. A value of 1 means
#'    no disaggregation. See terra::disaggregate for more details.
#' @param scen.path Path to the directory where the new scenatio will be created. If not given, the new scenario will be
#'    created in the same folder of the default scenario.
#' @param verbose Object type 'logical'. Should the progress messages be printed?
#' @param exclude.uses Exclude land use classes from the scenario. Currently, the change of urban and water classes cause errors in `SWAT+`; and therefore, should be avoided.
#'
#' @return
#' @export
#'
#' @examples
LUCscen <- function(scen.name,
                    textInOut.path,
                    hru.shape,
                    lc,
                    lookup.table,
                    change.yr,
                    skip.copy = FALSE,
                    jday = 1,
                    fact.diss = NULL,
                    scen.path = NULL,
                    verbose = TRUE,
                    exclude.uses = NULL){

  #-------------Checks--------------------------#
  if(!is.numeric(change.yr))
    stop("The object 'change.yr' is not numeric!")

  n <- terra::nlyr(lc)

  if(length(change.yr) != n - 1)
    stop("The object 'change.yr' has different dimensions than the expected changes!")


  # Checking where to save new scenario
  if(!is.null(scen.path)){
    scenarios <- scen.path
  } else {
    scenarios <- dirname(dirname(textInOut.path))
  }

  # Creating new folder to store the new
  new.folder <- file.path(scenarios, scen.name)

  if(!file.exists(new.folder))
    dir.create(new.folder)

    # Copying the defauld scenario inside new scenario
  if(!skip.copy){

    if(verbose){
      message("#------ Creating new scenario files in: ------------# \n")
      message(new.folder, "\n")
      file.copy(textInOut.path, new.folder, recursive = TRUE)
    }

  } else {

      if(verbose)
        message("#------ The file was not copied! ------------# \n")

  }

  if(verbose)
    message("#------ Starting the analysis of LU changes --------# \n")

  # Setting the path to the new path
  textInOut.new <- file.path(new.folder, basename(textInOut.path))

  # get decision table file
  scen.lu.file <- file.path(textInOut.new, "scen_lu.dtl")

  #---------------- Setting the conditions ------------------------#


  var       <- c(rep("year_cal", n - 1), "jday")
  obj       <- rep("null", n)
  obj_num   <- rep(0, n)
  lim_var   <- obj
  lim_op    <- rep ("-", n)
  lim_const <- c(sprintf('%.5f', change.yr), sprintf('%.5f', jday))

  # Setting number of conditions
  n.cond <- terra::nlyr(lc)

  # Creating the conditions section
  conditions <- data.frame(
    var = var,
    obj = obj,
    obj_num = obj_num,
    lim_var = lim_var,
    lim_op = lim_op,
    lim_const = lim_const
  )

  #---------------- Setting alternatives ---------------------------#
  n.alts    <- n - 1
  name.alts <- paste0("alt", 1:n.alts)

  # Creating an identity matrix to set the alternatives of LU change
  alts <- diag(n.alts)
  alts[alts == 1]   <- "="
  alts[alts == "0"] <- "-"

  # Transform the matrix to data frame and add the jday row
  alts <- rbind(alts, rep("=", ncol(alts)))
  alts <- data.frame(alts)
  names(alts) <- name.alts

  # Creating the alternatives section
  alternatives <- alts

  #---------------- Setting activities -----------------------------#
  # Apply function to analyse the changes according to HRU and land covers
  changes.df <- extract_changes(hru.shape, lc, fun = terra::modal, fact.diss = fact.diss, lookup.table = lookup.table)

  # Generating list to store the analysis per change
  activities.list <- list()

  # Setting vectod for outcomes
  outc.general <- rep("n", n.alts)

  # Creating new folder to store the shapefiles
  shape.folder <- file.path(new.folder, "Shapefiles_HRU_Changes")

  if(!file.exists(shape.folder))
    dir.create(shape.folder)

  # Read hru-data.hru
  hrus.file <- file.path(textInOut.new, "hru-data.hru")
  hrus      <- read.table(hrus.file, skip = 1, header = TRUE, sep = "")

  # keeping only character values of the changes
  changes.lum <- changes.df$LUM_ids[,2:ncol(changes.df$LUM_ids)]

  # Get changes
  for(i in 1:n.alts){
    
    # Storing id's of the HRUs
    ids <- changes.df$LU_values[,1]
    
    # Substracting the columns of the change 'i'
    subset.df        <- data.frame(changes.lum[,i], changes.lum[,i+1])
    names(subset.df) <- c("Before", "After")

    before <- subset.df[,1]
    after  <- subset.df[,2]

    # Subsetting HRUs according to NAs
    exclude <- which(is.na(before))
    exclude <- c(exclude, which(is.na(after)))
    exclude <- unique(exclude)

    if(length(exclude) > 0){
      before <- before[-exclude]
      after  <- after[-exclude]
      ids    <- ids[-exclude]
      
      if(verbose)
        Warning("There are NA values in the analysis, try increasing the 'fact.diss' parameter! \n")
    }
  
    # Subsetting HRUs according to exclude cases
    exclude <- which(before %in% exclude.uses)
    exclude <- c(exclude, which(after %in% exclude.uses))
    exclude <- unique(exclude)

    if(length(exclude) > 0){
      before <- before[-exclude]
      after  <- after[-exclude]
      ids    <- ids[-exclude]
    }

    # Cross-checking the land uses with the HRUs file
    lum.hrus <- hrus$lu_mgt[which(hrus$id %in% ids)]
    exclude  <- which(lum.hrus %in% exclude.uses)

    if(length(exclude) > 0){
      before <- before[-exclude]
      after  <- after[-exclude]
      ids    <- ids[-exclude]
    }

    # Warning of new clases!!
    new.elements <- after[which(!(after %in% before))]
    new.elements <- unique(new.elements)
    if(length(new.elements) > 0)
      warning("Warning! The land cover(s): '", new.elements, "' is(are) not contained in the initial conditions! \n",
            "Please make sure that it(they) is(are) defined in the 'landuse.lum', 'cntable.lum', 'plant.ini', and 'ovn_table.lum' files! \n",
            "In the case that the new land cover(s) is(are) included in 'landuse.lum' please omit this warning :)")

    # Subsetting HRUs that changed
    pos <- which(before != after)

    act_type <- rep("lu_change", length(pos))
    obj      <- rep("hru", length(pos))
    obj_num  <- ids[pos]
    name     <- after[pos]
    option   <- rep("null", length(pos))
    const    <- rep(sprintf('%.5f', 0), length(pos))
    const2   <- const
    fp       <- name

    outcome    <- outc.general
    outcome[i] <- "y"
    outcome    <- noquote(paste(outcome, collapse = "   "))
    outcome    <- rep(outcome, length(pos))

    # Creating the activities section
    activities.list[[i]] <- data.frame(
      act_type = act_type,
      obj = obj,
      obj_num = obj_num,
      name = name,
      option = option,
      const = const,
      const2 = const2,
      fp = fp,
      outcome = outcome
    )

    # Create shapefile of the changes

    # Getting numerical values of the HRUs of the shapefile
    hrus.nbr <- as.numeric(hru.shape$HRUS)

    # Matching them to the HRUs that changed
    subset   <- which(hrus.nbr %in% obj_num)

    # Subsetting shapefile
    hru.change <- hru.shape[subset,]

    # Getting attribute table of the HRUs shapefile and the specific HRUs
    values.hru <- terra::values(hru.change)
    hrus.shape <- as.numeric(values.hru$HRUS)

    # Creating and ordering the lum values for before and after
    before.hru <- data.frame(ids, before)
    after.hru  <- data.frame(ids, after)
    before.hru <- before.hru[match(hrus.shape, before.hru$ids),]
    after.hru  <- after.hru[match(hrus.shape, after.hru$ids),]

    # Assigning values to the shapefile
    values.hru$Before  <- before.hru$before
    values.hru$After   <- after.hru$after
    terra::values(hru.change) <- values.hru

    name.shape <- paste0(shape.folder, "/HRU_change_", change.yr[i])
    if(!file.exists(name.shape))
      terra::writeVector(hru.change, name.shape, overwrite = TRUE)

  }

  activities <- activities.list
  names(activities) <- paste0("LUchange_", change.yr)
  activities        <-  do.call("rbind", activities)

  n.acts <- nrow(activities)

  #---------------- Write Desision Table ---------------------------#

  # Read the decision tables file
  scen.lu    <- readLines(scen.lu.file)
  scen.lu[2] <- as.numeric(scen.lu[2]) +  1

  # Get the maximum numbers of lines
  max.lines <- length(scen.lu)

  # If the last line of the file is not blank then create it
  if(!scen.lu[max.lines] == ""){
    max.lines <-max.lines + 1
    scen.lu[max.lines] <- ""
  }

  # Write the firs line of the table
  scen.lu[max.lines + 1] <- "name                     conds      alts      acts  "

  # Creating the second line of the file
  name  <- stringi::stri_pad(scen.name, width = 29, side = "right")
  conds <- stringi::stri_pad(n.cond, width = 10, side = "right")
  alts  <- stringi::stri_pad(n.alts, width = 8, side = "right")
  acts  <- stringi::stri_pad(n.acts, width = 4, side = "right")
  scen.lu[max.lines + 2] <- paste0(name, conds, alts, acts)

  # resseting the max.lines
  nex.line <- length(scen.lu) + 1

  # Add title conditions and alternatives
  cond.names    <- "var                        obj   obj_num           lim_var            lim_op     lim_const"
  alt.names     <- unlist(lapply(names(alternatives), stringi::stri_pad, width = 9, side = "left"))
  alt.names     <- paste(paste0(alt.names), collapse = "")

  condalt.names <- paste(cond.names, alt.names)

  scen.lu[nex.line] <- condalt.names
  nex.line          <- nex.line + 1

  # Adding conditions and alternatives
  cond      <- conditions
  var       <- stringi::stri_pad(cond$var, width = 26, side = "right")
  obj       <- stringi::stri_pad(cond$obj, width = 13, side = "right")
  obj_num   <- stringi::stri_pad(cond$obj_num, width = 15, side = "right")
  lim_var   <- stringi::stri_pad(cond$lim_var, width = 21, side = "right")
  lim_op    <- stringi::stri_pad(cond$lim_op, width = 5, side = "right")
  lim_const <- stringi::stri_pad(cond$lim_const, width = 19, side = "right")

  alt        <- data.frame(apply(alternatives, 2, stringi::stri_pad,
                                 width = 11, side = "right"))
  names(alt) <- names(alternatives)

  condalt <- cbind(data.frame(
    var = var,
    obj = obj,
    obj_num = obj_num,
    lim_var = lim_var,
    lim_op = lim_op,
    lim_const = lim_const),
    alt
  )

  for(i in 1:nrow(condalt)){

    new.line <- length(scen.lu) + 1
    scen.lu[new.line] <- paste(condalt[i,], collapse = "")

  }

  # Creating header for activities
  new.line <- new.line + 1
  scen.lu[new.line] <- "act_typ                    obj   obj_num              name            option         const        const2                fp  outcome           "

  # Adding conditions and activities
  act      <- activities
  act_typ  <- stringi::stri_pad(act$act_type, width = 27, side = "right")
  obj      <- stringi::stri_pad(act$obj, width = 9, side = "right")
  obj_num  <- stringi::stri_pad(act$obj_num, width = 14, side = "right")
  name     <- stringi::stri_pad(act$name, width = 22, side = "right")
  option   <- stringi::stri_pad(act$option, width = 11, side = "right")
  const    <- stringi::stri_pad(act$const, width = 14, side = "right")
  const2   <- stringi::stri_pad(act$const2, width = 17, side = "right")
  fp       <- stringi::stri_pad(act$fp, width = 10, side = "right")
  out      <- act[,9:ncol(act)]

  activ <- data.frame(
    act_typ = act_typ,
    obj = obj,
    obj_num = obj_num,
    name = name,
    option = option,
    const = const,
    const2 = const2,
    fp = fp,
    out
  )

  for(i in 1:nrow(activ)){

    new.line <- length(scen.lu) + 1
    scen.lu[new.line] <- paste(activ[i,], collapse = "")

  }

  # Adding last line and saving the table
  scen.lu[new.line + 1] <- ""

  writeLines(scen.lu, scen.lu.file)


  #---------------- Write conditional.upd ---------------------------#

  cond.upd    <- paste("Conditional.upd: written by swatLUC:", Sys.time(), collapse = "")
  cond.upd[2] <- paste("1  !Ccond must match a decision table in scen_lu.dtl", collapse = "")
  cond.upd[3] <- paste("type        name            cond", collapse = "")

  type <- stringi::stri_pad("lu_change", width = 11, side = "right")
  name <- stringi::stri_pad(scen.name, width = 15, side = "right")
  cond <- scen.name

  cond.upd[4] <- paste(type, name, cond, collapse = "")
  cond.upd[5] <- ""

  cond.upd.file <- file.path(dirname(scen.lu.file), "conditional.upd")

  writeLines(cond.upd, cond.upd.file)

}

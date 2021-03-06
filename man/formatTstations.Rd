% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_stations.R
\name{formatTstations}
\alias{formatTstations}
\title{Format ground-based temperature data according to `SWAT+` requirements}
\usage{
formatTstations(dataTmax, dataTmin, metadata, output)
}
\arguments{
\item{dataTmax}{zoo object with the maximum temperature time series of all stations. The dates should have the format \%Y-\%m-\%d.}

\item{dataTmin}{zoo object with the minimum temperature time series of all stations. The dates should have the format \%Y-\%m-\%d.}

\item{metadata}{data.frame object containing the metadata of the precipitation stations. This data.frame must have five fields in the following specific order:
LONG: the longitude coordinates in WGS84; and ELEVATION: the elevation of the station in metres.}

\item{output}{character object with the path of the directory where the data will be stored.}
}
\value{
This function formats the ground-based precipitation stations data according to SWAT+ requirements.
}
\description{
Format ground-based temperature data according to `SWAT+` requirements
}

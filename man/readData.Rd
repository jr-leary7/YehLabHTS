% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readData.R
\name{readData}
\alias{readData}
\title{Read .xlsx HTS data}
\usage{
readData(parent.dir, file.name, sheet.name = 1, col.names = TRUE)
}
\arguments{
\item{parent.dir}{The directory within which the .xlsx file is located.}

\item{file.name}{The .xlsx file you wish to load.}

\item{sheet.name}{Which sheet to load.}

\item{col.names}{Should column names be loaded? Defaults to TRUE.}
}
\description{
This function reads in a \code{.xlsx} file, and outputs the name of the file being read before returning it. Useful for keeping track of data in large \code{for} loops.
}
\examples{
readData(parent.dir = "path/to/data/dir/", sheet.name = 1, file.name = "file.xlsx")
}

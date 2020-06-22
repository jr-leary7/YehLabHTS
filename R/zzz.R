#' Actions performed when `library(YehLabHTS)` is called.
#'
#' This file assigns our metadata and raw plate data objects to the global environment upon loading the package.

.onLoad <- function(libname, pkgname) {
  data("metadata_and_library_key",
       "raw_plates",
       package = pkgname,
       envir = parent.env(environment()))
  assign("results",
         results,
         envir = .GlobalEnv)
  assign("raw_list",
         raw_list,
         envir = .GlobalEnv)
  invisible()
 }

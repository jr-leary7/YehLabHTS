#' Actions performed when `library(YehLabHTS)` is called.
#'
#' This file assigns our metadata and raw plate data objects to the global environment upon loading the package.

.onLoad <- function(libname, pkgname) {
  data("normalized_data",
       package = pkgname,
       envir = parent.env(environment()))
  assign("drug_results", drug_results, .GlobalEnv)
  invisible()
 }

#' Read .xlsx HTS data
#'
#' This function reads in a .xlsx file, and outputs the name of the file being read before returning it. Useful for keeping track of data in large `for` loops.
#' @param parent.dir The directory within which the .xlsx file is located.
#' @param file.name The .xlsx file you wish to load.
#' @param col.names Should column names be loaded? Defaults to TRUE.
#' @import openxlsx
#' @export
#' @examples
#' readData(parent.dir = "path/to/data/dir/", file.name = "file.xlsx")

readData <- function(parent.dir, file.name, col.names = TRUE) {
  filepath <- sprintf("%s%s.xlsx",
                      parent.dir,
                      file.name)
  sprintf("Reading file: %s", file.name)
  df <- read.xlsx(filepath, colNames = col.names)
  return(df)
}

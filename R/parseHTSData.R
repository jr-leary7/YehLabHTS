#' Parse raw HTS plate data
#'
#' This functions imports HTS data in `.xlsx` format, parses it, calculates cell viability, and prepares the data for synergy analysis.
#' @param seed.use The random seed used to control stochasticity. Defaults to 629.
#' @import highSCREEN
#' @importFrom data.table rbindlist
#' @importFrom openxlsx read.xlsx
#' @export
#' @examples
#' parseHTSData()

parseHTSData <- function(seed.use = 629) {
  set.seed(629)  # not explicitly needed, but just in case
  # load metadata
  temp <- reformatMetadata()
  metadata <- temp[[1]]
  library_key <- temp[[2]]
  library_key_2016 <- temp[[3]]
  # load raw plate data into list
  HTS_data <- list()
  lib_sheets <- list()
  for (file in seq(metadata$Filename)) {
    # load raw HTS data
    raw_data <- readData(parent.dir = "./data/rawdata/",
                         file.name = metadata$Filename[file],
                         col.names = FALSE)
    library_sheet <- ifelse(metadata$LibrarySheet[1] == "LibraryKey",
                            "LibraryKey",
                            "LibraryKey2016")
    # save it for later -- might eventually wrap for loops in a subroutine then call lapply() on list
    HTS_data[[file]] <- raw_data
    lib_sheets[[file]] <- library_sheet
    # plate info
    plates_master <- rep(1:21, each = 16)
    plates_unique <- unique(plates_master)
    # normalize to DMSO on plate 1/2/3 depending on library drug dose

    # calculate viabilities

  }
}

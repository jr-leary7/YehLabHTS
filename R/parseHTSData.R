#' Parse raw HTS plate data
#'
#' This functions imports HTS data in `.xlsx` format, parses it, and prepares it for synergy analysis.
#' @param n.anchor.drugs The number of anchor drugs in the study. Defaults to 2.
#' @param plate.size The number of wells in each plate. Defaults to 384.
#' @param seed.use The random seed used to control stochasticity. Defaults to 629.
#' @import highSCREEN
#' @export

parseHTSData <- function(n.anchor.drugs = 2, plate.size = 384, seed.use = 629) {
  set.seed(629)
  metadata <- readData(parent.dir = "./data/",
                       file.name = "v3_YehLab_compound_library_synergy_screen_metadata",
                       col.names = TRUE)
  # remove unecessary metadata
  metadata <- metadata[-c(3, 4), ]
  HTS_data <- list()

  # create control map
  if (plate.size == 384) {
    control_map <- data.frame(X1 = c(rep("Control N", 8), rep("Control P", 8)),
                              X2 = rep("Control L", 16))
  }

  # load raw plate data into list
  for (file in seq(metadata$Filename)) {
    HTS_data[[file]] <- readData(parent.dir = "./data/rawdata/",
                                 file.name = metadata$Filename[file],
                                 sheet.name = "Merged 21 Plate Data File",
                                 col.names = FALSE)
  }
}

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
  set.seed(seed.use)  # not explicitly needed, but just in case
  # load metadata & library keys
  temp <- reformatMetadata()
  metadata <- temp[[1]]
  library_key <- temp[[2]]
  library_key_2016 <- temp[[3]]

  # load raw plate data into list & create matching library key list
  HTS_data <- list()
  lib_sheets <- list()
  for (file in seq(unique(metadata$Filename))) {
    # load raw HTS data
    raw_data <- readData(parent.dir = "./data/rawdata/",
                         file.name = metadata$Filename[file],
                         col.names = FALSE)
    # plate info
    raw_data$Plate_ID <- rep(NA, nrow(raw_data))
    raw_data$Plate_ID <- rep(1:21, each = 16)
    # little bit complicated if / else assignment due to how I'm using seq(), but I promise it works
    library_sheet <- ifelse(metadata[metadata$Filename == unique(metadata$Filename)[file], ]$LibrarySheet[1] == "LibraryKey",
                            "LibraryKey",
                            "LibraryKey2016")
    # save it for later
    HTS_data[[file]] <- raw_data
    lib_sheets[[file]] <- library_sheet
  }

  # create list containing library drug locations / doses & combine with anchor info
  drug_results <- list()
  for (i in seq(lib_sheets)) {
    if (lib_sheets[[i]] == "LibraryKey") {
      t <- library_key
    } else {
      t <- library_key_2016
    }
    # add anchor / score placeholder columns to library_key
    t$Anchor <- rep(NA, nrow(t))
    t$AnchorDose <- rep(NA, nrow(t))
    t$RawScore <- rep(NA, nrow(t))
    # create subset of metadata corresponding to anchor info raw `.xlsx` file
    metadata_sub <- metadata[metadata$Filename == unique(metadata$Filename)[i], ]
    # fill in anchor drug / dose / raw score
    for (j in seq(nrow(t))) {
      t$Anchor[j] <- metadata_sub[metadata_sub$Plate == t$Plate[j], ]$Anchor
      t$AnchorDose[j] <- metadata_sub[metadata_sub$Plate == t$Plate[j], ]$Dose
      HTS_sub <- HTS_data[[i]][HTS_data[[i]]$Plate_ID == t$Plate[j], ]
      t$RawScore[j] <- HTS_sub[t$Row[j], t$Column[j]]
    }
    drug_results[[i]] <- t
  }

  # normalize each raw value
  for (i in seq(drug_results)) {

  }
}








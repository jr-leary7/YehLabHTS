#' Parse raw HTS plate data
#'
#' This functions imports HTS data in `.xlsx` format, parses and normalizes it, calculates cell viability (bliss).
#' @param seed.use The random seed used to control stochasticity. Defaults to 629.
#' @importFrom data.table rbindlist
#' @importFrom openxlsx read.xlsx
#' @export
#' @examples
#' parseHTSData()

parseHTSData <- function(seed.use = 629) {
  set.seed(seed.use)  # not explicitly needed, but just in case
  # load metadata & library keys
  metadata <- results[[1]]
  library_key <- results[[2]]
  library_key_2016 <- results[[3]]

  # load raw plate data into list & create matching library key list
  HTS_data <- list()
  lib_sheets <- list()
  for (file in seq(unique(metadata$Filename))) {
    # load raw HTS data
    raw_data <- raw_list[[file]]
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
    # create temp var t to avoid editing original library key datasets
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

  # normalize each raw value according to DMSO on first plate of each
  for (i in seq(drug_results)) {
    # save drug_results[[i]] as temp var t (again) to avoid editing original data
    t <- drug_results[[i]]
    t$normalized <- NA  # pre-allocate length for memory reasons -- don't worry about it
    t$inhibition <- NA
    DMSO <- c(HTS_data[[1]][1:8, 1],
              HTS_data[[2]][1:8, 1],
              HTS_data[[3]][1:8, 1])
    mean_DMSO <- mean(DMSO)
    t$normalized <- t$RawScore / mean_DMSO * 100
    t$inhibition <- 100 - t$normalized
    drug_results[[i]] <- t
  }
}





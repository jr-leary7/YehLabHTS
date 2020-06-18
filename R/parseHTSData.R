#' Parse raw HTS plate data
#'
#' This functions imports HTS data in `.xlsx` format, parses it, and prepares it for synergy analysis.
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
  metadata <- readData(parent.dir = "./data/",
                       file.name = "v3_YehLab_compound_library_synergy_screen_metadata",
                       col.names = TRUE)
  # remove unecessary metadata
  metadata <- metadata[-c(3, 4), ]
  # load library sheets
  library_key <- readData(parent.dir = "./data/",
                          file.name = "v3_YehLab_compound_library_synergy_screen_metadata",
                          sheet.name = "LibraryKey",
                          col.names = TRUE)
  library_key_2016 <- readData(parent.dir = "./data/",
                               file.name = "v3_YehLab_compound_library_synergy_screen_metadata",
                               sheet.name = "LibraryKey2016",
                               col.names = TRUE)
  # reformat library datasets
  library_drugs <- unique(library_key$Compound)
  # create librarykey dataframe w/ compound, plate, row/col, & dose
  compound_list <- list()
  for (drug in seq(library_drugs)) {
    # make the generalizable later
    plates <- c(as.numeric(strsplit(library_key[drug, ]$Plate1, ",")[[1]]),
                as.numeric(strsplit(library_key[drug, ]$Plate2, ",")[[1]]),
                as.numeric(strsplit(library_key[drug, ]$Plate3, ",")[[1]]),
                as.numeric(strsplit(library_key[drug, ]$Plate4, ",")[[1]]),
                as.numeric(strsplit(library_key[drug, ]$Plate5, ",")[[1]]),
                as.numeric(strsplit(library_key[drug, ]$Plate6, ",")[[1]]))
    doses <- c(rep(library_key$`Dose1.(10uM)`[drug], length(strsplit(library_key$Plate1[1], ",")[[1]])),
               rep(library_key$`Dose2.(3uM)`[drug], length(strsplit(library_key$Plate1[2], ",")[[1]])),
               rep(library_key$`Dose3(1uM)`[drug], length(strsplit(library_key$Plate1[3], ",")[[1]])),
               rep(library_key$`Dose4.(300.nM)`[drug], length(strsplit(library_key$Plate1[4], ",")[[1]])),
               rep(library_key$`Dose5(100nM)`[drug], length(strsplit(library_key$Plate1[5], ",")[[1]])),
               rep(library_key$`Dose6(10nM)`[drug], length(strsplit(library_key$Plate1[6], ",")[[1]])))
    rows <- rep(library_key$Row[drug], length(doses))
    columns <- c(rep(library_key$Column1[drug], length(strsplit(library_key$Plate1[1], ",")[[1]])),
                 rep(library_key$Column2[drug], length(strsplit(library_key$Plate1[2], ",")[[1]])),
                 rep(library_key$Column3[drug], length(strsplit(library_key$Plate1[3], ",")[[1]])),
                 rep(library_key$Column4[drug], length(strsplit(library_key$Plate1[4], ",")[[1]])),
                 rep(library_key$Column5[drug], length(strsplit(library_key$Plate1[5], ",")[[1]])),
                 rep(library_key$Column6[drug], length(strsplit(library_key$Plate1[6], ",")[[1]])))
    compound <- rep(library_key$Compound[drug], length(doses))
    temp <- data.frame(Compound = compound,
                       Plate = plates,
                       Does = doses,
                       Row = rows,
                       Column = columns,
                       stringsAsFactors = FALSE)
    compound_list[[drug]] <- temp
  }
  # now coerce to a dataframe
  compound_df_librarykey <- rbindlist(compound_list)

  # repeat the previous operation for librarykey2016
  library_drugs <- unique(library_key_2016$Compound)
  compound_list <- list()
  for (drug in seq(library_drugs)) {
    # make this generalizable later
    plates <- c(as.numeric(strsplit(library_key_2016[drug, ]$Plate1, ",")[[1]]),
                as.numeric(strsplit(library_key_2016[drug, ]$Plate2, ",")[[1]]),
                as.numeric(strsplit(library_key_2016[drug, ]$Plate3, ",")[[1]]),
                as.numeric(strsplit(library_key_2016[drug, ]$Plate4, ",")[[1]]),
                as.numeric(strsplit(library_key_2016[drug, ]$Plate5, ",")[[1]]),
                as.numeric(strsplit(library_key_2016[drug, ]$Plate6, ",")[[1]]))
    doses <- c(rep(library_key_2016$Dose1[drug], length(strsplit(library_key_2016$Plate1[1], ",")[[1]])),
               rep(library_key_2016$Dose2[drug], length(strsplit(library_key_2016$Plate1[2], ",")[[1]])),
               rep(library_key_2016$Dose3[drug], length(strsplit(library_key_2016$Plate1[3], ",")[[1]])),
               rep(library_key_2016$Dose4[drug], length(strsplit(library_key_2016$Plate1[4], ",")[[1]])),
               rep(library_key_2016$Dose5[drug], length(strsplit(library_key_2016$Plate1[5], ",")[[1]])),
               rep(library_key_2016$Dose6[drug], length(strsplit(library_key_2016$Plate1[6], ",")[[1]])))
    columns <- c(rep(library_key_2016$Column1[drug], length(strsplit(library_key_2016$Plate1[1], ",")[[1]])),
                 rep(library_key_2016$Column2[drug], length(strsplit(library_key_2016$Plate1[2], ",")[[1]])),
                 rep(library_key_2016$Column3[drug], length(strsplit(library_key_2016$Plate1[3], ",")[[1]])),
                 rep(library_key_2016$Column4[drug], length(strsplit(library_key_2016$Plate1[4], ",")[[1]])),
                 rep(library_key_2016$Column5[drug], length(strsplit(library_key_2016$Plate1[5], ",")[[1]])),
                 rep(library_key_2016$Column6[drug], length(strsplit(library_key_2016$Plate1[6], ",")[[1]])))
    compound <- rep(library_key_2016$Compound[drug], length(doses))
    rows <- rep(library_key$Row[drug], length(doses))
    temp <- data.frame(Compound = compound,
                       Plate = plates,
                       Does = doses,
                       Row = rows,
                       Column = columns,
                       stringsAsFactors = FALSE)
    compound_list[[drug]] <- temp
  }
  # coerce to dataframe
  compound_df_librarykey2016 <- rbindlist(compound_list)

  # finally, create a long version of the metadata dataframe
  metadata_list <- list()
  for (drug in seq(nrow(metadata))) {
    # again, generalize this later
    plates <- c(as.numeric(strsplit(metadata$A1PlateRange0[drug], ",")[[1]]),
                as.numeric(strsplit(metadata$A1PlateRange1[drug], ",")[[1]]),
                as.numeric(strsplit(metadata$A1PlateRange2[drug], ",")[[1]]),
                as.numeric(strsplit(metadata$A1PlateRange3[drug], ",")[[1]]),
                as.numeric(strsplit(metadata$A2PlateRange0[drug], ",")[[1]]),
                as.numeric(strsplit(metadata$A2PlateRange1[drug], ",")[[1]]),
                as.numeric(strsplit(metadata$A2PlateRange2[drug], ",")[[1]]),
                as.numeric(strsplit(metadata$A2PlateRange3[drug], ",")[[1]]))
    doses <- c(rep(metadata$A1Dose0[drug], length(strsplit(metadata$A1PlateRange0[drug], ",")[[1]])),
               rep(metadata$A1Dose1[drug], length(strsplit(metadata$A1PlateRange1[drug], ",")[[1]])),
               rep(metadata$A1Dose2[drug], length(strsplit(metadata$A1PlateRange2[drug], ",")[[1]])),
               rep(metadata$A1Dose3[drug], length(strsplit(metadata$A1PlateRange3[drug], ",")[[1]])),
               rep(metadata$A2Dose0[drug], length(strsplit(metadata$A2PlateRange0[drug], ",")[[1]])),
               rep(metadata$A2Dose1[drug], length(strsplit(metadata$A2PlateRange1[drug], ",")[[1]])),
               rep(metadata$A2Dose2[drug], length(strsplit(metadata$A2PlateRange2[drug], ",")[[1]])),
               rep(metadata$A2Dose3[drug], length(strsplit(metadata$A2PlateRange3[drug], ",")[[1]])))
    anchors <- c(rep(metadata$Anchor1[drug], length(doses) / 2),
                 rep(metadata$Anchor2[drug], length(doses) / 2))
    filenames <- rep(metadata$Filename[drug], length(doses))
    lib_sheets <- rep(metadata$LibrarySheet[drug], length(doses))
    samp_info <- rep(metadata$SampleInfo[drug], length(doses))
    temp <- data.frame(Filename = filenames,
                       SampleInfo = samp_info,
                       LibrarySheet = lib_sheets,
                       Anchor = anchors,
                       Dose = doses,
                       Plate = plates,
                       stringsAsFactors = FALSE)
    metadata_list[[drug]] <- temp
  }
  # coerce to dataframe
  metadata_df <- rbindlist(metadata_list)

  # load raw plate data into list
  HTS_data <- list()
  lib_sheets <- list()
  for (file in seq(metadata$Filename)) {
    # load raw HTS data
    raw_data <- readData(parent.dir = "./data/rawdata/",
                         file.name = metadata$Filename[file],
                         col.names = FALSE)
    library_sheet <- ifelse(metadata$LibrarySheet[file] == "LibraryKey",
                            "LibraryKey",
                            "LibraryKey2016")
    # save it for later -- might eventually wrap for loops in a subroutine then call lapply() on list
    HTS_data[[file]] <- raw_data
    lib_sheets[[file]] <- library_sheet
    # plate info
    plates_master <- rep(1:21, each = 16)
    plates_unique <- unique(plates_master)


  }
}

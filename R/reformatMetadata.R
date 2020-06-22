#' Reformat metadata
#'
#' This function reformats the metadata and library sheets from out HTS experiments into a long format.
#' @importFrom openxlsx read.xlsx
#' @importFrom data.table rbindlist
#' @export
#' @examples
#' reformatMetadata()

reformatMetadata <- function() {
  metadata <- readData(parent.dir = "./data/",
                       file.name = "v3_YehLab_compound_library_synergy_screen_metadata",
                       col.names = TRUE)
  # remove unecessary metadata
  metadata <- metadata[-c(3, 4), ]  # files we don't want
  metadata <- na.omit(metadata)  # necessary to parse correctly
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
    # 7392 wells per `xlsx` file --> 336*24 - 216*21 = 7392
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
    temp <- data.table(Compound = compound,
                       Plate = plates,
                       Dose = doses,
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
    # make this generalizable to a different plate size / experimental design later
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
    temp <- data.table(Compound = compound,
                       Plate = plates,
                       Dose = doses,
                       Row = rows,
                       Column = columns,
                       stringsAsFactors = FALSE)
    compound_list[[drug]] <- temp
  }
  # coerce to dataframe
  compound_df_librarykey2016 <- rbindlist(compound_list)

  # finally, create a long version of the metadata dataframe
  metadata_list <- list()
  for (anchor in seq(nrow(metadata))) {
    # plates 1:3 always have dose=0 for anchor 1 & 2 (in this experiment)
    plates <- c(as.numeric(strsplit(metadata$A1PlateRange0[anchor], ",")[[1]]),
                as.numeric(strsplit(metadata$A1PlateRange1[anchor], ",")[[1]]),
                as.numeric(strsplit(metadata$A1PlateRange2[anchor], ",")[[1]]),
                as.numeric(strsplit(metadata$A1PlateRange3[anchor], ",")[[1]]),
                as.numeric(strsplit(metadata$A2PlateRange1[anchor], ",")[[1]]),
                as.numeric(strsplit(metadata$A2PlateRange2[anchor], ",")[[1]]),
                as.numeric(strsplit(metadata$A2PlateRange3[anchor], ",")[[1]]))
    doses <- c(rep(metadata$A1Dose0[anchor], length(strsplit(metadata$A1PlateRange0[anchor], ",")[[1]])),
               rep(metadata$A1Dose1[anchor], length(strsplit(metadata$A1PlateRange1[anchor], ",")[[1]])),
               rep(metadata$A1Dose2[anchor], length(strsplit(metadata$A1PlateRange2[anchor], ",")[[1]])),
               rep(metadata$A1Dose3[anchor], length(strsplit(metadata$A1PlateRange3[anchor], ",")[[1]])),
               rep(metadata$A2Dose1[anchor], length(strsplit(metadata$A2PlateRange1[anchor], ",")[[1]])),
               rep(metadata$A2Dose2[anchor], length(strsplit(metadata$A2PlateRange2[anchor], ",")[[1]])),
               rep(metadata$A2Dose3[anchor], length(strsplit(metadata$A2PlateRange3[anchor], ",")[[1]])))
    anchors <- c(rep("Anchor doses = 0", 3),
                 rep(metadata$Anchor1[anchor], (length(doses) - 3) / 2),
                 rep(metadata$Anchor2[anchor], (length(doses) - 3) / 2))
    filenames <- rep(metadata$Filename[anchor], length(doses))
    lib_sheets <- rep(metadata$LibrarySheet[anchor], length(doses))
    samp_info <- rep(metadata$SampleInfo[anchor], length(doses))
    temp <- data.table(Filename = filenames,
                       SampleInfo = samp_info,
                       LibrarySheet = lib_sheets,
                       Anchor = anchors,
                       Dose = doses,
                       Plate = plates,
                       stringsAsFactors = FALSE)
    metadata_list[[anchor]] <- temp
  }
  # coerce to dataframe
  metadata_df <- rbindlist(metadata_list)

  # save raw plate data in a list so that it can be a global environment variable later
  raw_list <- list()
  for (file in unique(metadata$Filename)) {
    raw_list[[file]] <- readData(parent.dir = "./data/rawdata/",
                                 file.name = file,
                                 col.names = FALSE)
  }
  # save raw plate data
  names(raw_list) <- unique(metadata$Filename)
  save(raw_list, file = "./data/raw_plates.RData")

  # save metadata and library_key sheets
  results <- list(metadata_df, compound_df_librarykey, compound_df_librarykey2016)
  save(results, file = "./data/metadata_and_library_key.RData")

  return("Parsing completed successfully !")
}

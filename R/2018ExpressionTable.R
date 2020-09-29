library(data.table)
library(tidyverse)
library(data.table)
library(YehLabHTS)
library(dplyr)
setwd("~/YehLabHTS")

seqMetadata <- readData(parent.dir = "./data/",
                          file.name = "Sequenced_PDX_CAF_lines",
                          col.names = TRUE)


# cellLinesExact <- c(P100422T1p7m879, P140227T1, P100422T1p7m299, P100422T1pm862, P130411T1p2m862, P140710N1, P170119T1)

# P140710N1 was also in the screen but removed from analysis due to
# a cell line mix-up: P140710N1 matches P140227T1 for STR
cellLines <- c("P100422", "P130411T1", "P140227", "P170119T1")
#filter the metadata file for PDXcell type
seqMetadata <- filter(seqMetadata, Line %in% cellLines, Type == "PDXcells")

#removing a seq file on an organoid
seqMetadata <- seqMetadata[-c(5),]

filenames <- seqMetadata$fpkm_file_header
#expression subset
expSubset <- subset(Yeh_Salmon[["ex"]], select = c(filenames))


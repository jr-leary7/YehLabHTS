library(YehLabHTS)
library(data.table) #to run reformatMetadata
library(tidyverse)

setwd("~/YehLabHTS")

classMetadata <- readData(parent.dir = "./data/",
                     file.name = "Combined Screen Compound List 160912",
                     col.names = TRUE)

kinaseMetadata <- filter(classMetadata, Class == "Kinase")
libraryKinases <- as.vector(kinaseMetadata$Product.Name)
anchorKinases <- c("Dasatinib", "MLN8237", "Dinaciclib", "Anchor doses = 0")

# 176 libdrugs x 42 plates (2 lib drugs per plate, 21 plates ) = 7,392
# 97 kinase drugs x 2 per plate x 21 = 4072
# ISSUE:
# there is a space after the drug: LDN-212854 in the original metadata sheet
# not properly getting indexed even when a space is added to the "classMetadata" sheet
# drug_results[[1]][3152,]
# kinaseMetadata[51,]
# setdiff(libraryKinases, kinase_results[[3]]$Compound)

# kinase_results is a list of viability and associated data for kinase inhibitors only
kinase_results <- list()
for (i in 1:length(drug_results)) {

  kinase_results[[i]] <- drug_results[[i]] %>%
    filter(Compound %in% libraryKinases & Anchor %in% anchorKinases)
}



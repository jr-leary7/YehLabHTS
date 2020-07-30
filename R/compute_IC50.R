library(YehLabHTS)
library("data.table")
library(drexplorer)
library(reshape2)
library(ggplot2)
library("tidyverse")
library(anchors)

calculateIC50 <- function(drug_results) {
  ic_results <- data.frame()
  #every library and anchor drug combination
  #pairs <- unique(paste(drug_results[[h]]$Compound, drug_results[[h]]$Anchor))


  for(h in 1:length(drug_results)) {
    #anch drugs needs to start at 2 to avoid "anchor doses = 0"
    for (i in 2:length(anchor_drugs)) {
      for(j in 1:length(library_drugs)) {
        library_drugs <- unique(drug_results[[h]]$Compound)
        anchor_drugs <- unique(drug_results[[h]]$Anchor)

        #makes a subset for each unique anchor/lib pair
        unique_pair_subset <- drug_results[[h]][ which(drug_results[[h]]$Compound == library_drugs[[j]] & drug_results[[h]]$Anchor == anchor_drugs[[i]])]
        unique_pair_subset <- setorder(unique_pair_subset, Dose, AnchorDose)
        compound_dose <- unique(unique_pair_subset$Dose)
        anchor_dose <- unique(unique_pair_subset$AnchorDose)

        #matrix for one anchor/lib pair with compound doses along the y, anchor doses on the x, as a funct of percent inhibition
        inhibitionMatrix <- matrix(as.numeric(unique_pair_subset$inhibition), length(anchor_dose), length(compound_dose), byrow = F)

        compound_dose <- append(compound_dose, 0)
        compound_dose <- sort(compound_dose)
        anchor_dose <- append(anchor_dose, 0)
        anchor_dose <- sort(anchor_dose)

        #add in a row to the matrix for the DMSO anchor dose
        zero_dose_anchor <- drug_results[[h]][ which(drug_results[[h]]$Compound == library_drugs[[j]] & drug_results[[h]]$Anchor == "Anchor doses = 0")]
        zero_dose_anchor <- zero_dose_anchor[order(Dose),]
        inhibitionMatrix <- rbind(zero_dose_anchor$inhibition, inhibitionMatrix)
        inhibitionMatrix <- cbind(0, inhibitionMatrix)

        colnames(inhibitionMatrix) <- compound_dose
        rownames(inhibitionMatrix) <- anchor_dose

        #working calculates IC50 of lib compound only (anchor dose=0)
        fit <- drFit(drMat = cbind(as.numeric(colnames(inhibitionMatrix)), 1-inhibitionMatrix[1, ]/100), modelName = "sigEmax", alpha=1, fitCtr=T)
        ic50 <- as.numeric(computeIC(fit, log.d = FALSE, interpolation = F))
        ic25 <- as.numeric(computeIC(fit, percent = 0.25, log.d = FALSE, interpolation = F))
        if(all(inhibitionMatrix[1, ]/100) > .5) ic50 = NA
        if(all(inhibitionMatrix[1, ]/100) > .5) ic25 = NA
        #print(ic50)

        #original
        #calculates IC50 of lib compound only (anchor dose=0)
        # fit <- drFit(drMat = cbind(as.numeric(colnames(inhibitionMatrix)[c(1:7)]),1-inhibitionMatrix[1, c(1:7)]/100), modelName = "sigEmax", alpha=1, fitCtr=T)
        # ic50 <- as.numeric(computeIC(fit, log.d=FALSE, interpolation=F))
        # if(all(1-inhibitionMatrix[1,c(1:7)]/100 > 0.5)) ic50 = NA
        # print(ic50)

        ic_results <- rbind(ic_results, c(anchor_drugs[[i]], library_drugs[[j]], ic50, ic25))

      }
    }
  }

}
colnames(ic_results) <- c("Anchor Drug","Library Drug","ic50", "ic25")
save(ic_results, file = "ic_results.Rdata")
load("ic_results.Rdata")

######plotting and visualizing IC50 ranges

#total number of IC50 values calculated for 1 drug file (102)
#352 - sum(is.na(ic50_results$ic50))

#range is 1.1 e-7 to Inf
#range(ic50_results$ic50, na.rm = T)

#replace the Inf values to NA
ic_results <- replace.value(ic_results, c("ic50","ic25"), from = Inf, to= NA, verbose = FALSE)
ic_results <- as.numeric(ic_results)
#ic50Data is ic50 results with only numbers (removes Inf and NA)
ic50Data <- ic_results[!is.na(as.numeric(as.character(ic_results$ic50))),]
#ic50Data <- ic50Data[!grepl("Inf", ic50Data$ic50),]
ic50Data <- setorder(ic50Data, ic50)
ic50Data <- as.numeric(ic50Data$ic50)
is.na(ic_results) <-sapply(ic_results, is.infinite)
#ic50Data <- as.data.frame(ic50Data)

#make all 0s into NA
ic_results[ic_results == 0] <- NA
ic25Data <- ic_results[ic25]
#remove NA values
ic25Data <- ic_results$ic25[complete.cases(ic_results$ic25),]
ic25Data <- ic_results[!is.na(as.numeric(as.character(ic_results$ic25))),]
#ic25Data <- ic25Data[!grepl("Inf", ic25Data$ic50),]
ic25Data <- setorder(ic25Data, ic25)
ic25Data <- as.numeric(ic25Data$ic25)
is.na(ic_results) <-sapply(ic_results, is.infinite)

ic50_results2 <- ic50_results2 %>%
  group_by(`Library Drug`) %>%
  max(ic50)

#total lib/anchor unique combinations for 15 files = 5280
# 2586
#histogram of all combined files and their calculatable IC50s (2586)
hist(ic50Data, freq = T, xlim = c(1e-9, 1.25e-5), breaks = 1000)
abline(v = median(ic50Data), col = "red", lwd = 2)

hist(ic25Data, freq = T, xlim = c(-1e-9, 5e-7), breaks = 100)
abline(v = median(ic50Data), col = "red", lwd = 2)

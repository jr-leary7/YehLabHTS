library(YehLabHTS)
library("data.table")
library(drexplorer)
library(reshape2)
library(ggplot2)
library("tidyverse")

calculateIC50 <- function(drug_results) {
  ic50_results <- data.frame()
  for(h in 1:length(drug_results)) {
    #anch drugs needs to start at 2 to avoid "anchor doses = 0"
    for (i in 2:length(anchor_drugs)) {
      for(j in 1:length(library_drugs)) {
        #every library and anchor drug combination
        pairs <- unique(paste(drug_results[[h]]$Compound, drug_results[[h]]$Anchor))
        library_drugs <- unique(drug_results[[h]]$Compound)
        anchor_drugs <- unique(drug_results[[h]]$Anchor)

        #makes a subset for each unique anchor/lib pair
        unique_pair_subset <- drug_results[[h]][ which(drug_results[[h]]$Compound == library_drugs[[j]] & drug_results[[h]]$Anchor == anchor_drugs[[i]])]
        unique_pair_subset <- setorder(unique_pair_subset, Dose, AnchorDose)
        compound_dose <- unique(unique_pair_subset$Dose)
        anchor_dose <- unique(unique_pair_subset$AnchorDose)

        #matrix for one anchor/lib pair with compound doses along the y, anchor doses on the x
        viabilityMatrix <- matrix(as.numeric(unique_pair_subset$inhibition), length(anchor_dose), length(compound_dose), byrow = F)

        compound_dose <- append(compound_dose, 0)
        compound_dose <- sort(compound_dose)
        anchor_dose <- append(anchor_dose, 0)
        anchor_dose <- sort(anchor_dose)

        #add in a row to the matrix for the DMSO anchor dose
        zero_dose_anchor <- drug_results[[h]][ which(drug_results[[h]]$Compound == library_drugs[[j]] & drug_results[[h]]$Anchor == "Anchor doses = 0")]
        zero_dose_anchor <- zero_dose_anchor[order(Dose),]
        viabilityMatrix <- rbind(zero_dose_anchor$inhibition, viabilityMatrix)
        viabilityMatrix <- cbind(0, viabilityMatrix)

        colnames(viabilityMatrix) <- compound_dose
        rownames(viabilityMatrix) <- anchor_dose

        #calculates IC50 of lib compound only (anchor dose=0)
        fit <- drFit(drMat = cbind(as.numeric(colnames(viabilityMatrix)[c(1:7)]),1-viabilityMatrix[1, c(1:7)]/100), modelName = "sigEmax", alpha=1, fitCtr=T)
        ic50 <- as.numeric(computeIC(fit, log.d=FALSE, interpolation=F))
        if(all(1-viabilityMatrix[1,c(1:7)]/100 > 0.5)) ic50 = NA
        print(ic50)

        ic50_results <- rbind(ic50_results, c(anchor_drugs[[i]], library_drugs[[j]], ic50))
        colnames(ic50_results) <- c("Anchor Drug","Library Drug","ic50")
      }
    }
  }
}


######plotting and visualizing IC50 ranges
#for each file there are 2 of the same library drugs

#total number of IC50 values calculated for 1 drug file (102)
352 - sum(is.na(ic50_results$ic50))

#range is 1.1 e-7 to Inf
range(ic50_results$ic50, na.rm = T)

#ic50Data is ic50 results with only numbers (removes Inf and NA)
ic50Data <- ic50_results[!is.na(as.numeric(as.character(ic50_results$ic50))),]
ic50Data <- ic50Data[!grepl("Inf", ic50Data$ic50),]
ic50Data <- setorder(ic50Data, ic50)
ic50Data <- as.numeric(ic50Data$ic50)
ic50Data <- as.data.frame(ic50Data)

hist(ic50Data, freq = T, xlim = c(0, 1e-05), breaks = 20)
abline(v = median(ic50Data), col = "red", lwd = 2)



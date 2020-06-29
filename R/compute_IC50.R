calculateIC50 <- function(drug_results) {
  ic50_results <- data.frame()
  #every library and anchor drug combination
  pairs <- unique(paste(drug_results[[1]]$Compound, drug_results[[1]]$Anchor))
  library_drugs <- unique(drug_results[[1]]$Compound)
  anchor_drugs <- unique(drug_results[[1]]$Anchor)

  for(h in 2:length(drug_results)) {
    #anch drugs needs to start at 2 to avoid "anchor doses = 0"
    for (i in 2:length(anchor_drugs)) {
      for(j in 1:length(library_drugs)) {
        #makes a subset for each unique anchor/lib pair
        subset <- drug_results[[h]][ which(drug_results[[h]]$Compound == library_drugs[[j]] & drug_results[[h]]$Anchor == anchor_drugs[[i]])]
        subset <- setorder(subset, Dose, AnchorDose)
        compound_dose <- unique(subset$Dose)
        anchor_dose <- unique(subset$AnchorDose)

        #matrix for one anchor/lib pair with compound doses along the y, anchor doses on the x
        matrix <- matrix(as.numeric(subset$inhibition), length(anchor_dose), length(compound_dose), byrow = F)

        compound_dose <- append(compound_dose, 0)
        compound_dose <- sort(compound_dose)
        anchor_dose <- append(anchor_dose, 0)
        anchor_dose <- sort(anchor_dose)

        #add in a row to the matrix for the DMSO anchor dose
        zero_dose_anchor <- drug_results[[h]][ which(drug_results[[h]]$Compound == library_drugs[[j]] & drug_results[[h]]$Anchor == "Anchor doses = 0")]
        zero_dose_anchor <- zero_dose_anchor[order(Dose),]
        matrix <- rbind(zero_dose_anchor$inhibition, matrix)
        matrix <- cbind(0, matrix)

        colnames(matrix) <- compound_dose
        rownames(matrix) <- anchor_dose

        #calculates IC50 of lib compound only (anchor dose=0)
        fit <- drFit(drMat = cbind(as.numeric(colnames(matrix)[c(1:7)]),1-matrix[1, c(1:7)]/100), modelName = "sigEmax", alpha=1, fitCtr=T)
        ic50 <- as.numeric(computeIC(fit, log.d=FALSE, interpolation=F))
        if(all(1-matrix[1,c(1:7)]/100 > 0.5)) ic50 = NA
        print(ic50)

        ic50_results <- rbind(ic50_results, c(anchor_drugs[[i]], library_drugs[[j]], ic50))
        colnames(ic50_results) <- c("Anchor Drug","Library Drug","ic50")
      }
    }
  }
}

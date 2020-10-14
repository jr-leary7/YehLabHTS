library(data.table)
library(tidyverse)
library(data.table)
library(dplyr)
library(projectR) #PCA, NMF
library(ggplot2)
library(gplots)
library('CoGAPS') #NMF
library(scFeatureFilter)
library(NormExpression)
library(dendextend)
library(SCnorm)
library(edgeR)

setwd("~/YehLabHTS")
load("/Volumes/Jen Jen Yeh Lab/RNAseq/RData/Yeh_Salmon.RData")

seqMetadata <- readData(parent.dir = "./data/",
                          file.name = "Sequenced_PDX_CAF_lines",
                          col.names = TRUE)

# P140710N1 was also in the screen but removed from analysis due to
# a cell line mix-up: P140710N1 matches P140227T1 for STR
cellLines <- c("P100422", "P130411T1", "P140227", "P170119T1")

#filter the metadata file for PDXcell type
invitroMetadata <- filter(seqMetadata, Line %in% cellLines, Type == "PDXcells")
invivoMetadata <- filter(seqMetadata, Line %in% cellLines)

#removing a seq file on an organoid
invitroMetadata <- invitroMetadata[-c(5),]

CreateExpressionMatrix <- function(Metadata) {
  filenames <- Metadata$fpkm_file_header

  #expression subset
  expSubset <- subset(Yeh_Salmon[["ex"]], select = c(filenames))

  # open brackets retain the row names (genes)
  #ExpSubset[] <- sapply(expSubset[], as.numeric)
}

invitroExpression <- CreateExpressionMatrix(invitroMetadata)
invitroExpression <- data.matrix(invitroExpression, rownames.force = NA)
#invitroExpression[] <- sapply(invitroExpression[], as.numeric)

# log transform data
loginvitroExpression <- log(invitroExpression +1)
plot(loginvitroExpression)

##########
# NMF using CoGAPS

NMFResult <- function(expressionMatrix){
  # create model parameters object
  params <- new("CogapsParams")

  # set the value for a specific parameter
  params <- setParam(params, "nPatterns", 2)

  # run CoGAPS with specified model parameters
  # cogapsResult: sampleFactors = A matrix, featureLoadings = P matrix
  result <- CoGAPS(expressionMatrix, params, nIterations=1000)
}
invitroResult <- NMFResult(loginvitroExpression)
# plots the Cogapsresult
plot(invitroResult)

# GAPS function isn't found ?
#resultGAPS <- GAPS(loginvitroExpression, unc=0.01, isPercentError=FALSE, numPatterns = 3)

# heatmap crashes R- probably too large
#plotPatternMarkers(result, loginvitroExpression)

#########
# ProjectR, returns projectionPatterns = relative weights for samples in feature space

pc.loginvitroExpression <- prcomp(t(loginvitroExpression))
pcVAR <- round(((pc.loginvitroExpression$sdev)^2/sum(pc.loginvitroExpression$sdev^2))*100,2)
dPCA <- data.frame(cbind(pc.loginvitroExpression$x,pc.loginvitroExpression))
# ISSUE: getting dPCA into continuous values. transforming/lappy is freezing R Studio
transform(dPCA, PC1 = as.numeric(levels(PC1)))

pPCA <- ggplot(dPCA, aes(x = PC1, y = PC2)) +
  geom_point()
pPCA


# can't get this first part to work for projectR: Obtaining CoGAPS patterns to project.
AP <- get(as.matrix(expSubset)) #CoGAPS run data
AP <- expSubset$mean

pNMF <- heatmap.2(as.matrix(logexpSubset), Rowv = FALSE)

# data/loadings names = rownames(data/loadings), NP = NULL to use full matrix, full =false to return only the new pattern object
# I think loadings is the std matrix- what is used in CoGAPS
PCAexpSubset <- projectR(data = expSubset, loadings = pc.expSubset, dataNames = rownames(data), loadingsNames = rownames(loadings), NP = NULL, full = false)

group <- c('')
d <- DGEList(counts = expSubset, group = group)
TMM <- calcNormFactors(expSubset, method="TMM")
expSubset.AUCVCs1 <- gridAUCVC(data = expSubset, dataType = "bk", HG7=
                                 bkRNA18_factors$HG7, ERCC= bkRNA18_factors$ERCC, TN=bkRNA18_factors$TN,
                               TC=bkRNA18_factors$TC, CR=bkRNA18_factors$CR, NR=bkRNA18_factors$NR,
                               DESeq=bkRNA18_factors$DESeq, UQ=bkRNA18_factors$UQ,
                               TMM=bkRNA18_factors$TMM, TU= 1, GAPDH=bkRNA18_factors$GAPDH,
                               nonzeroRatios= c(0.7, 0.8, 0.9, 1));  )

expSubset.cvs <- calculate_cvs(expSubset, max_zeros = 1)
nIter <- 5000
results <- CoGAPS(expSubset, expSubset,
                  +GStoGenes=GSets,
                  + nFactor=3,
                  + nEquil=nIter, nSample=nIter,
                  + plot=FALSE)

######
# ProjectR example data- row= samples, col= genes
data(p.RNAseq6l3c3t)
#principle components of dataset
pc.RNAseq6l3c3t<-prcomp(t(p.RNAseq6l3c3t))
#variance of the PC
pcVAR <- round(((pc.RNAseq6l3c3t$sdev)^2/sum(pc.RNAseq6l3c3t$sdev^2))*100,2)
#x= transposed sample=rows, PCs = col
dPCA <- data.frame(cbind(pc.RNAseq6l3c3t$x,pd.RNAseq6l3c3t))
#plot pca

setCOL <- scale_colour_manual(values = c("blue","black","red"), name="Condition:")
setFILL <- scale_fill_manual(values = c("blue","black","red"),guide = FALSE)
setPCH <- scale_shape_manual(values=c(23,22,25,25,21,24),name="Cell Line:")
pPCA <- ggplot(dPCA, aes(x=PC1, y=PC2, colour=ID.cond, shape=ID.line,
                         fill=ID.cond)) +
  geom_point(aes(size=days),alpha=.6)+
  setCOL + setPCH + setFILL +
  scale_size_area(breaks = c(2,4,6), name="Day") +
  theme(legend.position=c(0,0), legend.justification=c(0,0),
        legend.direction = "horizontal",
        panel.background = element_rect(fill = "white",colour=NA),
        legend.background = element_rect(fill = "transparent",colour=NA),
        plot.title = element_text(vjust = 0,hjust=0,face="bold")) +
  labs(title = "PCA of hPSC PolyA RNAseq",
       x=paste("PC1 (",pcVAR[1],"% of varience)",sep=""),
       y=paste("PC2 (",pcVAR[2],"% of varience)",sep=""))

data(p.ESepiGen4c1l)
PCA2ESepi <- projectR(data = p.ESepiGen4c1l$mRNA.Seq,loadings=pc.RNAseq6l3c3t,
                      full=TRUE, dataNames=map.ESepiGen4c1l[["GeneSymbols"]])

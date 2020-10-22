library(xlsx)
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
library(fastICA)

setwd("~/YehLabHTS")
load("/Volumes/Jen Jen Yeh Lab/RNAseq/RData/Yeh_Salmon.RData")
load("/Volumes/Jen Jen Yeh Lab/RNAseq/RData/Yeh_Salmon_caf_naf_with_fails.RData")


# P140710N1 was also in the screen but removed from analysis due to
# a cell line mix-up: P140710N1 matches P140227T1 for STR
cellLines <- c("100422", "130411", "140227", "170119")
treatment <- c("untreated", "")
type <- c("PDXcells")

CreateExpressionMatrix <- function(Lines, sampInfo, type, treatment, RData) {
   combinedFilenames <- c()
   filename <- c()
   for (i in 1:length(Lines)){
     # subset Yeh_Salmon samp info by cell line
     lineSubset <- sampInfo[grep(Lines[[i]], sampInfo$filename), ]
     # filter by untreated and unmixed samples
     lineSubset <- filter(lineSubset, Type %in% type, Treatment_orig %in% treatment, Treatment1 %in% treatment, line != "mix")

     # filename for one line
     filename <- as.character(lineSubset$filename)

     # list of all filenames for every line
     combinedFilenames <- c(combinedFilenames, filename)

     # index by filename to create expression subset (expSub)
     expSub <- subset(RData[["ex"]], select = c(combinedFilenames))

   }
   return(expSub)
}
# Need to remove scramble sgRNA samples *
invitroExpression <- CreateExpressionMatrix(Lines = cellLines, sampInfo = Yeh_Salmon$sampInfo,
                                            type = c("CAF", "PDXcells"), treatment = treatment, RData = Yeh_Salmon)

# expression dataframe of the lines in any type (ie human, organoid, pdx, cell line...)
totalExpression <- CreateExpressionMatrix(Lines = cellLines, sampInfo = Yeh_Salmon$sampInfo,
                                          type = unique(Yeh_Salmon[["sampInfo"]]$Type), treatment = treatment, RData = Yeh_Salmon)

# cell line exp data from another R data (Yeh_Salmon_caf_naf_with_fails)
invitroCAFNAFSalmon <- CreateExpressionMatrix(Lines = cellLines, sampInfo = Yeh_Salmon_caf_naf_with_fails$sampInfo,
                                              type = c("CAF", "PDXcells"), treatment = treatment, RData = Yeh_Salmon_caf_naf_with_fails)

# final dataframe of all cell line expression data
invitroExpression <- cbind(invitroExpression, invitroCAFNAFSalmon)

# total is no different than invitro
# totalCAFNAFSalmon <- CreateExpressionMatrix(Lines = cellLines, sampInfo = Yeh_Salmon_caf_naf_with_fails$sampInfo,
#                                          type = unique(Yeh_Salmon_caf_naf_with_fails[["sampInfo"]]$Type), treatment = treatment,
#                                          RData = Yeh_Salmon_caf_naf_with_fails)


invitroExpression <- data.matrix(invitroExpression, rownames.force = NA)
#invitroExpression[] <- sapply(invitroExpression[], as.numeric)

# log transform data
loginvitroExpression <- log(invitroExpression +1)
# plot(loginvitroExpression)

# normalization
lib.size <- sum(d$sample1)
scale.factors <- calcNormFactors(loginvitroExpression, method = "TMM")
RPKM <- rpkm(scale.factors)
norm.data <- t(t(TMM$sample1)/(scale.factors*lib.size))

##########
# NMF using CoGAPS

NMFResult <- function(expressionMatrix){
  # create model parameters object
  params <- new("CogapsParams")

  # set the value for a specific parameter
  params <- setParam(params, "nPatterns", 3)

  # run CoGAPS with specified model parameters
  # cogapsResult: sampleFactors = A matrix, featureLoadings = P matrix
  result <- CoGAPS(expressionMatrix, params, nIterations=1000)
}

invitroResult <- NMFResult(loginvitroExpression)
# plots the Cogapsresult. Samples 1-5 are P100422 replicates, samples 6-8= P130411, 9-10= 140227, 11 = 170119
# Pattern 1 = shared between all data sets?
dev.off()
plot(invitroResult)
axis(labels = TRUE, side = 1)
ggplot(invitroResult) + geom_point()

# GAPS function isn't found ?
#resultGAPS <- GAPS(loginvitroExpression, unc=0.01, isPercentError=FALSE, numPatterns = 3)

# heatmap crashes R- probably too large
#plotPatternMarkers(result, loginvitroExpression)

#########
# ProjectR, PCA plot
# returns projectionPatterns = relative weights for samples in feature space
CreateProjectionPattern <- function(expressionMatrix) {
  pc.expressionMatrix <- prcomp(t(expressionMatrix))
  pcVAR <- round(((pc.expressionMatrix$sdev)^2/sum(pc.expressionMatrix$sdev^2))*100,2)
  dPCA <- data.frame(pc.expressionMatrix$x)
}

invitroPCA <- CreateProjectionPattern(loginvitroExpression)

# PCA plot. The 4 replicate 130411 cell lines don't cluster like I would expect
pPCA <- ggplot(invitroPCA, aes(x = PC1, y = PC2)) +
  geom_point() + geom_text(aes(label=row.names(invitroPCA)))
pPCA

#####
# ICA plots using fastICA
# X = expmatrix, n.comp = # of components to extract
invitroICA <- fastICA(loginvitroExpression, n.comp = 2, alg.typ = "parallel", fun = "logcosh", alpha = 1,
             method = "R", row.norm = FALSE, maxit = 200,
             tol = 0.0001, verbose = TRUE)
par(mfrow = c(1, 3))
plot(invitroICA$X, main = "Pre-processed data")
plot(invitroICA$X %*% invitroICA$K, main = "PCA components")
plot(invitroICA$S, main = "ICA components")

#####
# sample data for ICA
S <- matrix(runif(10000), 5000, 2)
A <- matrix(c(1, 1, -1, 3), 2, 2, byrow = TRUE)
X <- S %*% A
a <- fastICA(X, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1,
             method = "C", row.norm = FALSE, maxit = 200,
             tol = 0.0001, verbose = TRUE)
par(mfrow = c(1, 3))
plot(a$X, main = "Pre-processed data")
plot(a$X %*% a$K, main = "PCA components")
plot(a$S, main = "ICA components")

#####
# can't get this projectR for NMF to work
# stuck on first part: Obtaining CoGAPS patterns to project.
AP <- get(as.matrix(expSubset)) #CoGAPS run data
AP <- expSubset$mean

pNMF <- heatmap.2(as.matrix(logexpSubset), Rowv = FALSE)

# data/loadings names = rownames(data/loadings), NP = NULL to use full matrix, full =false to return only the new pattern object
# I think loadings is the std matrix- what is used in CoGAPS
PCAexpSubset <- projectR(data = expSubset, loadings = pc.expSubset, dataNames = rownames(data), loadingsNames = rownames(loadings), NP = NULL, full = false)

group <- c('')
d <- DGEList(counts = expSubset, group = group)

TMM <- calcNormFactors(loginvitroExpression, method="TMM")
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
adPCA <- data.frame(cbind(pc.RNAseq6l3c3t$x,pd.RNAseq6l3c3t))
bdPCA <- data.frame(cbind(pc.RNAseq6l3c3t$x))
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
pPCA
data(p.ESepiGen4c1l)
PCA2ESepi <- projectR(data = p.ESepiGen4c1l$mRNA.Seq,loadings=pc.RNAseq6l3c3t,
                      full=TRUE, dataNames=map.ESepiGen4c1l[["GeneSymbols"]])

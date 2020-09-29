#run the library and anchor_drugs first, then run everything inside the calculate IC50 function
#can just run load("ic_results.Rdata") without calculating IC everytime
#Issues:
#calculating IC25, currently the IC25 is bigger than IC50 but it should be the other way around
setwd("~/Documents/Yeh Lab/HTS_Raw_Data/Scripts")
library(YehLabHTS)
library("data.table")
library(drexplorer)
library(reshape2)
library(ggplot2)
library("tidyverse")
library(anchors)
library(tidyverse)
library(hrbrthemes)
library(dbplyr)
library(scales)
calculateIC50 <- function(drug_results) {
  ic_results <- data.frame()
  #every library and anchor drug combination
  #pairs <- unique(paste(drug_results[[h]]$Compound, drug_results[[h]]$Anchor))

  anchor_drugs <- ()
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

######plotting and visualizing IC50 ranges

#replace the Inf values to NA
ic_results <- replace.value(ic_results, c("ic50","ic25"), from = Inf, to= NA, verbose = FALSE)

save(ic_results, file = "ic_results.Rdata")
load("/Users/madisonjenner/YehLabHTS/ic_results.Rdata")
load("ic_results.Rdata")

#removes ic50 lines that contain NA
ic50Data <- ic_results[!is.na(as.numeric(as.character(ic_results$ic50))),]
ic50df <- ic_results[!is.na(as.data.frame(ic_results$ic50)),]
ic50df <- as.data.frame(transform(ic50df, ic50 = as.numeric(ic50)))

#converts the ic50 column into numeric
ic50Data <- transform(ic50Data, ic50 = as.numeric(ic50))
ic50Data <- ic50Data[order(ic50Data$ic50),]

#get rid of duplicate ic50 values (ex lib + both anchors on a plate)
rmduplicate <- ic50Data %>% distinct(ic50, .keep_all = T)

ic50grouped <- rmduplicate %>%
  group_by(Library.Drug) %>%
  filter(ic50 > 0 & ic50 <= 1.6e-5)
ic50grouped <- head(ic50grouped, 100)
ggplot(ic50grouped, aes( x=Library.Drug, y=log10(ic50), fill = Library.Drug))+#group=Library.Drug)) +
  geom_boxplot()+ geom_jitter()+
  theme_bw()+
  labs(title = "Distribution of IC50 values within each compound", x= "Library Drugs", y = "Log10 Concentration")+
  theme(axis.text.x = element_text(angle = 90, size=12, hjust = 1), text = element_text(size =14), plot.title = element_text(hjust = 0.5),
                axis.text.y = element_text(size=12), axis.line = element_line(size = 0.8, color = "black"), panel.border = element_blank(),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#ylim(0, 3e-5)

ic50DataUnique <- ic50Data
ic50DataUnique <- ic50DataUnique[!duplicated(ic50DataUnique$Library.Drug), ]

#total lib/anchor unique combinations for 15 files = 5280
# 2586
#histogram of all combined files and their calculatable IC50s (2586)
hist(log(ic50Data$ic50),
     main = "Histogram of IC50 values from 2018 drug screen",
     xlab = "Log Concentration (M)",
     freq = T,
     xlim = c(-20, 0),
     ylim = c(0, 1200))
     #breaks = c(4e-9, 100e-9,300e-9, 1e-6, 3e-6, 10e-6, 0.0009))
# xlim = c(4.6e-9, 2e-5), 10e-9, 100e-9, 300e-9, 1e-6, 3e-6, 10e-6
abline(v = c(-8, -7, -6.5, -6, -5.5, -5), col = "red", lwd = 2)
abline(v = median(ic50Data$ic50), col = "red", lwd = 2)

hist.data = hist(ic50Data$ic50, plot=F)
hist.data$counts[hist.data$counts>0] <- log(hist.data$counts[hist.data$counts>0], 10)
plot(hist.data)

#max: 2e-4, min: 4.6e-9
res <- with(ic50DataUnique, hist(ic50[ic50 > 0 & ic50 < 2e-5],
   main = "IC50 values of KCGS drugs",
   xlab = "IC50 (M)",
   freq = T,
   #labels = T,
   #counts = T,
   breaks = c(0,10e-9, 100e-9, 300e-9, 1e-6, 3e-6, 10e-6, 2e-5),
   #xlim = range(breaks),
   ylim = c(0, 600)))
abline(v = c(10e-9, 100e-9, 300e-9, 1e-6), col = "red", lwd = 2)
abline(v = c(1e-6, 3e-6, 10e-6), col = "red", lwd = 2)
abline(v = c(10e-9, 100e-9, 300e-9, 1e-6, 3e-6, 10e-6), col = "red", lwd = 2)
#gives the frequency counts for each bin
res$counts
bins <- c(10e-9, 100e-9, 300e-9, 1e-6, 3e-6, 10e-6)
uniquehist<-with(ic50DataUnique, hist(ic50[ic50 > 0 & ic50 <= 1e-6],
          main = "IC50 Values of Unique Drugs",
          xlab = "IC50 (M)",
          freq = T,
          right = T, #include.lowest = T,
          #xaxt="n",
          breaks = 20,
          #breaks = bins, #seq(0, 1.5e-5, 1e-6),
          xlim = c(4e-9, 1e-6),
          ylim = c(0, 80),
          labels = T))
          #counts = F))
axis(1, at=seq(0, 2e-6, 1e-7))

#max 1.5e-5 (-4.7) min 4.6e-9
uniquehist <-ic50DataUnique %>%
  filter(ic50 > 0 & ic50 <= 1.6e-5) %>%
  ggplot(aes(x = log10(ic50) )) +
  labs(title = "IC50 values of 134 KCGS compounds", x= "Log Concentration (M)", y = "Density") +
  geom_histogram(aes(y =..density..), color = "black", fill = "gray", breaks = c(-8.3, -8, -7, -6.5, -6, -5.5, -5, -4.8)) +
  geom_density(alpha=.2, fill="#FF6666") +
  scale_y_continuous(breaks=seq(0, 1, .1))+
  theme_bw()+
  theme(text = element_text(size =14), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12), axis.line = element_line(size = 0.8, color = "black"), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
uniquehist

ic50DataUnique %>%
  filter(ic50 > 0 & ic50 <= 1.6e-5) %>%
  ggplot(aes(x = log10(ic50) )) +
  labs(title = "IC50 values of 134 KCGS compounds", x= "Log Concentration (M)", y= "Counts") +
  geom_histogram(aes(), color = "black", fill = "gray", breaks = c(-8.3, -8, -7, -6.5, -6, -5.5, -5, -4.8)) +
  #geom_density(alpha=.2, fill="#FF6666") +
  #scale_y_continuous(breaks=seq(0, 1, .1))+
  theme_bw()+
  theme(text = element_text(size =14), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12), axis.line = element_line(size = 0.8, color = "black"), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#density plot and hist
alldrugs <- rmduplicate %>%
  filter(ic50 > 0 & ic50 <= 2e-4) %>%
  ggplot(aes(x = log10(ic50) )) +
    labs(title = "IC50 values of KCGS compounds", x= "Log Concentration (M)", y = "Density") +
    geom_histogram(aes(y =..density..), color = "black", fill = "gray")+ #, breaks = c(-8.3, -8, -7, -6.5, -6, -5.5, -5, -3.6)) +
    geom_density(alpha=.2, fill="#FF6666") +
    scale_y_continuous(breaks=seq(0, 1, .1)) + #,  limits = c(0, 1), expand = c(0,-8))+
    expand_limits(x= -8.2, y=0)+
    theme_bw()+
    theme(text = element_text(size =16), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14), axis.line = element_line(size = 0.8, color = "black"), panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
alldrugs
#geom_vline(xintercept = c(-8, -7, -6.5, -6, -5.5, -5))

#just hist
rmduplicate %>%
  filter(ic50 > 0 & ic50 <= 2e-4) %>%
  ggplot(aes(x = log10(ic50) )) +
  labs(title = "IC50 values of KCGS compounds", x = "Log Concentration (M)", y= "Counts") +
  geom_histogram(aes(), color = "black", fill = "gray")+ #breaks = c(-8.3, -8, -7, -6.5, -6, -5.5, -5, -4.8))+
  #geom_density(alpha=.2, fill="#FF6666") +
  #scale_x_discrete(limits = c(-8, -7, -6.5, -6, -5.5, -5))+
  scale_y_continuous(breaks= seq(0, 250, 50), limits = c(0, 250), expand = c(0,-8.2))+ #,  +
  expand_limits(x= -8.2, y=0)+
  theme_bw()+
  theme(text = element_text(size =16, face="bold"), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14), axis.line = element_line(colour = "black"), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())



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

hist(ic25Data, freq = T, xlim = c(-1e-9, 5e-7), breaks = 100)
abline(v = median(ic50Data), col = "red", lwd = 2)

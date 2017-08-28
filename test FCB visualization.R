setwd("C:/Users/giudicev/Documents/Fluorescent Cell Barcoding/Analysis by R studio/For R analysis/BCD+stain")
rm(list=ls())
library(flowCore)
library(flowClust)
library(flowViz)
library(flowWorkspace)
library(ggcyto)
library(flowType)

filenames <- c("BCD+stain_HC3_003")
comp.matrix <- c("Compensation-1")

###Reading and compensation
BCD_stain <- read.FCS(paste(filenames[1], ".fcs", sep=""), transformation=FALSE, alter.names = T)
summary(BCD_stain)
comp.mat <- read.csv(file=paste(comp.matrix[1], ".csv", sep=""), row.names=1)
comp.mat <- as.matrix(comp.mat)
comp.mat
file_comps <- compensate(BCD_stain, comp.mat)
###Asinh transformation
tf <- transformList(from=colnames(file_comps)[4:10], tfun=asinh)
BCD_stain_trans<-tf %on% file_comps
###Filtering data using FCB dye channels
BCD_gate <- rectangleGate(filterId = "FluoRegion1", "BUV.396.A"=c(4,10), "APC.Cy7.A"=c(5,11))
BCD_filter = Subset(BCD_stain_trans, filter(BCD_stain_trans,BCD_gate))
###Clustering
FCB_clust <- flowClust(BCD_filter, varNames = c("BUV.396.A", "APC.Cy7.A"), K=9, B=50)
###Plotting the 9 barcoded clusters
par(mar=c(6,6,4,1)+.1)
plot1 <- plot(FCB_clust, data=BCD_filter, level=0.8, z.cutoff=0, xlim=c(4,10), ylim=c(5,11), 
     pch=16, cex=0.3, cex.axis=2, xlab="DyLight 350", ylab="DyLight 800", cex.lab=2, las=1)
axis(1, at=plot1, labels = FALSE, lwd = 3)
axis(2, at=plot1, labels = FALSE, lwd = 3)
axis(3, at=plot1, labels = FALSE, lwd = 3, tck=0)
axis(4, at=plot1, labels = FALSE, lwd = 3, tck =0)
###Extract clusters as separate objects
FCB_pop <- split(BCD_filter, FCB_clust, population=list(p1=1,p2=2,p3=3,p4=4,p5=5,p6=6,p7=7,p8=8,p9=9))
pop1 <- FCB_pop$p1
pop2 <- FCB_pop$p2
pop3 <- FCB_pop$p3
pop4 <- FCB_pop$p4
pop5 <- FCB_pop$p5
pop6 <- FCB_pop$p6
pop7 <- FCB_pop$p7
pop8 <- FCB_pop$p8
pop9 <- FCB_pop$p9
FCB_fs <- c(pop1, pop2, pop3, pop4, pop5, pop6, pop7, pop8, pop9) ###Combine object
###By using first population to get a prior, filtering and clustering for CD4 population
CD4gate <- rectangleGate(filterId = "FluoRegion2", "BV605.A"=c(4,10), "APC.A"=c(0,8))
CD4_filter = Subset(pop1, filter(pop1, CD4gate))
CD4_clust <- flowClust(CD4_filter, varNames = c("BV605.A", "APC.A"), K=3, B=100)
###Plotting CD4 vs CD3
par(mar=c(6,6,4,1)+.1)
plot(CD4_clust, data=CD4_filter, level=0.8, z.cutoff=0, xlim=c(4,10), ylim=c(2,8), 
     pch=16, cex=0.9, cex.axis=2, xlab="CD3-BV605", ylab="CD4-APC", cex.lab=2, las=1)
axis(1, at=plot1, labels = FALSE, lwd = 3)
axis(2, at=plot1, labels = FALSE, lwd = 3)
axis(3, at=plot1, labels = FALSE, lwd = 3, tck=0)
axis(4, at=plot1, labels = FALSE, lwd = 3, tck =0)
###Filter CD3+CD4+ cells for Vb usage
CD4_cells <- split(CD4_filter, CD4_clust, population=list(CD4=c(2)))
CD4_pop <- CD4_cells$CD4
Vbgate <- rectangleGate(filterId = "FluoRegion3", "FITC.A"=c(4,6), "PE.A"=c(4,6))
Vb_filter_CD4 = Subset(CD4_pop, filter(CD4_pop, Vbgate))
###Cluster Vb populations
Vb_clust_CD4 <- flowClust(Vb_filter_CD4, varNames = c("FITC.A", "PE.A"), K=4, B=100)
###Plot Vb populations
par(mar=c(6,6,4,1)+.1)
plot(Vb_clust_CD4, data=Vb_filter_CD4, level=0.8, z.cutoff=0, 
     pch=19, cex=0.9, cex.axis=2, las=1, cex.lab=1.5, ylab="")
axis(1, at=plot1, labels = FALSE, lwd = 3)
axis(2, at=plot1, labels = FALSE, lwd = 3)
axis(3, at=plot1, labels = FALSE, lwd = 3, tck=0)
axis(4, at=plot1, labels = FALSE, lwd = 3, tck =0)
###Get the prior for the analysis of all 9 populations
prior1 <- flowClust2Prior(Vb_clust_CD4, kappa=1)

###Filter Barcoded population for CD3 and CD8 expression
CD8gate <- rectangleGate(filterId = "FluoRegion4", "BV605.A"=c(4,10), "PE.Cy5.A"=c(0,10))
CD8_filter = Subset(pop1, filter(pop1, CD8gate))
###Cluster in order to identify CD3+CD8+ cells
CD8_clust <- flowClust(CD8_filter, varNames = c("BV605.A", "PE.Cy5.A"), K=3, B=100)
###Plot data
par(mar=c(6,6,4,1)+.1)
plot(CD8_clust, data=CD8_filter, level=0.8, z.cutoff=0, xlim=c(4,10), ylim=c(2,10), 
     pch=16, cex=0.9, cex.axis=2, xlab="CD3-BV605", ylab="CD8-PE-Cy5", cex.lab=2, las=1)
axis(1, at=plot1, labels = FALSE, lwd = 3)
axis(2, at=plot1, labels = FALSE, lwd = 3)
axis(3, at=plot1, labels = FALSE, lwd = 3, tck=0)
axis(4, at=plot1, labels = FALSE, lwd = 3, tck =0)
###Extract CD3+CD8+ cells for Vb usage
CD8_cells <- split(CD8_filter, CD8_clust, population=list(CD8=c(3)))
CD8_pop <- CD8_cells$CD8
Vbgate <- rectangleGate(filterId = "FluoRegion5", "FITC.A"=c(4,6), "PE.A"=c(4,6))
Vb_filter_CD8 = Subset(CD8_pop, filter(CD8_pop, Vbgate))
###Cluster Vb populations
Vb_clust_CD8 <- flowClust(Vb_filter_CD8, varNames = c("FITC.A", "PE.A"), K=4, B=100)
###Plot Vb populations
par(mar=c(6,6,4,1)+.1)
plot(Vb_clust_CD8, data=Vb_filter_CD8, level=0.8, z.cutoff=0, 
     pch=19, cex=0.9, cex.axis=2, las=1, ylim=c(4,6), ylab="")
axis(1, at=plot1, labels = FALSE, lwd = 3)
axis(2, at=plot1, labels = FALSE, lwd = 3)
axis(3, at=plot1, labels = FALSE, lwd = 3, tck=0)
axis(4, at=plot1, labels = FALSE, lwd = 3, tck =0)
###Get the prior for the analysis of all 9 population
prior2 <- flowClust2Prior(Vb_clust_CD8, kappa=1)

#Repeat the analysis for all 9 barcoded populations
  for (ii in 1:9){
    ###CD3+CD4+ cells
    CD4gate <- rectangleGate(filterId = "FluoRegion2", "BV605.A"=c(4,10), "APC.A"=c(0,8))
    CD4_filter = Subset(FCB_pop[[ii]], filter(FCB_pop[[ii]], CD4gate))
    CD4_clust <- flowClust(CD4_filter, varNames = c("BV605.A", "APC.A"), K=3, B=100)
    PropMarkers <- 9:10 ###Phenotypic analysis using CD3 and CD4
    MFIMarkers <- PropMarkers
    MarkerNames <- c('FSC.A','FSC.H','SSC.A','CD8','DyL350','DyL800','Vb.FITC','Vb.PE', 'CD4', 'CD3') ###Define markers names
    Res <- flowType(CD4_filter, PropMarkers, MFIMarkers, 'kmeans', MarkerNames);
    MFIs=Res@MFIs; ###Calculate MFI
    Proportions=Res@CellFreqs; ###Calculate proportions of each population
    Proportions <- Proportions / max(Proportions) ###Calculate proportions on total
    names(Proportions) <- unlist(lapply(Res@PhenoCodes,
                                        function(x){return(decodePhenotype(
                                          x,Res@MarkerNames[PropMarkers],
                                          Res@PartitionsPerMarker))})) ###Match populations and markers names
    ###Draw graph bar with results
    index=order(Proportions,decreasing=TRUE)[2:9]
    par(mar=c(15,6,4,1)+.1)
    par(lwd=2)
    bp=barplot(Proportions[index], axes=FALSE, names.arg=FALSE, col = "navajowhite1", ylim = c(0,1))
    axis(2, line=0.0001, cex.axis=1.7, las=1, lwd = 2);
    axis(1, at=bp+0.05, labels=names(Proportions[index]), 
         las=2, side=1, line=0.3, cex.axis=1.7, lwd = 2)
    ###Vb usage on CD3+CD4+ cells
    CD4_cells <- split(CD4_filter, CD4_clust, population=list(CD4=2))
    CD4_pop <- CD4_cells$CD4
    Vbgate <- rectangleGate(filterId = "FluoRegion3", "FITC.A"=c(4,6), "PE.A"=c(4,6))
    Vb_filter_CD4 = Subset(CD4_pop, filter(CD4_pop, Vbgate))
    PropMarkers <- 7:8 ###Phenotypic analysis using Vb-FITC Vb-PE
    MFIMarkers <- PropMarkers 
    MarkerNames <- c('FSC.A','FSC.H','SSC.A','CD8','DyL350','DyL800','Vb.FITC','Vb.PE', 'CD4', 'CD3') ###Define markers names
    Res <- flowType(Vb_filter_CD4, PropMarkers, MFIMarkers, 'kmeans', MarkerNames);
    MFIs=Res@MFIs; ###Calculate MFI
    Proportions=Res@CellFreqs; ###Calculate proportions of each population
    Proportions <- Proportions / max(Proportions) ###Calculate proportions on total
    names(Proportions) <- unlist(lapply(Res@PhenoCodes,
                                        function(x){return(decodePhenotype(
                                          x,Res@MarkerNames[PropMarkers],
                                          Res@PartitionsPerMarker))})) ###Match populations and markers names
    ###Draw graph bar with results
    index=order(Proportions,decreasing=TRUE)[2:9]
    par(mar=c(15,6,4,1)+.1)
    par(lwd=2)
    bp=barplot(Proportions[index], axes=FALSE, names.arg=FALSE, col = "navajowhite1", ylim = c(0,1))
    axis(2, line=0.0001, cex.axis=1.7, las=1, lwd=2);
    axis(1, at=bp+0.05, labels=names(Proportions[index]), 
         las=2, side=1, line=0.3, cex.axis=1.7,lwd=2)
    
    ###CD3+CD8+ cells
    CD8gate <- rectangleGate(filterId = "FluoRegion4", "BV605.A"=c(0,10), "PE.Cy5.A"=c(-4,12))
    CD8_filter = Subset(FCB_pop[[ii]], filter(FCB_pop[[ii]], CD8gate))
    CD8_clust <- flowClust(CD8_filter, varNames = c("BV605.A", "PE.Cy5.A"), K=3, B=100)
    PropMarkers <- c(4,10) ###Phenotypic analysis using CD3 and CD8
    MFIMarkers <- PropMarkers
    MarkerNames <- c('FSC.A','FSC.H','SSC.A','CD8','DyL350','DyL800','Vb.FITC','Vb.PE', 'CD4', 'CD3')
    Res <- flowType(CD8_filter, PropMarkers, MFIMarkers, 'kmeans', MarkerNames);
    MFIs=Res@MFIs; ###Calculate MFI
    Proportions=Res@CellFreqs; ###Calculate proportions of each population
    Proportions <- Proportions / max(Proportions) ###Calculate proportions on total
    names(Proportions) <- unlist(lapply(Res@PhenoCodes,
                                        function(x){return(decodePhenotype(
                                          x,Res@MarkerNames[PropMarkers],
                                          Res@PartitionsPerMarker))})) ###Match populations and markers names
    ###Draw graph bar with results
    index=order(Proportions,decreasing=TRUE)[2:9]
    par(mar=c(15,6,4,1)+.1)
    par(lwd=2)
    bp=barplot(Proportions[index], axes=FALSE, names.arg=FALSE, col = "navajowhite1", ylim = c(0,1))
    axis(2, line=0.0001, cex.axis=1.7, las=1, lwd=2);
    axis(1, at=bp+0.05, labels=names(Proportions[index]), 
         las=2, side=1, line=0.3, cex.axis=1.7, lwd=2)
    ###Vb usage on CD3+CD4+ cells
    CD8_cells <- split(CD8_filter, CD8_clust, population=list(CD8=3))
    CD8_pop <- CD8_cells$CD8
    Vbgate <- rectangleGate(filterId = "FluoRegion3", "FITC.A"=c(4,8), "PE.A"=c(4,8))
    Vb_filter_CD8 = Subset(CD8_pop, filter(CD8_pop, Vbgate))
    PropMarkers <- 7:8 ###Phenotypic analysis using Vb-FITC Vb-PE
    MFIMarkers <- PropMarkers
    MarkerNames <- c('FSC.A','FSC.H','SSC.A','CD8','DyL350','DyL800','Vb.FITC','Vb.PE', 'CD4', 'CD3')
    Res <- flowType(Vb_filter_CD8, PropMarkers, MFIMarkers, 'kmeans', MarkerNames);
    MFIs=Res@MFIs; ###Calculate MFI
    Proportions=Res@CellFreqs; ###Calculate proportions of each population
    Proportions <- Proportions / max(Proportions) ###Calculate proportions on total
    names(Proportions) <- unlist(lapply(Res@PhenoCodes,
                                        function(x){return(decodePhenotype(
                                          x,Res@MarkerNames[PropMarkers],
                                          Res@PartitionsPerMarker))})) ###Match populations and markers names
    ###Draw graph bar with results
    index=order(Proportions,decreasing=TRUE)[2:9]
    par(mar=c(15,6,4,1)+.1)
    par(lwd=2)
    bp=barplot(Proportions[index], axes=FALSE, names.arg=FALSE, col = "navajowhite1", ylim = c(0,1))
    axis(2, line=0.0001, cex.axis=1.7, las=1, lwd=2);
    axis(1, at=bp+0.05, labels=names(Proportions[index]), 
         las=2, side=1, line=0.3, cex.axis=1.7, lwd=2)
}

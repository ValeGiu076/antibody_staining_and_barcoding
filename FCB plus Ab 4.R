setwd("C:/Users/giudicev/Documents/Fluorescent Cell Barcoding/Analysis by R studio/For R analysis/BCD+stain")
rm(list=ls())
library(flowCore)
library(flowClust)
library(flowViz)
library(flowWorkspace)
library(ggcyto)
library(flowType)

filenames <- c("BCD+stain_HC3_003")
comp.matrix <- c("Compensation-3")

###read the data
BCD_stain <- read.FCS(paste(filenames[1], ".fcs", sep=""), transformation=FALSE, alter.names = T)
summary(BCD_stain)
####Let us remove debris
nodebris <- rectangleGate(filterId="no Debris", "FSC.A"=c(20000,Inf), "SSC.A"=c(0,Inf))
BCD_stain_nodebris <- Subset(BCD_stain, filter(BCD_stain,nodebris))
comp.mat <- read.csv(file=paste(comp.matrix[1], ".csv", sep=""), row.names=1)
comp.mat <- as.matrix(comp.mat)
comp.mat
file_comps <- compensate(BCD_stain_nodebris, comp.mat)
###Let us transform the data
tf <- transformList(from=colnames(file_comps)[4:10], tfun=asinh)
BCD_stain_trans<-tf %on% file_comps

######Let us cluster for FCB dyes
BCD_gate <- rectangleGate(filterId = "FluoRegion1", "BUV.396.A"=c(0,12), "APC.Cy7.A"=c(0,12))
BCD_filter = Subset(BCD_stain_trans, filter(BCD_stain_trans,BCD_gate))
FCB_clust <- flowClust(BCD_filter, varNames = c("BUV.396.A", "APC.Cy7.A"), K=9, B=100)

####For cluster assignment
#head(exprs(BCD_stain_trans_filter))
result<-cbind(exprs(BCD_filter)[, c("BUV.396.A", "APC.Cy7.A", "APC.A", "PE.Cy5.A", "FITC.A", "PE.A")], FCB_clust@z,
              matrix(rep("unknown", length(FCB_clust@label)), ncol=1))
colnames(result)<-c("BUV.396.A", "APC.Cy7.A", "APC.A", "PE.Cy5.A", "FITC.A", "PE.A", paste("cluster", 1:9, sep=""), "clusterAssigned")
for(jj in 1:length(FCB_clust@label)){
  if(!is.na(result[jj, "cluster1"])){result[jj, "clusterAssigned"]<-as.numeric(which.max(result[jj, paste("cluster", 1:9, sep="")]))}
}
write.csv(result, file=paste(filenames[1], "_withPrior_result.csv", sep=""), row.names=FALSE)

png(width=2000, height=2000, res=200, file=paste(filenames[1], "_clust.png", sep=""))
par(mar=c(6,6,4,1)+.1)
plot1 <- plot(FCB_clust, data=BCD_filter, level=0.8, z.cutoff=0, xlim=c(0,12), ylim=c(0,12), 
              pch=16, cex=0.3, cex.axis=2, xlab="DyLight 350", ylab="DyLight 800", cex.lab=2, las=1)
axis(1, at=plot1, labels = FALSE, lwd = 3)
axis(2, at=plot1, labels = FALSE, lwd = 3)
axis(3, at=plot1, labels = FALSE, lwd = 3, tck=0)
axis(4, at=plot1, labels = FALSE, lwd = 3, tck =0)
plot1
dev.off()
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
fs <- as(FCB_fs, "flowSet") ###Read data as flowset 
data1 <- summary(fs)
write.csv(data1, file =paste(filenames[1],"_for_gate.csv", sep = ""))

#####CD3+ cells identification
CD3gate <- rectangleGate("BV605.A"=c(5.5,10), "SSC.A"=c(0,1500))
BCD.CD3.filter <- filter(fs, CD3gate)
CD3_pop <- Subset(fs, BCD.CD3.filter)
png(width=2000, height=2000, res=200, file=paste(filenames[1], "_CD3.png", sep=""))
par(mar=c(6,6,4,1)+.1)
plot6 <- ggcyto(fs, aes(x="BV605.A", y="SSC.A")) + geom_hex(bins=128)
plot6 <- plot6 + ggtitle("CD3+ T cells") + xlim(c(0,10)) + ylim(c(0,4e3)) + geom_gate(CD3gate)
plot6
dev.off()

###CD4+ cells identification
CD4_gate <- rectangleGate("APC.A"=c(1.1,6.9), "SSC.A"=c(232,1500))
BCD.CD4.filter <- filter(CD3_pop, CD4_gate)
CD4pos <- Subset(CD3_pop, filter(CD3_pop, CD4_gate))
png(width=2000, height=2000, res=200, file=paste(filenames[1], "_CD4.png", sep=""))
plot4 <- ggcyto(CD3_pop, aes(x="APC.A", y="SSC.A")) + geom_hex(bins=128)
plot4 <- plot4 + ggtitle("CD4+ T cells") + xlim(c(-10,10)) + ylim(c(0,4000)) + geom_gate(CD4_gate)
plot4
dev.off()

###CD8+ cells identification
CD8_gate <- rectangleGate("PE.Cy5.A"=c(5.5,8.3), "SSC.A"=c(232,1500))
BCD.CD8.filter <- filter(CD3_pop, CD8_gate)
CD8pos <- Subset(CD3_pop, filter(CD3_pop, CD4_gate))
png(width=2000, height=2000, res=200, file=paste(filenames[1], "_CD8.png", sep=""))
plot5 <- ggcyto(CD3_pop, aes(x="PE.Cy5.A", y="SSC.A")) + geom_hex(bins=128)
plot5 <- plot5 + ggtitle("CD8+ T cells") + xlim(c(0,10)) + ylim(c(0,4000)) + geom_gate(CD8_gate)
plot5
dev.off()

###Vb usage gates
Vb9 <- rectangleGate("PE.A"=c(5,7.8), "FITC.A"=c(4.5,5.2))
Vb17 <- rectangleGate("PE.A"=c(5,7.8), "FITC.A"=c(5.2,7.7))
Vb16 <- rectangleGate("PE.A"=c(3.2,5), "FITC.A"=c(5.2,7.7))

##Vb usage on CD4+ T cells
CD4_Vb9 <- filter(CD4pos,Vb9)
CD4_Vb17 <- filter(CD4pos, Vb17)
CD4_Vb16 <- filter(CD4pos, Vb16)
png(width=2000, height=2000, res=200, file=paste(filenames[1], "_CD4_Vb.png", sep=""))
plot2 <- ggcyto(CD4pos, aes(x="FITC.A", y="PE.A")) + geom_hex(bins=128)
plot2 <- plot2 + ggtitle("Vb usage on CD4+ T cells") + xlim(c(0,10)) + ylim(c(0,10))
plot2 <- plot2 + geom_gate(Vb9) + geom_gate(Vb17) + geom_gate(Vb16)
plot2
dev.off()

##Vb usage on CD8+ T cells
CD8_Vb9 <- filter(CD8pos,Vb9)
CD8_Vb17 <- filter(CD8pos, Vb17)
CD8_Vb16 <- filter(CD8pos, Vb16)
png(width=2000, height=2000, res=200, file=paste(filenames[1], "_CD8_Vb.png", sep=""))
plot3 <- ggcyto(CD8pos, aes(x="FITC.A", y="PE.A")) + geom_hex(bins=128)
plot3 <- plot3 + ggtitle("Vb usage on CD8+ T cells") + xlim(c(0,10)) + ylim(c(0,10))
plot3 <- plot3 + geom_gate(Vb9) + geom_gate(Vb17) + geom_gate(Vb16)
plot3
dev.off()

###Let us write the results
CD3data <- summary(BCD.CD3.filter)
V1 <- CD3data$V1@p
V2 <- CD3data$V2@p
V3 <- CD3data$V3@p
V4 <- CD3data$V4@p
V5 <- CD3data$V5@p
V6 <- CD3data$V6@p
V7 <- CD3data$V7@p
V8 <- CD3data$V8@p
V9 <- CD3data$V9@p
dataCD3 <- cbind(V1, V2, V3, V4, V5, V6, V7, V8, V9)

resultsCD4 <- lapply(head(BCD.CD4.filter, 9), summary)
V1 <- resultsCD4$V1@p
V2 <- resultsCD4$V2@p
V3 <- resultsCD4$V3@p
V4 <- resultsCD4$V4@p
V5 <- resultsCD4$V5@p
V6 <- resultsCD4$V6@p
V7 <- resultsCD4$V7@p
V8 <- resultsCD4$V8@p
V9 <- resultsCD4$V9@p
dataCD4 <- cbind(V1, V2, V3, V4, V5, V6, V7, V8, V9)

resultsCD8 <- lapply(head(BCD.CD8.filter, 9), summary)
V1 <- resultsCD8$V1@p
V2 <- resultsCD8$V2@p
V3 <- resultsCD8$V3@p
V4 <- resultsCD8$V4@p
V5 <- resultsCD8$V5@p
V6 <- resultsCD8$V6@p
V7 <- resultsCD8$V7@p
V8 <- resultsCD8$V8@p
V9 <- resultsCD8$V9@p
dataCD8 <- cbind(V1, V2, V3, V4, V5, V6, V7, V8, V9)

results_CD4_Vb9 <- summary(CD4_Vb9)
V1 <- results_CD4_Vb9$V1@p
V2 <- results_CD4_Vb9$V2@p
V3 <- results_CD4_Vb9$V3@p
V4 <- results_CD4_Vb9$V4@p
V5 <- results_CD4_Vb9$V5@p
V6 <- results_CD4_Vb9$V6@p
V7 <- results_CD4_Vb9$V7@p
V8 <- results_CD4_Vb9$V8@p
V9 <- results_CD4_Vb9$V9@p
dataVb9_on_CD4 <- cbind(V1, V2, V3, V4, V5, V6, V7, V8, V9)

results_CD4_Vb17 <- summary(CD4_Vb17)
V1 <- results_CD4_Vb17$V1@p
V2 <- results_CD4_Vb17$V2@p
V3 <- results_CD4_Vb17$V3@p
V4 <- results_CD4_Vb17$V4@p
V5 <- results_CD4_Vb17$V5@p
V6 <- results_CD4_Vb17$V6@p
V7 <- results_CD4_Vb17$V7@p
V8 <- results_CD4_Vb17$V8@p
V9 <- results_CD4_Vb17$V9@p
dataVb17_on_CD4 <- cbind(V1, V2, V3, V4, V5, V6, V7, V8, V9)

results_CD4_Vb16 <- summary(CD4_Vb16)
V1 <- results_CD4_Vb16$V1@p
V2 <- results_CD4_Vb16$V2@p
V3 <- results_CD4_Vb16$V3@p
V4 <- results_CD4_Vb16$V4@p
V5 <- results_CD4_Vb16$V5@p
V6 <- results_CD4_Vb16$V6@p
V7 <- results_CD4_Vb16$V7@p
V8 <- results_CD4_Vb16$V8@p
V9 <- results_CD4_Vb16$V9@p
dataVb16_on_CD4 <- cbind(V1, V2, V3, V4, V5, V6, V7, V8, V9)

results_CD8_Vb9 <- summary(CD8_Vb9)
V1 <- results_CD8_Vb9$V1@p
V2 <- results_CD8_Vb9$V2@p
V3 <- results_CD8_Vb9$V3@p
V4 <- results_CD8_Vb9$V4@p
V5 <- results_CD8_Vb9$V5@p
V6 <- results_CD8_Vb9$V6@p
V7 <- results_CD8_Vb9$V7@p
V8 <- results_CD8_Vb9$V8@p
V9 <- results_CD8_Vb9$V9@p
dataVb9_on_CD8 <- cbind(V1, V2, V3, V4, V5, V6, V7, V8, V9)

results_CD8_Vb17 <- summary(CD8_Vb17)
V1 <- results_CD8_Vb17$V1@p
V2 <- results_CD8_Vb17$V2@p
V3 <- results_CD8_Vb17$V3@p
V4 <- results_CD8_Vb17$V4@p
V5 <- results_CD8_Vb17$V5@p
V6 <- results_CD8_Vb17$V6@p
V7 <- results_CD8_Vb17$V7@p
V8 <- results_CD8_Vb17$V8@p
V9 <- results_CD8_Vb17$V9@p
dataVb17_on_CD8 <- cbind(V1, V2, V3, V4, V5, V6, V7, V8, V9)

results_CD8_Vb16 <- summary(CD8_Vb16)
V1 <- results_CD8_Vb16$V1@p
V2 <- results_CD8_Vb16$V2@p
V3 <- results_CD8_Vb16$V3@p
V4 <- results_CD8_Vb16$V4@p
V5 <- results_CD8_Vb16$V5@p
V6 <- results_CD8_Vb16$V6@p
V7 <- results_CD8_Vb16$V7@p
V8 <- results_CD8_Vb16$V8@p
V9 <- results_CD8_Vb16$V9@p
dataVb16_on_CD8 <- cbind(V1, V2, V3, V4, V5, V6, V7, V8, V9)

results_all <- rbind(dataCD3, dataCD4, dataCD8, dataVb9_on_CD4,dataVb17_on_CD4,dataVb16_on_CD4,dataVb9_on_CD8,dataVb17_on_CD8,dataVb16_on_CD8)
rownames(results_all) <- c("CD3+", "CD4+", "CD8+", "CD4+Vb9+", "CD4+Vb17+", "CD4+Vb16+",
                           "CD8+Vb9+", "CD8+Vb17+", "CD8+Vb16+")
write.csv(results_all, file =paste(filenames[1],"_results_Vb_usage.csv", sep = ""))


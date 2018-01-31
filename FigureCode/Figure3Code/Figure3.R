library(ggplot2)
library(cowplot)
library(gridExtra)
rm(list = ls())

VertCount <- c(105, 143)
VertMain <- c("Supermatrix", "Coalescent")
VertMyData <- data.frame(VertMain,VertCount)
VertMyData$VertMain <- factor(VertMyData$VertMain, levels = c("Supermatrix", "Coalescent"))
pVertMain <- ggplot(VertMyData, aes(VertMyData$VertMain,VertMyData$VertCount)) + ylim(0,248)
pVertMain = pVertMain + geom_bar(stat = "identity") + xlab("Topology") + ylab("Gene Count") + ggtitle("Vertebrate Two Topology Test")


VertCount2 <- c(91, 144, 13)
VertMain2 <- c("Supermatrix", "Coalescent", "Alternatives")
VertMyData2 <- data.frame(VertMain2,VertCount2)
VertMyData2$VertMain2 <- factor(VertMyData2$VertMain2, levels = c("Supermatrix", "Coalescent", "Alternatives"))
pVertMain2 <- ggplot(VertMyData2, aes(VertMyData2$VertMain2,VertMyData2$VertCount2)) + ylim(0,248)
pVertMain2 = pVertMain2 + geom_bar(stat = "identity") + xlab("Edge") + ylab("Gene Count") + ggtitle("Vertebrate MGWE")
pVertMain2

CarCount <- c(623, 614)
CarMain <- c("Supermatrix", "Coalescent")
CarMyData <- data.frame(CarMain,CarCount)
CarMyData$CarMain <- factor(CarMyData$CarMain, levels = c("Supermatrix", "Coalescent"))
pCarMain <- ggplot(CarMyData, aes(CarMyData$CarMain,CarMyData$CarCount)) + ylim(0,1237)
pCarMain = pCarMain + geom_bar(stat = "identity") + xlab("Topology") + ylab("Gene Count") + ggtitle("Carnivory Two Topology Test")
pCarMain

CarCount2 <- c(499, 466, 272)
CarMain2 <- c("Supermatrix", "Coalescent", "Alternatives")
CarMyData2 <- data.frame(CarMain2,CarCount2)
CarMyData2$CarMain2 <- factor(CarMyData2$CarMain2, levels = c("Supermatrix", "Coalescent", "Alternatives"))
pCarMain2 <- ggplot(CarMyData2, aes(CarMyData2$CarMain2,CarMyData2$CarCount2)) + ylim(0,1237)
pCarMain2 = pCarMain2 + geom_bar(stat = "identity") + xlab("Edge") + ylab("Gene Count") + ggtitle("Carnivory MGWE")
pCarMain2

plot_grid(pVertMain,pVertMain2,pCarMain,pCarMain2,labels=c("A","B","C","D"),ncol=2,nrow=2)

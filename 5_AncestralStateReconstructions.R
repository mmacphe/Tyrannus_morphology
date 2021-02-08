#Mapping ancestral traits
rm(list=ls())

require(phytools)

phy<-read.tree("Tyrannus_phylogeny.tre")
plotTree(phy)
is.binary(phy)
phy<-multi2di(phy) #resolve polytomies
is.binary(phy)

Tyrannus.data<-read.csv("Tyrannus morphology + PCA avg.csv", row.names = 1)
species<-rownames(Tyrannus.data)
Tyrannus.data<-cbind(species,Tyrannus.data)
View(Tyrannus.data) #to see that column no. for Bill PC2

### Check that tree$tip.label is the same as otu_avg
phy$tip.label[!phy$tip.label %in% Tyrannus.data$species] #should return: charactor(0)

#BillPC2
png(file="ASR_BillPPC2.png",width=6.5,height=5.5,units="in",res=500)
BillPC2<-Tyrannus.data$BillPC2
names(BillPC2)<- rownames(Tyrannus.data)
View(BillPC2)
obj<-contMap(phy, BillPC2)
plot(obj,legend=0.7*max(nodeHeights(phy)), fsize=c(0.7,0.9))
dev.off()

Tyrannus.data<-as.matrix(Tyrannus.data)[,12] #change to a vector; 12 is the Bill PC2 column
View(Tyrannus.data) #this does have rownames
class(Tyrannus.data) #returns: "character"
fit<-fastAnc(phy, as.numeric(Tyrannus.data), vars=TRUE, CI=TRUE)
fit
fit$CI[1,]
range(as.numeric(Tyrannus.data))
phenogram(phy, BillPC2, fsize=0.6, spread.costs=c(1,0))

#Bill Length/Tarsus
png(file="ASR_Bill_Length.png",width=6.5,height=5.5,units="in",res=500)
BillLength<-Tyrannus.data$BL.Average/Tyrannus.data$Tarsus.Average
names(BillLength)<- rownames(Tyrannus.data)
View(BillLength)
obj<-contMap(phy, BillLength, plot=FALSE)
plot(obj,legend=0.7*max(nodeHeights(phy)), fsize=c(0.7,0.9))
dev.off()

#Bill Width/Tarsus
png(file="ASR_Bill_Width.png",width=6.5,height=5.5,units="in",res=500)
BillWidth<-Tyrannus.data$BW.Average/Tyrannus.data$Tarsus.Average
names(BillWidth)<- rownames(Tyrannus.data)
View(BillWidth)
obj<-contMap(phy, BillWidth, plot=FALSE)
plot(obj,legend=0.7*max(nodeHeights(phy)), fsize=c(0.7,0.9))
dev.off()

#Bill Depth/Tarsus
png(file="ASR_Bill_Depth.png",width=6.5,height=5.5,units="in",res=500)
BillDepth<-Tyrannus.data$BD.Average/Tyrannus.data$Tarsus.Average
names(BillDepth)<- rownames(Tyrannus.data)
View(BillDepth)
obj<-contMap(phy, BillDepth, plot=FALSE)
plot(obj,legend=0.7*max(nodeHeights(phy)), fsize=c(0.7,0.9))
dev.off()

#Wing Cord/Tarsus
png(file="ASR_Wing_Cord.png",width=6.5,height=5.5,units="in",res=500)
WingCord<-Tyrannus.data$WC.Average/Tyrannus.data$Tarsus.Average
names(WingCord)<- rownames(Tyrannus.data)
View(WingCord)
obj<-contMap(phy, WingCord, plot=FALSE)
plot(obj,legend=0.7*max(nodeHeights(phy)), fsize=c(0.7,0.9))
dev.off()

#Kipp's Index/Tarsus
png(file="ASR_Kipps_Index.png",width=6.5,height=5.5,units="in",res=500)
KippsIndex<-Tyrannus.data$Kipp.s.Average/Tyrannus.data$Tarsus.Average
names(KippsIndex)<- rownames(Tyrannus.data)
View(KippsIndex)
obj<-contMap(phy, KippsIndex, plot=FALSE)
plot(obj,legend=0.7*max(nodeHeights(phy)), fsize=c(0.7,0.9))
dev.off()

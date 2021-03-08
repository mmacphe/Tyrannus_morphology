#Mapping ancestral traits
rm(list=ls())

require(phytools)

### Set source directory to the folder this file came from within RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

phy<-read.tree('./Output Files/Tyrannus_phylogeny.tre')
plotTree(phy)
is.binary(phy)
phy<-multi2di(phy) #resolve polytomies
is.binary(phy)

Tyrannus.data<-read.csv('./Output Files/Tyrannus morphology + PCA avg.csv', row.names = 1)
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

### Migration strategy
### Reload the dataset to save as a factor instead of a vector
Tyrannus.data<-read.csv('./Output Files/Tyrannus morphology + PCA avg.csv', row.names = 1)
species<-rownames(Tyrannus.data)
Tyrannus.data<-cbind(species,Tyrannus.data)

### Check that tree$tip.label is the same as otu_avg
phy$tip.label[!phy$tip.label %in% Tyrannus.data$species] #should return: charactor(0)


### Use diversitree, the mk-n model to reconstruct marginal ancestral states
require(diversitree)

#Change character to numeric, because this is the format diversitree accepts
char1<-as.numeric(Migration)
names(char1)<-names(Migration)

#Create mkn model. Use equal probabilites for all states at the root
lik.mn<-make.mkn(phy,char1,k=3, control=list(root=ROOT.EQUI))
#Look at the different parameters of the model
argnames(lik.mn)
#Create model using constraints. 1=migratory, 2=partial, 3=sedentary. q12 is the transition rate from 1 to 2, and ~ indicates equivalency. 
lik.mkn.base<-constrain(lik.mn, q21~q13, q23~q32, q12~q31)
#Create a starting point for the search
p.mkn<-starting.point.musse(phy,3)
#Fit the model
fit.mkn<-find.mle(lik.mkn.base,p.mkn[argnames(lik.mkn.base)])
#Model parameters
fit.mkn[1:2]
#Export marginal ancestral reconstruction at the notes of the tree
st<-t(asr.marginal(lik.mkn.base,coef(fit.mkn)))

### Plot the results with the probabilities of the different states at each noce
#Vector of colours
co<-c("black","gray40","gray70")
plot(phy,type="p",FALSE,label.offset=0.6,cex=0.6,no.margin=TRUE,edge.color="gray32", tip.color = "gray32")
tiplabels(pch=21,bg=co[as.numeric(char1),col="gray32", cex=1,adj=0.6])
nodelabels(pie=st, piecol=co, cex=0.5, col="gray32")
legend("bottomleft",legend=levels(Migration), pch=20, col=co, bty="n", text.col="gray32", cex=0.8, pt.cex=2)


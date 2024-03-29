rm(list=ls())

require(phytools)
require(RColorBrewer)
require(png)

### Set source directory to the folder this file came from within RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Read in average morphology data from .csv ###
morpho<-read.csv('./Output Files/Tyrannus_Adults_data.csv', row.names=2)

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#morpho<-read.csv('./Output Files/Tyrannus_Females_Adults_data.csv', row.names=2)
#OR
#morpho<-read.csv('./Output Files/Tyrannus_Males_Adults_data.csv', row.names=2)

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

### Read in phylogeny 
phy<-read.tree('./Output Files/Tyrannus_phylogeny.tre')

### Check that tree$tip.label is the same as morpho
phy$tip.label[!phy$tip.label %in% rownames(morpho)]

### Bill Phylogenetic PCA ###
bill_pca<-phyl.pca(phy,morpho[,3:5],method="lambda")
bill_pca

### Find the % variance explained by each eigenvector 
diag(bill_pca$Eval)/sum(bill_pca$Eval)*100

### Bring in morphology data set for all individuals to add points onto figure
morpho_whole<-read.csv('./Output Files/Tyrannus Adults morphology data.csv', row.names = 1)

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#morpho_whole<-read.csv('./Output Files/Tyrannus Adults Females morphology data.csv', row.names=1)
#OR
#morpho_whole<-read.csv('./Output Files/Tyrannus Adults Males morphology data.csv', row.names=1)

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

### Match tip.labels with phylogeny tip.labels ###
morpho_whole$tip.label
names(table(morpho_whole$tip.label))

phy$tip.label

phy$tip.label[!phy$tip.label %in% morpho_whole$tip.label] #should return character(0)
morpho_whole$tip.label[! morpho_whole$tip.label %in% phy$tip.label] #should return character(0)

### Split the morphology data sets up into categories by tip.label
morpho_split<-split(morpho_whole,morpho_whole$tip.label)
names(morpho_split)
phy$tip.label

phy$tip.label[!phy$tip.label %in% names(morpho_split)] #should return character(0)
names(morpho_split)[! names(morpho_split) %in% phy$tip.label] #should return character (0)

morpho_trim<-lapply(morpho_split,function(x) x[colnames(x) %in% c("BL.Average","BW.Average","BD.Average", "Kipp.s.Distance", "WC.Average", "Tail", "Tarsus.Average", "name_num")])

bill_mat<-matrix(as.numeric(as.matrix(do.call(rbind, lapply(morpho_trim,function(x) x[,1:3])))),nrow=nrow(do.call(rbind,morpho_trim))) 

bill_scores<-phytools:::scores(bill_pca,newdata=bill_mat)
rownames(bill_scores)<-rownames(do.call(rbind, lapply(morpho_trim,function(x) x[,1:3])))
bill_scores_df<-data.frame(bill_scores)
bill_scores_df$otu<-factor(sapply(strsplit(rownames(bill_scores_df),"[.]"),function(x) x[1]))
bill_scores_split<-split(bill_scores_df, bill_scores_df$otu)
bill_scores_avg<-t(sapply(bill_scores_split,function(x) c(mean(na.omit(x[,1])),mean(na.omit(x[,2])))))

bill_scores_avg<-as.data.frame(bill_scores_avg)
morpho$BillPC1<-bill_scores_avg$V1
morpho$BillPC2<-bill_scores_avg$V2

### Flight Feather Phylogenetic PCA ###
body_pca<-phyl.pca(phy,morpho[,7:9],method="lambda")
body_pca

### Find the % variance explained by each eigenvector 
diag(body_pca$Eval)/sum(body_pca$Eval)*100

body_mat<-matrix(as.numeric(as.matrix(do.call(rbind, lapply(morpho_trim,function(x) x[,5:7])))),nrow=nrow(do.call(rbind,morpho_trim))) 
body_scores<-phytools:::scores(body_pca,newdata=body_mat)
rownames(body_scores)<-rownames(do.call(rbind, lapply(morpho_trim,function(x) x[,5:7])))
body_scores_df<-data.frame(body_scores)
body_scores_df$otu<-factor(sapply(strsplit(rownames(body_scores_df),"[.]"),function(x) x[1]))
body_scores_split<-split(body_scores_df, body_scores_df$otu)
body_scores_avg<-t(sapply(body_scores_split,function(x) c(mean(na.omit(x[,1])),mean(na.omit(x[,2])))))

body_scores_avg<-as.data.frame(body_scores_avg)
morpho$BodyPC1<-body_scores_avg$V1

morpho$BodyPC2<-body_scores_avg$V2

### Body PCA without T. forficatus or T. savana bc their tails are so long ###
phy_trim<-drop.tip(phy,c(grep("forficatus",phy$tip.label),grep("savana",phy$tip.label)))
otu_avg_trim<-morpho[rownames(morpho) %in% phy_trim$tip.label,]
body_pca_trim<-phyl.pca(phy_trim,otu_avg_trim[,7:9],method="lambda")
body_pca_trim

### Find the % variance explained by each eigenvector 
diag(body_pca_trim$Eval)/sum(body_pca_trim$Eval)*100

body_mat_trim<-matrix(as.numeric(as.matrix(do.call(rbind, lapply(morpho_trim[-c(grep("forficatus",names(morpho_trim)),grep("savana",names(morpho_trim)))],function(x) x[,5:7])))),nrow=nrow(do.call(rbind,morpho_trim[-c(grep("forficatus",names(morpho_trim)),grep("savana",names(morpho_trim)))]))) 

body_scores_trim<-phytools:::scores(body_pca_trim,newdata=body_mat_trim)
rownames(body_scores_trim)<-rownames(do.call(rbind, lapply(morpho_trim[-c(grep("forficatus",names(morpho_trim)),grep("savana",names(morpho_trim)))],function(x) x[,5:7])))
body_scores_trim_df<-data.frame(body_scores_trim)
body_scores_trim_df$otu<-factor(sapply(strsplit(rownames(body_scores_trim_df),"[.]"),function(x) x[1]))
body_scores_split_trim<-split(body_scores_trim_df, body_scores_trim_df$otu)
body_scores_trim_avg<-t(sapply(body_scores_split_trim,function(x) c(mean(na.omit(x[,1])),mean(na.omit(x[,2])))))

### Write.csv of morpho that includes all OTUs ###
write.csv(morpho, file="Tyrannus morphology + PCA avg_Adults.csv")

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#write.csv(morpho, file="Tyrannus morphology + PCA avg_Females_Adults.csv")
#OR
#write.csv(morpho, file="Tyrannus morphology + PCA avg_Males_Adults.csv")

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

### Bring in cv summary dataset 
cv_summary<-read.csv('./Output Files/cv_summary_Adults_table.csv', row.names = 2)

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#cv_summary<-read.csv('./Output Files/cv_summary_Females_Adults_table.csv', row.names = 2)
#OR
#cv_summary<-read.csv('./Output Files/cv_summary_Males_Adults_table.csv', row.names = 2)

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

### Coefficient of Variation (CV) Bill Phylogenetic PCA ###
CVbill_pca<-phyl.pca(phy,cv_summary[,3:5],method="lambda")
CVbill_pca

### Find the % variance explained by each eigenvector 
diag(CVbill_pca$Eval)/sum(CVbill_pca$Eval)*100

CVbill_mat<-matrix(as.numeric(as.matrix(do.call(rbind, lapply(morpho_trim,function(x) x[,1:3])))),nrow=nrow(do.call(rbind,morpho_trim))) 
CVbill_scores<-phytools:::scores(CVbill_pca,newdata=CVbill_mat)
rownames(CVbill_scores)<-rownames(do.call(rbind, lapply(morpho_trim,function(x) x[,1:3])))
CVbill_scores_df<-data.frame(CVbill_scores)
CVbill_scores_df$otu<-factor(sapply(strsplit(rownames(CVbill_scores_df),"[.]"),function(x) x[1]))
CVbill_scores_split<-split(CVbill_scores_df, CVbill_scores_df$otu)
CVbill_scores_cv<-t(sapply(CVbill_scores_split,function(x) c(mean(na.omit(x[,1])),mean(na.omit(x[,2])))))

CVbill_scores_cv<-as.data.frame(CVbill_scores_cv)
cv_summary$CVBillPC1<-CVbill_scores_cv$V1
cv_summary$CVBillPC2<-CVbill_scores_cv$V2

### CV Flight Feather Phylogenetic PCA ###
CVbody_pca<-phyl.pca(phy,cv_summary[,7:9],method="lambda")

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing females, use one of the following lines of code instead of the above line, then proceed.
#CVbody_pca<-phyl.pca(phy,na.omit(cv_summary[,6:8]),method="lambda") #there is a missing value in the tail column

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

CVbody_pca

### Find the % variance explained by each eigenvector 
diag(CVbody_pca$Eval)/sum(CVbody_pca$Eval)*100

CVbody_mat<-matrix(as.numeric(as.matrix(do.call(rbind, lapply(morpho_trim,function(x) x[,5:7])))),nrow=nrow(do.call(rbind,morpho_trim))) 
CVbody_scores<-phytools:::scores(CVbody_pca,newdata=CVbody_mat)
rownames(CVbody_scores)<-rownames(do.call(rbind, lapply(morpho_trim,function(x) x[,5:7])))
CVbody_scores_df<-data.frame(CVbody_scores)
CVbody_scores_df$otu<-factor(sapply(strsplit(rownames(CVbody_scores_df),"[.]"),function(x) x[1]))
CVbody_scores_split<-split(CVbody_scores_df, CVbody_scores_df$otu)
CVbody_scores_cv<-t(sapply(CVbody_scores_split,function(x) c(mean(na.omit(x[,1])),mean(na.omit(x[,2])))))

CVbody_scores_cv<-as.data.frame(CVbody_scores_cv)
cv_summary$CVBodyPC1<-CVbody_scores_cv$V1
cv_summary$CVBodyPC2<-CVbody_scores_cv$V2

write.csv(cv_summary, file="cv_summary_Adults.csv")

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#write.csv(cv_summary, file="cv_summary_Females_Adults.csv")
#OR
#write.csv(cv_summary, file="cv_summary_Males_Adults.csv")

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###


### CV Body PCA without T. forficatus or T. savana bc their tails are so long ###
phy_trim<-drop.tip(phy,c(grep("forficatus",phy$tip.label),grep("savana",phy$tip.label)))
otu_cv_trim<-cv_summary[rownames(cv_summary) %in% phy_trim$tip.label,]
CVbody_pca_trim<-phyl.pca(phy_trim,otu_cv_trim[,7:9],method="lambda")
CVbody_pca_trim

### Find the % variance explained by each eigenvector 
diag(CVbody_pca_trim$Eval)/sum(CVbody_pca_trim$Eval)*100

CVbody_mat_trim<-matrix(as.numeric(as.matrix(do.call(rbind, lapply(morpho_trim[-c(grep("forficatus",names(morpho_trim)),grep("savana",names(morpho_trim)))],function(x) x[,5:7])))),nrow=nrow(do.call(rbind,morpho_trim[-c(grep("forficatus",names(morpho_trim)),grep("savana",names(morpho_trim)))]))) 

CVbody_scores_trim<-phytools:::scores(CVbody_pca_trim,newdata=CVbody_mat_trim)
rownames(CVbody_scores_trim)<-rownames(do.call(rbind, lapply(morpho_trim[-c(grep("forficatus",names(morpho_trim)),grep("savana",names(morpho_trim)))],function(x) x[,5:7])))
CVbody_scores_trim_df<-data.frame(CVbody_scores_trim)
CVbody_scores_trim_df$otu<-factor(sapply(strsplit(rownames(CVbody_scores_trim_df),"[.]"),function(x) x[1]))
CVbody_scores_split_trim<-split(CVbody_scores_trim_df, CVbody_scores_trim_df$otu)
CVbody_scores_trim_cv<-t(sapply(CVbody_scores_split_trim,function(x) c(mean(na.omit(x[,1])),mean(na.omit(x[,2])))))

### Build Ancestral State Reconstruction objects for phylogeny in figure
is.binary(phy)
phy<-multi2di(phy) #resolve polytomies
is.binary(phy)

### Migration strategy
### Load the dataset to save as a factor (continuous var as vector below)
Tyrannus.data<-read.csv('./Output Files/Tyrannus morphology + PCA avg_Adults.csv', row.names = 1)

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#Tyrannus.data<-read.csv('./Output Files/Tyrannus morphology + PCA avg_Females_Adults.csv', row.names = 1)
#OR
#Tyrannus.data<-read.csv('./Output Files/Tyrannus morphology + PCA avg_Males_Adults.csv', row.names = 1)

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

species<-rownames(Tyrannus.data)
Tyrannus.data<-cbind(species,Tyrannus.data)

### Check that tree$tip.label is the same as otu_avg
phy$tip.label[!phy$tip.label %in% Tyrannus.data$species] #should return: charactor(0)

### Use diversitree, the mk-n model to reconstruct marginal ancestral states
require(diversitree)

Migration<-factor(Tyrannus.data$Strategy)
names(Migration)<-rownames(Tyrannus.data)

#Change character to numeric, because this is the format diversitree accepts
char1<-as.numeric(Migration)
names(char1)<-names(Migration)

is.ultrametric(phy)
phy<-force.ultrametric(phy,method="nnls")

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

WC<-Tyrannus.data$WC.Average
names(WC)<- rownames(Tyrannus.data)
View(WC)

obj<-contMap(phy, WC)
n<-length(obj$cols) #get the length of the color ramp
obj$cols[1:n]<-grey(0:(n-1)/(n-1)) #change the color ramp to a grey scale

### Build a figure to show PPCA results next to the phylogeny ###
## Set up color labels ##
coldf<-data.frame(row.names=sort(phy$tip.label),col=rep(NA,Ntip(phy)),pch=rep(NA,Ntip(phy)))

species<-unique(sapply(strsplit(rownames(coldf),"_"),function(x) paste(x[c(1,2)],collapse="_")))
ssp_count<-table(sapply(strsplit(rownames(coldf),"_"),function(x) paste(x[c(1,2)],collapse="_")))

sp_col1<-rainbow(n=(length(species)+2))[1:length(species)]
sp_col2<-c(brewer.pal(12,"Set3"),"#D2691E")
sp_col2[2]<-"#663399"

col_vec<-vector()
pch_vec<-vector()
for(i in 1:length(ssp_count)){
  col_vec<-c(col_vec,rep(sp_col2[i],ssp_count[i]))
  pch_vec <-c(pch_vec,c(21:25,3:4)[1:ssp_count[i]])
}
coldf$col<-col_vec
coldf$pch<-pch_vec

mig<-read.csv("Tyrannus_subspecies_MigrationStrategies.csv", row.names = 1)
coldf2<-merge(coldf, mig, by=0)
rownames(coldf)<-coldf2$Row.names
coldf$Row.names<-NULL

coldf$Strategy<-factor(coldf2$Strategy,levels=c("sedentary","partial","migratory")) #Convert to a factor

### Use outline to indicate migratory strategy, no outline = non migratory, dashed= partial, complete = migrant ###
coldf$outlinelwd<-c(0.5,1.5,2.5)[as.numeric(coldf$Strategy)]
coldf$outlinecolor<-c("gray70","gray40","black")[as.numeric(coldf$Strategy)]

### Fix pch 3 and 4 so they plot correctly on phylogeny ###
coldf$outlinecolor[coldf$pch %in% c(3,4)]<-coldf$col[coldf$pch %in% c(3,4)]
coldf$col[coldf$pch %in% c(3,4)]<-NA

png(file="Tyrannus_phylogenetic_PCA_Adults.png",width=7.5,height=8,units="in",res=500)

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#png(file="Tyrannus_phylogenetic_PCA_Females_Adults.png",width=7.5,height=8,units="in",res=500)
#OR
#png(file="Tyrannus_phylogenetic_PCA_Males_Adults.png",width=7.5,height=8,units="in",res=500)

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

layout(matrix(c(1,1,1,5,2,3,4,5),nrow=4), widths=c(2.75,1,1,1,3.75))
layout.show(n=5)

par(mar=c(2,0.5,0,14))

plot(obj, legend=FALSE, ftype="off", mar=c(5.1,0.4,0.4,13), fsize=c(0.7,0.9))
add.color.bar(2.5,obj$cols,title="Wing Chord Length",
              lims=obj$lims,digits=0,prompt=FALSE,x=0.15,
              y=2.5-0.08*(Ntip(obj$tree)-1),lwd=4,fsize=0.8,subtitle="")

co<-c("black","gray40","gray70") #Vector of colours for migration strategy
nodelabels(pie=st, piecol=co, cex=1.5, col="gray32") #plot the probabilities of each state at each node

par(xpd=NA)

format_tl<-phy$tip.label
format_tl<-gsub("_"," ",format_tl)
format_tl<-gsub("Tyrannus","T.",format_tl)

text(x=rep(5.2,Ntip(phy)),y=1:Ntip(phy),label=format_tl,adj=c(0,0.5),font=3)

points(x=rep(4.95,Ntip(phy)),y=1:Ntip(phy),cex=2,pch=coldf[phy$tip.label,]$pch,bg=coldf[phy$tip.label,]$col, col=coldf[phy$tip.label,]$outlinecolor,lwd=coldf[phy$tip.label,]$outlinelwd,lty=3)

### Time axis ###
axisPhylo(1)
text(x=4.7,y=-2.5,label="mya",cex=1,pch=22)

### Add legend ###
table(coldf$Strategy)#Counts for each migratory strategy, used in legend
points(x=rep(0.25,3),y=c(2,3,4),pch=21,bg="gray85",lwd=c(0.5,1.5,2.5),col=c("gray70","gray40","black"),cex=2)
text(x=rep(0.5,3),y=c(2,3,4),label=c("sedentary (15)","partial migrant (5)","migratory (8)"),adj=c(0,0.5))

### Add figure label ###
text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.05,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.05,label=LETTERS[1],font=2,cex=1.5)

### PCA Panel plots ###
par(mar=c(2.75,2.75,1.5,1.5))

plot(bill_scores_df[,1],bill_scores_df[,2],bg=gsub("NA60",NA,paste0(coldf[bill_scores_df$otu,]$col,"60")),col= gsub("NA60","#663399",paste0(coldf[bill_scores_df$otu,]$col,"60")),pch= coldf[bill_scores_df$otu,]$pch,axes=F,xlab="",ylab="")

points(bill_scores_avg[,1], bill_scores_avg[,2],pch=coldf[rownames(bill_scores_avg),]$pch,bg= coldf[rownames(bill_scores_avg),]$col,col=coldf[rownames(bill_scores_avg),]$outlinecolor,cex=2,lwd=coldf$outlinelwd)

box()
title(main="Bill PPCA")

axis(1,mgp=c(0,0,0),tck=-0.025,cex.axis=0.45)
axis(2,mgp=c(0,0.2,0),tck=-0.025,cex.axis=0.45)

mtext(at=-2, text=paste0("short, narrow, shallow"), side=1, line=0.5, cex=0.45)
mtext(at=11, text=paste0("long, wide, deep"), side=1, line=0.5, cex=0.45)
mtext(at=1.75, text=paste0("long, narrow, shallow"), side=2, line=0.7, cex=0.45)
mtext(at=-2.75, text=paste0("short, wide, deep"), side=2, line=0.7, cex=0.45)

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines of code instead of the above lines, then proceed.
#mtext(at=-4, text=paste0("short, narrow, shallow"), side=1, line=0.5, cex=0.45)
#mtext(at=11, text=paste0("long, wide, deep"), side=1, line=0.5, cex=0.45)
#mtext(at=1.75, text=paste0("long, narrow, shallow"), side=2, line=0.7, cex=0.45)
#mtext(at=-2, text=paste0("short, wide, deep"), side=2, line=0.7, cex=0.45)
#OR
#mtext(at=-10, text=paste0("short, narrow, shallow"), side=1, line=0.5, cex=0.45)
#mtext(at=2, text=paste0("long, wide, deep"), side=1, line=0.5, cex=0.45)
#mtext(at=1.75, text=paste0("long, narrow, shallow"), side=2, line=0.7, cex=0.45)
#mtext(at=-2.75, text=paste0("short, wide, deep"), side=2, line=0.7, cex=0.45)

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

mtext(text=paste0("PC1 (",round(100*(diag(bill_pca$Eval)[1]/sum(diag(bill_pca$Eval))),2),"%)"),side=1,line=1.3,cex=0.55)
mtext(text=paste0("PC2 (",round(100*(diag(bill_pca$Eval)[2]/sum(diag(bill_pca$Eval))),2),"%)"),side=2,line=1.3,cex=0.55)

text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.075,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.1,label=LETTERS[2],font=2,cex=1.5)

### Feather PPCA All data ###
plot(body_scores_df[,1],body_scores_df[,2],bg=gsub("NA60",NA,paste0(coldf[body_scores_df$otu,]$col,"60")),col= gsub("NA60","#663399",paste0(coldf[body_scores_df$otu,]$col,"60")),pch= coldf[body_scores_df$otu,]$pch,axes=F,xlab="",ylab="")

### add rectangle to indicate zoom for panel D ###
body_scores_df2<-subset(body_scores_df, !(otu=="Tyrannus_forficatus"|otu=="Tyrannus_savana_circumdatus"|otu=="Tyrannus_savana_savana"|otu=="Tyrannus_savana_monachus_CA"|otu=="Tyrannus_savana_monachus_SA"|otu=="Tyrannus_savana_sanctaemartae"))

range(na.omit(body_scores_df2[,1])) #Use these values to draw box in Feather PPCA plot 2 for zooming in
range(na.omit(body_scores_df2[,2])) #Use these values to draw box in Feather PPCA plot 2 for zooming in
rect(xleft= 22.49359,ybottom=-29.45303,xright= 58.11388,ytop= 18.04108,lty=3) ### see range body scores avg 2 below to set these properly

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#rect(xleft= 12.29256,ybottom=-30.44410,xright= 45.24831,ytop= 12.25531,lty=3) ### see range body scores avg 2 below to set these properly
#OR
#rect(xleft= 31.12105,ybottom=-27.00028,xright= 66.87721,ytop= 18.12568,lty=3) ### see range body scores avg 2 below to set these properly

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

box()
title(main="Feather PPCA")

points(body_scores_avg[,1], body_scores_avg[,2],pch=coldf[rownames(body_scores_avg),]$pch,bg= coldf[rownames(body_scores_avg),]$col,col=coldf[rownames(body_scores_avg),]$outlinecolor,cex=2,lwd=coldf$outlinelwd)

axis(1,mgp=c(0,0.25,0),tck=-0.025,cex.axis=0.45)
axis(2,mgp=c(0,0.25,0),tck=-0.025,cex.axis=0.45)

mtext(at=-160, text=paste0("long feathers"), side=1, line=0.7, cex=0.45)
mtext(at=50, text=paste0("short feathers"), side=1, line=0.7, cex=0.45)
mtext(at=-25, text=paste0("long feathers"), side=2, line=0.7, cex=0.45)
mtext(at=10, text=paste0("short feathers"), side=2, line=0.7, cex=0.45)

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines of code instead of the above lines, then proceed.
#mtext(at=-125, text=paste0("long feathers"), side=1, line=0.7, cex=0.45)
#mtext(at=25, text=paste0("short feathers"), side=1, line=0.7, cex=0.45)
#mtext(at=-25, text=paste0("long feathers"), side=2, line=0.7, cex=0.45)
#mtext(at=10, text=paste0("short feathers"), side=2, line=0.7, cex=0.45)#OR
#OR
#mtext(at=-160, text=paste0("long feathers"), side=1, line=0.7, cex=0.45)
#mtext(at=50, text=paste0("short feathers"), side=1, line=0.7, cex=0.45)
#mtext(at=-25, text=paste0("long feathers"), side=2, line=0.7, cex=0.45)
#mtext(at=10, text=paste0("short feathers"), side=2, line=0.7, cex=0.45)

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

mtext(text=paste0("PC1 (",round(100*(diag(body_pca$Eval)[1]/sum(diag(body_pca$Eval))),2),"%)"),side=1,line=1.3,cex=0.55)
mtext(text=paste0("PC2 (",round(100*(diag(body_pca$Eval)[2]/sum(diag(body_pca$Eval))),2),"%)"),side=2,line=1.3,cex=0.55)

## Figure panel label ###
text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.075,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.1,label=LETTERS[3],font=2,cex=1.5)

### Feather PPCA Zoomed in ###
plot(body_scores_df2[,1],body_scores_df2[,2],bg=gsub("NA60",NA,paste0(coldf[body_scores_df2$otu,]$col,"60")),col= gsub("NA60","#663399",paste0(coldf[body_scores_df2$otu,]$col,"60")),pch= coldf[body_scores_df2$otu,]$pch,axes=F,xlab="",ylab="")

box()
title(main="Feather PPCA zoomed in")

body_scores_avg2<-body_scores_avg[-c(15,20,21,22,23,24),] 

points(body_scores_avg2[,1], body_scores_avg2[,2],pch=coldf[rownames(body_scores_avg2),]$pch,bg= coldf[rownames(body_scores_avg2),]$col,col=coldf[rownames(body_scores_avg2),]$outlinecolor,cex=2,lwd=coldf[rownames(body_scores_avg2),]$outlinelwd)

axis(1,mgp=c(0,0,0),tck=-0.025,cex.axis=0.45)
axis(2,mgp=c(0,0.25,0),tck=-0.025,cex.axis=0.45)

mtext(at=27.5, text=paste0("long feathers"), side=1, line=0.5, cex=0.45)
mtext(at=57.5, text=paste0("short feathers"), side=1, line=0.5, cex=0.45)
mtext(at=-25, text=paste0("long feathers"), side=2, line=0.7, cex=0.45)
mtext(at=7.75, text=paste0("short feathers"), side=2, line=0.7, cex=0.45)

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines of code instead of the above lines, then proceed.
#mtext(at=14, text=paste0("long feathers"), side=1, line=0.5, cex=0.45)
#mtext(at=40, text=paste0("short feathers"), side=1, line=0.5, cex=0.45)
#mtext(at=-25, text=paste0("long feathers"), side=2, line=0.7, cex=0.45)
#mtext(at=5, text=paste0("short feathers"), side=2, line=0.7, cex=0.45)
#OR
#mtext(at=35, text=paste0("long feathers"), side=1, line=0.5, cex=0.45)
#mtext(at=60, text=paste0("short feathers"), side=1, line=0.5, cex=0.45)
#mtext(at=-25, text=paste0("long wings, short tails"), side=2, line=0.7, cex=0.45)
#mtext(at=7.75, text=paste0("short wings, long tails"), side=2, line=0.7, cex=0.45)

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

mtext(text=paste0("PC1 (",round(100*(diag(body_pca$Eval)[1]/sum(diag(body_pca$Eval))),2),"%)"),side=1,line=1.3,cex=0.55)
mtext(text=paste0("PC2 (",round(100*(diag(body_pca$Eval)[2]/sum(diag(body_pca$Eval))),2),"%)"),side=2,line=1.3,cex=0.55)

## Figure panel label ###
text(x=par("usr")[1]+diff(c(par("usr")[1],par("usr")[2]))*0.075,y=par("usr")[4]-diff(c(par("usr")[3],par("usr")[4]))*0.1,label=LETTERS[4],font=2,cex=1.5)

### Bird Pictures Panel ###
#par(mar=c(2,0.5,0,14))

require(magick)

# Crop images
#img<-image_crop(example, "100x150+50")#crop out width:100px and height:150px starting +50px from the left

T_savana<-magick::image_read("50763510127_317508235f_o.jpg")
T_savana<-image_scale(T_savana, "200") # width: 200px
T_savana<-image_crop(T_savana, "200x300+0+20")
T_savana

T_caudifasciatus<-magick::image_read("T_caudifasciatus.png")
T_caudifasciatus<-image_scale(T_caudifasciatus,"600")
T_caudifasciatus<-image_crop(T_caudifasciatus, "300x280+120+10")
T_caudifasciatus

T_crassirostris<-magick::image_read("50250535916_399b2a0dea_o.jpg")
T_crassirostris<-image_scale(T_crassirostris,"275")
T_crassirostris<-image_crop(T_crassirostris,"250x280+40+0")
T_crassirostris

T_cubensis<-magick::image_read("Giant Kingbird 2 HR.jpg")
T_cubensis<-image_scale(T_cubensis,"450")
T_cubensis<-image_crop(T_cubensis,"210x280+190")
T_cubensis

photos<-image_append(c(T_savana, T_caudifasciatus, T_crassirostris, T_cubensis))

par(mar=c(0,0,0,0))

plot(photos)

points(x=180, y=260, cex=2, pch=coldf[phy$tip.label,]$pch[26],bg=coldf[phy$tip.label,]$col[26], col=coldf[phy$tip.label,]$outlinecolor[26],lwd=coldf[phy$tip.label,]$outlinelwd[26],lty=3)
points(x=480, y=260, cex=2, pch=coldf[phy$tip.label,]$pch[21],bg=coldf[phy$tip.label,]$col[21], col=coldf[phy$tip.label,]$outlinecolor[21],lwd=coldf[phy$tip.label,]$outlinelwd[21],lty=3)
points(x=720, y=260, cex=2, pch=coldf[phy$tip.label,]$pch[14],bg=coldf[phy$tip.label,]$col[14], col=coldf[phy$tip.label,]$outlinecolor[14],lwd=coldf[phy$tip.label,]$outlinelwd[14],lty=3)
points(x=930, y=260, cex=2, pch=coldf[phy$tip.label,]$pch[4],bg=coldf[phy$tip.label,]$col[4], col=coldf[phy$tip.label,]$outlinecolor[4],lwd=coldf[phy$tip.label,]$outlinelwd[4],lty=3)

text(x=20, y=260, paste0("E"), font=2, cex=1.5)
text(x=220, y=260, paste0("F"),font=2, cex=1.5)
text(x=520, y=260, paste0("G"),font=2, cex=1.5)
text(x=750, y=260, paste0("H"),font=2, cex=1.5)

dev.off()

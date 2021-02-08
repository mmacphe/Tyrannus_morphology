rm(list=ls())

require(phytools)
require(RColorBrewer)

### Read in average morphology data from .csv ###
morpho<-read.csv("Tyrannus_data.csv", row.names=2)

### Read in phylogeny 
phy<-read.tree("Tyrannus_phylogeny.tre")

### Check that tree$tip.label is the same as morpho
phy$tip.label[!phy$tip.label %in% rownames(morpho)]

### Bill Phylogenetic PCA ###
bill_pca<-phyl.pca(phy,morpho[,2:4],method="lambda")
bill_pca

### Bring in morphology data set for all individuals to add points onto figure
morpho_whole<-read.csv("Tyrannus morphology data.csv", row.names = 1)

### Match tip.labels with phylogeny tip.labels ###
morpho_whole$tip.label
names(table(morpho_whole$tip.label))

morpho_whole$tip.label<-gsub(" ","_",morpho_whole$tip.label)
names(table(morpho_whole$tip.label))

phy$tip.label

phy$tip.label[!phy$tip.label %in% morpho_whole$tip.label]
morpho_whole$tip.label[! morpho_whole$tip.label %in% phy$tip.label]

### Split the morphology data sets up into categories by tip.label
morpho_split<-split(morpho_whole,morpho_whole$tip.label)
names(morpho_split)
phy$tip.label

phy$tip.label[!phy$tip.label %in% names(morpho_split)]
names(morpho_split)[! names(morpho_split) %in% phy$tip.label]

morpho_trim<-lapply(morpho_split,function(x) x[colnames(x) %in% c("BL.Average","BW.Average","BD.Average", "Kipp.s.Average", "WC.Average", "Tail", "Tarsus.Average", "name_num")])

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
body_pca<-phyl.pca(phy,morpho[,6:8],method="lambda")
body_pca

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
body_pca_trim<-phyl.pca(phy_trim,otu_avg_trim[,6:8],method="lambda")
body_pca_trim

body_mat_trim<-matrix(as.numeric(as.matrix(do.call(rbind, lapply(morpho_trim[-c(grep("forficatus",names(morpho_trim)),grep("savana",names(morpho_trim)))],function(x) x[,5:7])))),nrow=nrow(do.call(rbind,morpho_trim[-c(grep("forficatus",names(morpho_trim)),grep("savana",names(morpho_trim)))]))) 

body_scores_trim<-phytools:::scores(body_pca_trim,newdata=body_mat_trim)
rownames(body_scores_trim)<-rownames(do.call(rbind, lapply(morpho_trim[-c(grep("forficatus",names(morpho_trim)),grep("savana",names(morpho_trim)))],function(x) x[,5:7])))
body_scores_trim_df<-data.frame(body_scores_trim)
body_scores_trim_df$otu<-factor(sapply(strsplit(rownames(body_scores_trim_df),"[.]"),function(x) x[1]))
body_scores_split_trim<-split(body_scores_trim_df, body_scores_trim_df$otu)
body_scores_trim_avg<-t(sapply(body_scores_split_trim,function(x) c(mean(na.omit(x[,1])),mean(na.omit(x[,2])))))

#I don't think we do the phylANOVA on the trimmed dataset, and some extra code will be needed to match them to the right OTU bc there are 6 fewer OTUs in body_scores_trim_avg
# View(bill_scores_avg)
# body_scores_trim_avg<-as.data.frame(body_scores_trim_avg)
# morpho$Body_trim_PC1<-body_scores_trim_avg$V1
# morpho$BillPC2<-bill_scores_avg$V2
# View(morpho)

write.csv(morpho, file="Tyrannus morphology + PCA avg.csv")

### Bring in cv summary dataset 
cv_summary<-read.csv("cv_summary_table.csv", row.names = 2)

### Coefficient of Variation (CV) Bill Phylogenetic PCA ###
CVbill_pca<-phyl.pca(phy,cv_summary[,2:4],method="lambda")
CVbill_pca

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
CVbody_pca<-phyl.pca(phy,cv_summary[,6:8],method="lambda")
CVbody_pca

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

write.csv(cv_summary, file="cv_summary.csv")

### CV Body PCA without T. forficatus or T. savana bc their tails are so long ###
phy_trim<-drop.tip(phy,c(grep("forficatus",phy$tip.label),grep("savana",phy$tip.label)))
otu_cv_trim<-cv_summary[rownames(cv_summary) %in% phy_trim$tip.label,]
CVbody_pca_trim<-phyl.pca(phy_trim,otu_cv_trim[,6:8],method="lambda")
CVbody_pca_trim

CVbody_mat_trim<-matrix(as.numeric(as.matrix(do.call(rbind, lapply(morpho_trim[-c(grep("forficatus",names(morpho_trim)),grep("savana",names(morpho_trim)))],function(x) x[,5:7])))),nrow=nrow(do.call(rbind,morpho_trim[-c(grep("forficatus",names(morpho_trim)),grep("savana",names(morpho_trim)))]))) 

CVbody_scores_trim<-phytools:::scores(CVbody_pca_trim,newdata=CVbody_mat_trim)
rownames(CVbody_scores_trim)<-rownames(do.call(rbind, lapply(morpho_trim[-c(grep("forficatus",names(morpho_trim)),grep("savana",names(morpho_trim)))],function(x) x[,5:7])))
CVbody_scores_trim_df<-data.frame(CVbody_scores_trim)
CVbody_scores_trim_df$otu<-factor(sapply(strsplit(rownames(CVbody_scores_trim_df),"[.]"),function(x) x[1]))
CVbody_scores_split_trim<-split(CVbody_scores_trim_df, CVbody_scores_trim_df$otu)
CVbody_scores_trim_cv<-t(sapply(CVbody_scores_split_trim,function(x) c(mean(na.omit(x[,1])),mean(na.omit(x[,2])))))

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
  pch_vec <-c(pch_vec,c(21:25,21:22)[1:ssp_count[i]])
}
coldf$col<-col_vec
coldf$pch<-pch_vec

mig<-read.csv("Tyrannus_subspecies_MigrationStrategies.csv", row.names = 1)
coldf2<-merge(coldf, mig, by=0)
pch2<-c(21,22,23)[as.numeric(coldf2$Strategy)]

coldf2$pch2<-c(21,22,23)[as.numeric(coldf2$Strategy)]

png(file="Tyrannus_phylogenetic_PCA.png",width=6.5,height=5.5,units="in",res=500)

layout(matrix(c(1,1,1,2,3,4),nrow=3),widths=c(2.25,1))
layout.show(n=4)

par(mar=c(2,0,0,14))

plot(phy,cex=0.8,show.tip.label=F)
phy$edge[,2][phy$edge[,2]<=Ntip(phy)]

par(xpd=NA)

format_tl<-phy$tip.label
format_tl<-gsub("_"," ",format_tl)
format_tl<-gsub("Tyrannus","T.",format_tl)

text(x=rep(5.2,Ntip(phy)),y=1:Ntip(phy),label=format_tl,adj=c(0,0.5),font=3)

#text(x=rep(4.8,Ntip(phy)),y=1:Ntip(phy),label=morpho[phy$tip.label,]$Strategy,cex=1,pch=22)

#points(x=rep(4.8,Ntip(phy)),y=1:Ntip(phy),bg=coldf2[phy$tip.label,]$col,levels=c("sedentary","partial","migratory"),cex=1.75,pch=c(21,22,23)[as.numeric(coldf2[phy$tip.label,]$Strategy)])

points(x=rep(4.8,Ntip(phy)),y=1:Ntip(phy),bg=c("white","gray","black")[factor(morpho[phy$tip.label,]$Strategy,levels=c("sedentary","partial","migratory"))],cex=1.75,pch=22)
points(x=rep(5.05,Ntip(phy)),y=1:Ntip(phy), pch=c(21,22,23)[as.numeric(coldf2[phy$tip.label,]$Strategy)],bg=coldf2[phy$tip.label,]$col,cex=1.75)

### Time axis ###
axisPhylo(1)
text(x=5.1,y=-0.95,label="mya",cex=1,pch=22)

par(mar=c(2.75,2.75,1.5,1.5))

plot(bill_scores_df[,1],bill_scores_df[,2],bg=paste0(coldf[bill_scores_df$otu,]$col,"90"),col= paste0(coldf[bill_scores_df$otu,]$col,90),pch= coldf[bill_scores_df$otu,]$pch,axes=F,xlab="",ylab="")
points(bill_scores_avg[,1], bill_scores_avg[,2],pch=coldf$pch,col="black",bg=paste0(coldf$col,95),cex=2)
box()
title(main="Bill PPCA")

axis(1,mgp=c(0,0.25,0),tck=-0.025,cex.axis=0.75)
axis(2,mgp=c(0,0.25,0),tck=-0.025,cex.axis=0.75)
mtext(text=paste0("PC1 (",round(100*(diag(bill_pca$Eval)[1]/sum(diag(bill_pca$Eval))),2),"%)"),side=1,line=1.3,cex=0.55)
mtext(text=paste0("PC2 (",round(100*(diag(bill_pca$Eval)[2]/sum(diag(bill_pca$Eval))),2),"%)"),side=2,line=1.3,cex=0.55)

plot(body_scores_df[,1],body_scores_df[,2],bg=paste0(coldf[body_scores_df$otu,]$col,"90"),col= paste0(coldf[body_scores_df$otu,]$col,90),pch= coldf[body_scores_df$otu,]$pch,axes=F,xlab="",ylab="")
title(main="Feather PPCA")

points(body_scores_avg[,1], body_scores_avg[,2],pch=coldf$pch,col="black",bg=paste0(coldf$col,95),cex=2)
box()
axis(1,mgp=c(0,0.25,0),tck=-0.025,cex.axis=0.75)
axis(2,mgp=c(0,0.25,0),tck=-0.025,cex.axis=0.75)
mtext(text=paste0("PC1 (",round(100*(diag(body_pca$Eval)[1]/sum(diag(body_pca$Eval))),2),"%)"),side=1,line=1.3,cex=0.55)
mtext(text=paste0("PC2 (",round(100*(diag(body_pca$Eval)[2]/sum(diag(body_pca$Eval))),2),"%)"),side=2,line=1.3,cex=0.55)

body_scores_df2<-subset(body_scores_df, !(otu=="Tyrannus_forficatus"|otu=="Tyrannus_savana_circumdatus"|otu=="Tyrannus_savana_savana"|otu=="Tyrannus_savana_monachus_CA"|otu=="Tyrannus_savana_monachus_SA"|otu=="Tyrannus_savana_sanctaemartae"))
                                          
plot(body_scores_df2[,1],body_scores_df2[,2],bg=paste0(coldf[body_scores_df2$otu,]$col,"90"),col= paste0(coldf[body_scores_df2$otu,]$col,90),pch= coldf[body_scores_df2$otu,]$pch,axes=F,xlab="",ylab="", xlim=c(20,65), ylim=c(-30,20))
title(main="Feather PPCA zoomed")

body_scores_avg2<-body_scores_avg[-c(15,20,21,22,23,24),] 

points(body_scores_avg2[,1], body_scores_avg2[,2],pch=coldf$pch,col="black",bg=paste0(coldf$col,95),cex=2)
box()
axis(1,mgp=c(0,0.25,0),tck=-0.025,cex.axis=0.75)
axis(2,mgp=c(0,0.25,0),tck=-0.025,cex.axis=0.75)
mtext(text=paste0("PC1 (",round(100*(diag(body_pca$Eval)[1]/sum(diag(body_pca$Eval))),2),"%)"),side=1,line=1.3,cex=0.55)
mtext(text=paste0("PC2 (",round(100*(diag(body_pca$Eval)[2]/sum(diag(body_pca$Eval))),2),"%)"),side=2,line=1.3,cex=0.55)

dev.off()

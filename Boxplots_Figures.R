rm(list=ls())

require(phytools)

setwd("~/Downloads/Tyrannus/4_phylANOVA")

phy<-read.tree("Tyrannus_phylogeny.tre")

### Bring in morphology dataset with average morphology values for each OTU and PC scores
otu_avg<-read.csv("Tyrannus morphology + PCA avg.csv", row.names = 1)
species<-rownames(otu_avg)
otu_avg<-cbind(species,otu_avg)
View(otu_avg)

### Check that tree$tip.label is the same as otu_avg
phy$tip.label[!phy$tip.label %in% otu_avg$species]

### phylANOVA on OTU morphometric means while controlling for body size ###
response_var<-colnames(otu_avg)[c(3,4,5,7,8,9)]

### Set up predictor variable for phylANOVA ###
migration<-otu_avg$Strategy
names(migration)<-rownames(otu_avg)

### Set up tarsus for size correction #
otu_avg$Tarsus.Average <-as.numeric(otu_avg$Tarsus.Average)
tarsus.length<-otu_avg$Tarsus.Average
names(tarsus.length)<-rownames(otu_avg)

#Create empty list to store output of each for loop iteration
phyl.resid_output<-list()
phylANOVA_output<-list() 

#Conduct the phylogenetic ANOVA he following code works without assuming y is the same order as tree$tip.label
for(i in 1:length(response_var)){
	print(paste0("Character ",i," -- ",response_var[i]))#Report character
	character_of_interest<-otu_avg[,response_var[i]] #Extract character of interest (coi)
	names(character_of_interest)<-rownames(otu_avg) #Add names to vector of coi
	phyl.resid_output[[i]]<-phyl.resid(phy,character_of_interest,tarsus.length,method="lambda") #Perform phylogenetic residuals
	phylANOVA_output[[i]]<-phylANOVA(phy, migration, phyl.resid_output[[i]]$resid[,1])
}

names(phyl.resid_output)<-response_var
names(phylANOVA_output)<-response_var 

### Write out output ###
#sink("phylANOVA_posthoc_output.txt")
sink("phylANOVA_output.txt")
for(i in 1:length(phylANOVA_output)){
  cat(paste0("Character ",i," -- ",response_var[i]))
  cat("\n\n")
  print(phylANOVA_output[[i]])
  cat("#############")
  cat("\n\n")
}
sink()

### Conduct the phylogenetic ANOVA for PC scores that are not size-corrected ###
phy$tip.label[!phy$tip.label %in% rownames(otu_avg)]
response_var2<-colnames(otu_avg)[c(11,12,13,14)]
#migration<-otu_avg$Strategy
#names(migration)<-rownames(otu_avg)

phylANOVA_PCoutput<-list() 
#Question: How to label y? (the tip.labels so that it doesn't assume same order as tree$tip.label)
for(i in 1:length(response_var2)){
  print(paste0("Character ",i," -- ",response_var2[i]))#Report character
  character_of_interest<-otu_avg[,response_var2[i]] #Extract character of interest (coi)
  names(character_of_interest)<-rownames(otu_avg) #Add names to vector of coi
  phylANOVA_PCoutput[[i]]<-phylANOVA(tree=phy, x=migration, y=character_of_interest) 
}

names(phylANOVA_PCoutput)<-response_var2 

### Write out output ###
sink("phylANOVA_output_PCscores.txt")

for(i in 1:length(phylANOVA_PCoutput)){
	cat(paste0("Character ",i," -- ",response_var2[i]))
	cat("\n\n")
	print(phylANOVA_PCoutput[[i]])
	cat("#############")
	cat("\n\n")
}
sink()

### Repeat for CV data: Bring in cv_summary_table_v4
otu_cv_avg<-read.csv("cv_summary_v5.csv", row.names = 1) 
species<-rownames(otu_cv_avg)
otu_cv_avg<-cbind(species,otu_cv_avg)
View(otu_cv_avg)

### Check that tree$tip.label is the same as otu_avg
phy$tip.label[!phy$tip.label %in% otu_cv_avg$species]

response_var3<-colnames(otu_cv_avg)[c(3,4,5,6,7,8,9,11,12,13,14)]

### Set up predictor variable for phylANOVA ###
migration<-otu_cv_avg$Strategy
names(migration)<-rownames(otu_cv_avg)

#Create empty list to store output of each for loop iteration
phylANOVA_CVoutput<-list() 

### Conduct the phylogenetic ANOVA on coefficients of variation ###
for(i in 1:length(response_var3)){
  print(paste0("Character ",i," -- ",response_var3[i]))#Report character
  character_of_interest<-otu_cv_avg[,response_var3[i]] #Extract character of interest (coi)
  names(character_of_interest)<-rownames(otu_cv_avg) #Add names to vector of coi
  phylANOVA_CVoutput[[i]]<-phylANOVA(phy, migration, y=character_of_interest)
}

names(phylANOVA_CVoutput)<-response_var3

### Write out output ###
sink("phylANOVA_cv_output.txt")

for(i in 1:length(phylANOVA_CVoutput)){
  cat(paste0("Character ",i," -- ",response_var3[i]))
  cat("\n\n")
  print(phylANOVA_CVoutput[[i]])
  cat("#############")
  cat("\n\n")
}
sink()

### Bar plots of phenotypic characters (not size-corrected i.e., outputs of phylresid) ###
add_label_legend <- function(pos = "topleft", label, ...) {
    legend(pos, label, bty = "n", ...)
}

png(file="Phenotypes and Variation1.png",width=7,height=5.5,units="in",res=500)

par(mfrow=c(4,3),mar = c(3, 5, 0.5, 0.1))
for(i in 1:12) {
  boxplot(otu_avg$WC.Average~migration, data=otu_avg, axes=FALSE, ylab="Wing Cord (mm)")
  axis(2)
  add_label_legend("topright", paste0("(a)"))
  
  boxplot(otu_avg$Kipp.s.Average~migration, data=otu_avg, axes=FALSE, ylab="Kipp's Index (mm)", xlab="")
  axis(2)
  add_label_legend("topright", paste0("(b)"))
  
  boxplot(otu_avg$Tail~migration, data=otu_avg, axes=FALSE, ylab="Tail Length (mm)")
  axis(2)
  add_label_legend("topright", paste0("(c)"))
  
  boxplot(otu_avg$BL.Average~migration, data=otu_avg, axes=FALSE, ylab="Bill Length (mm)")
  axis(2)
  add_label_legend("topright", paste0("(d)"))
  
  boxplot(otu_avg$BW.Average~migration, data=otu_avg, axes=FALSE, ylab="Bill Width (mm)")
  axis(2)
  add_label_legend("topright", paste0("(e)"))
  
  boxplot(otu_avg$BD.Average~migration, data=otu_avg, axes=FALSE, ylab="Bill Depth (mm)")
  axis(2)
  add_label_legend("topright", paste0("(f)"))
  
  boxplot(otu_cv_avg$WC.Average~migration, data=otu_cv_avg, axes=FALSE, ylab="Wing Cord Variation")
  axis(2)
  add_label_legend("topright", paste0("(g)"))
  
  boxplot(otu_cv_avg$Kipp.s.Average~migration, data=otu_cv_avg, axes=FALSE, ylab="Kipp's Index Variation")
  axis(2)
  add_label_legend("topright", paste0("(h)"))
  
  boxplot(otu_cv_avg$Tail~migration, data=otu_cv_avg, axes=FALSE, ylab="Tail Length Variation")
  axis(2)
  add_label_legend("topright", paste0("(i)"))
  
  boxplot(otu_cv_avg$BL.Average~migration, data=otu_cv_avg, axes=FALSE, ylab="Bill Length Variation")
  axis(2)
  axis(1, at=1:3, labels=c("migratory","partial","sedentary"), tick=FALSE)
  add_label_legend("topright", paste0("(j)"))
  
  boxplot(otu_cv_avg$BW.Average~migration, data=otu_cv_avg, axes=FALSE, ylab="Bill Width Variation")
  axis(2)
  axis(1,at=1:3, labels=c("migratory","partial","sedentary"), tick=FALSE)
  add_label_legend("topright", paste0("(k)"))
  
  boxplot(otu_cv_avg$BD.Average~migration, data=otu_cv_avg, axes=FALSE, ylab="Bill Depth Variation")
  axis(2)
  axis(1, at=1:3, labels=c("migratory","partial","sedentary"), tick=FALSE)
  add_label_legend("topright", paste0("(", letters[i], ")"))
}
dev.off()

### Build residuals dataframe for a figure (like above) that is size and phylogeny corrected
phyl_resid<-as.data.frame(phyl.resid_output$BL.Average,row.names = phy$tip.label)
phyl_resid$BL<-phyl_resid$resid
phyl_resid2<-as.data.frame(phyl.resid_output$BW.Average, row.names=phy$tip.label)
phyl_resid$BW<-phyl_resid2$resid
phyl_resid3<-as.data.frame(phyl.resid_output$BD.Average, row.names=phy$tip.label)
phyl_resid$BD<-phyl_resid3$resid
phyl_resid4<-as.data.frame(phyl.resid_output$Kipp.s.Average, row.names=phy$tip.label)
phyl_resid$KI<-phyl_resid4$resid
phyl_resid5<-as.data.frame(phyl.resid_output$WC.Average, row.names=phy$tip.label)
phyl_resid$WC<-phyl_resid5$resid
phyl_resid6<-as.data.frame(phyl.resid_output$Tail, row.names=phy$tip.label)
phyl_resid$TL<-phyl_resid6$resid

residuals<-merge(phyl_resid, migration, by=0, row.names=TRUE) #the taxonomic order can be checked visually by comparin residuals to phyl.resid_output$ coi
residuals$migration<-residuals$y

png(file="Phenotypes_Residuals_Variation.png", width=7, height=10, units="in", res=500)
par(mfrow=c(6,3),mar = c(0.4, 5, 0.7, 0.1))
for(i in 1:18) {
  boxplot(otu_avg$WC.Average~migration, data=otu_avg, axes=FALSE, ylab="Wing Cord")
  axis(2)
  mtext(at=c(2), text="Phenotype (mm)", side=3, line=-0.5, cex=0.65)
  add_label_legend("topleft", paste0("A"))
  
  boxplot(residuals$WC~migration, data=residuals, axes=FALSE, ylab="")
  axis(2)
  mtext(at=2, text="Residuals", side=3, line=-0.5, cex=0.65)
  add_label_legend("topleft", paste0("B"))
  
  boxplot(otu_cv_avg$WC.Average~migration, data=otu_cv_avg, axes=FALSE, ylab="")
  axis(2)
  mtext(at=2, text="Variation", side=3, line=-0.5, cex=0.65)
  add_label_legend("topleft", paste0("C"))
  
  
  boxplot(otu_avg$Kipp.s.Average~migration, data=otu_avg, axes=FALSE, ylab="Kipp's Index", xlab="")
  axis(2)
  add_label_legend("topleft", paste0("D"))
  
  boxplot(residuals$KI~migration, data=residuals, axes=FALSE, ylab="", xlab="")
  axis(2)
  add_label_legend("topleft", paste0("E"))
  
  boxplot(otu_cv_avg$Kipp.s.Average~migration, data=otu_cv_avg, axes=FALSE, ylab="")
  axis(2)
  add_label_legend("topleft", paste0("F"))
  
  
  boxplot(otu_avg$Tail~migration, data=otu_avg, axes=FALSE, ylab="Tail Length")
  axis(2)
  add_label_legend("topleft", paste0("G"))
  
  boxplot(residuals$TL~migration, data=residuals, axes=FALSE, ylab="")
  axis(2)
  add_label_legend("topleft", paste0("H"))
  
  boxplot(otu_cv_avg$Tail~migration, data=otu_cv_avg, axes=FALSE, ylab="")
  axis(2)
  add_label_legend("topleft", paste0("I"))
  
  boxplot(otu_avg$BL.Average~migration, data=otu_avg, axes=FALSE, ylab="Bill Length")
  axis(2)
  add_label_legend("topleft", paste0("J"))
  
  boxplot(residuals$BL~migration, data=residuals, axes=FALSE, ylab="")
  axis(2)
  add_label_legend("topleft", paste0("K"))
  
  boxplot(otu_cv_avg$BL.Average~migration, data=otu_cv_avg, axes=FALSE, ylab="")
  axis(2)
  add_label_legend("topleft", paste0("L"))
  
  boxplot(otu_avg$BW.Average~migration, data=otu_avg, axes=FALSE, ylab="Bill Width")
  axis(2)
  add_label_legend("topleft", paste0("M"))
  
  boxplot(residuals$BW~migration, data=residuals, axes=FALSE, ylab="")
  axis(2)
  add_label_legend("topleft", paste0("N"))
  
  boxplot(otu_cv_avg$BW.Average~migration, data=otu_cv_avg, axes=FALSE, ylab="")
  axis(2)
  add_label_legend("topleft", paste0("O"))
  
  boxplot(otu_avg$BD.Average~migration, data=otu_avg, axes=FALSE, ylab="Bill Depth")
  axis(2)
  axis(1, tick=FALSE)
  mtext(at=c(1:3), text=c("migratory", "partial", "sedentary"), side=1, line=-0.55, cex=0.65)
  add_label_legend("topleft", paste0("P"))
  
  boxplot(residuals$BD~migration, data=residuals, axes=FALSE, ylab="")
  axis(2)
  axis(1, tick=FALSE)
  mtext(at=c(1:3), text=c("migratory", "partial", "sedentary"), side=1, line=-0.55, cex=0.65)
  add_label_legend("topleft", paste0("Q"))
  
  boxplot(otu_cv_avg$BD.Average~migration, data=otu_cv_avg, axes=FALSE, ylab="")
  axis(2)
  axis(1, tick=FALSE)
  mtext(at=c(1:3), text=c("migratory", "partial", "sedentary"), side=1, line=-0.55, cex=0.65) 
  add_label_legend("topleft", paste0("", LETTERS[i], ""))
 }
dev.off()

png(file="PPC scores boxplots.png",width=7,height=5.5,units="in",res=500)
par(mfrow=c(2,2),mar = c(3, 5, 0.5, 0.1))
for(i in 1:4) {
  boxplot(otu_avg$BillPC1~migration, data=otu_avg, axes=FALSE, ylab="Bill PPC1")
  axis(2)
  add_label_legend("topright", paste0("A"))
  
  boxplot(otu_avg$BillPC2~migration, data=otu_avg, axes=FALSE, ylab="Bill PPC2")
  axis(2)
  text(x=2, y=2.71, "______________", cex=1.8)
  text(x=2, y=2.50, "p=0.035", cex=0.8)
  add_label_legend("topright", paste0("B"))
  
  boxplot(otu_avg$BodyPC1~migration, data=otu_avg, axes=FALSE, ylab="Body PPC1")
  axis(2)
  axis(1, at=1:3, labels=c("migratory","partial","sedentary"), tick=FALSE)
  add_label_legend("topright", paste0("C"))
  
  boxplot(otu_avg$BodyPC2~migration, data=otu_avg, axes=FALSE, ylab="Body PPC2")
  axis(2)
  axis(1, at=1:3, labels=c("migratory","partial","sedentary"), tick=FALSE)
  add_label_legend("topright", paste0("D"))
}
dev.off()

### Trying to plot residuals
par(mfrow=c(4,3))
par(mar=c(2,3,0.5,0.5))
for(i in 1:length(phyl.resid_output)){
 	coi_bp<-data.frame(migration,phyl.resid_output[[i]]$resid)
 	coi_bp$migration<-factor(coi_bp$migration,levels=c("S","P","M"))
 	colnames(coi_bp)[2]<-response_var[i]
 	boxplot(formula(paste0(response_var[i],"~migration")),data=coi_bp)
 }
# 

rm(list=ls())

require(phytools)

### Set source directory to the folder this file came from within RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Bring in morphology datasets with average morphology values for each OTU and PC scores
otu_avg<-read.csv('./Output Files/Tyrannus morphology + PCA avg.csv', row.names = 1)
otu_cv_avg<-read.csv('./Output Files/cv_summary.csv', row.names = 1) 
residuals<-read.csv('./Output Files/phylANOVA_tarsus-corrected_residuals.csv', row.names = 1)

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing sexes separately, use the following lines of code instead of the above 3 lines, then proceed.
#otu_avg<-read.csv('./Output Files/Tyrannus morphology + PCA avg_Females.csv', row.names = 1)
#otu_cv_avg<-read.csv('./Output Files/cv_summary_Females.csv', row.names = 1) 
#residuals<-read.csv('./Output Files/phylANOVA_tarsus-corrected_residuals_Females.csv', row.names = 1)
### OR
#otu_avg<-read.csv('./Output Files/Tyrannus morphology + PCA avg_Males.csv', row.names = 1)
#otu_cv_avg<-read.csv('./Output Files/cv_summary_Males.csv', row.names = 1) 
#residuals<-read.csv('./Output Files/phylANOVA_tarsus-corrected_residuals_Males.csv', row.names = 1)

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

### Box plots of phenotypic characters (not size-corrected i.e., outputs of phylresid) ###
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

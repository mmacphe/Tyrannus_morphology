rm(list=ls())

require(ggplot2)
require(tidyr)
require(cowplot)
require(gridExtra)
require(gridGraphics)

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

### Box plots of phenotypic characters ###
### Make the series of boxplots showing residuals of morphometrics
png(file="Phenotype Residuals.png",width=7,height=5.5,units="in",res=500)

responsevariable.labs<-c("Bill Length", "Bill Width", "Bill Depth", "Kipp's Index", "Wing Cord", "Tail Length")
names(responsevariable.labs)<- c("BL", "BW", "BD", "KI", "WC", "TL")

p<-residuals %>%
  pivot_longer(BL:TL, names_to = "responsevariable", values_to = "residuals") %>%
  ggplot(aes(y=residuals, x=migration)) +
  geom_boxplot() +
  facet_wrap(vars(responsevariable), ncol=3, labeller=labeller(responsevariable = responsevariable.labs)) +
  background_grid(major="none", minor="none") +
  panel_border() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y= "Residuals") 

draw_plot(plot, x = 0, y = 0, width = 1, height = 1)
ggdraw() +  
  draw_plot(p,0,0,1,1) +
  draw_plot_label(c("A", "B", "C", "D", "E", "F"), c(0.075,0.385,0.69,0.075,0.385,0.69), c(0.99,0.99,0.99,0.50,0.50,0.50), size=15)
dev.off()

### Make the box plots of variation in morphometrics
png(file="Phenotype Variation.png",width=7,height=5.5,units="in",res=500)
otu_cv_avg<-subset(otu_cv_avg, select=-c(Tarsus.Average)) #Remove tarsus from the df to be able to call BL:TL below

responsevariable.labs<-c("Bill Length", "Bill Width", "Bill Depth", "Kipp's Index", "Wing Cord", "Tail Length")
names(responsevariable.labs)<- c("BL.Average", "BW.Average", "BD.Average", "Kipp.s.Average", "WC.Average", "Tail")

p<-otu_cv_avg %>%
  pivot_longer(BL.Average:Tail, names_to="responsevariable", values_to = "variation") %>% 
  ggplot(aes(y=variation, x=Strategy)) +
  geom_boxplot() +
  facet_wrap(vars(responsevariable), ncol=3, labeller=labeller(responsevariable = responsevariable.labs)) +
  background_grid(major="none", minor="none") +
  panel_border() +
  labs(x="Migration Strategy", y= "Variation")

draw_plot(plot, x = 0, y = 0, width = 1, height = 1)
ggdraw() +  
  draw_plot(p,0,0,1,1) +
  draw_plot_label(c("G", "H", "I", "J", "K", "L"), c(0.075,0.385,0.69,0.075,0.385,0.69), c(0.99,0.99,0.99,0.53,0.53,0.53), size=15)
dev.off()

### Build a single figure from the two .pngs above
require(png)

png(file="Phenotype Residuals and Variation.png", width=7, height=5.5, units = "in", res=500)
rl<-lapply(list("./Output Files/Phenotype Residuals.png", "./Output Files/Phenotype Variation.png"), png::readPNG)
gl<-lapply(rl, grid::rasterGrob)
do.call(gridExtra::grid.arrange, gl)

dev.off()

### Build boxplots of PPCA scores
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

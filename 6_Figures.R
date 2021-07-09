rm(list=ls())

require(ggplot2)
require(tidyr)
require(cowplot)
require(gridExtra)
require(gridGraphics)

### Set source directory to the folder this file came from within RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Bring in morphology datasets with average morphology values for each OTU and PC scores
otu_avg<-read.csv('./Output Files/Tyrannus morphology + PCA avg_Adults.csv', row.names = 1)
otu_cv_avg<-read.csv('./Output Files/cv_summary_Adults.csv', row.names = 1) 
residuals<-read.csv('./Output Files/phylANOVA_tarsus-corrected_residuals.csv', row.names = 1)

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing sexes separately, use the following lines of code instead of the above 3 lines, then proceed.
otu_avg<-read.csv('./Output Files/Tyrannus morphology + PCA avg_Females_Adults.csv', row.names = 1)
otu_cv_avg<-read.csv('./Output Files/cv_summary_Females_Adults.csv', row.names = 1) 
residuals<-read.csv('./Output Files/phylANOVA_tarsus-corrected_residuals_Females.csv', row.names = 1)
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

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing sexes separately, use the following lines of code instead of the above 3 lines, then proceed.
png(file="Phenotype Residuals_Females.png",width=7,height=5.5,units="in",res=500)
### OR
#png(file="Phenotype Residuals_Males.png",width=7,height=5.5,units="in",res=500)

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

responsevariable.labs<-c("Bill Length", "Bill Width", "Bill Depth", "Kipp's Distance", "Wing Length", "Tail Length")
names(responsevariable.labs)<- c("BL", "BW", "BD", "KD", "WC", "TL")

col<-c("black","gray40","gray70","black","gray40","gray70","black","gray40","gray70","black","gray40","gray70","black","gray40","gray70","black","gray40","gray70")

p<-residuals %>%
  pivot_longer(BL:TL, names_to = "responsevariable", values_to = "residuals") %>%
  ggplot(aes(y=residuals, x=migration)) + 
  geom_boxplot(color=col) +
  facet_wrap(vars(responsevariable), ncol=3, labeller=labeller(responsevariable = responsevariable.labs),scales="free_y") +
  panel_border() +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text=element_text(size=15)) +
  scale_y_continuous(expand=c(0.02,10)) +
  labs(y= "Residuals") 

draw_plot(plot, x = 0, y = 0, width = 1, height = 1)
ggdraw() +  
  draw_plot(p,0,0,1,1) +
  draw_plot_label(c("A", "B", "C", "D", "E", "F"), c(0.075,0.395,0.715,0.075,0.395,0.715), c(0.99,0.99,0.99,0.50,0.50,0.50), size=15) +
  draw_plot_label(c("A", "A", "A"), c(0.115, 0.20, 0.285), c(0.935, 0.935, 0.935), size=10) + #label within plot A
  draw_plot_label(c("A", "A", "A"), c(0.43, 0.522, 0.607), c(0.935, 0.935, 0.935), size=10) + #label within plot B
  draw_plot_label(c("A", "A", "A"), c(0.755, 0.84, 0.926), c(0.935, 0.935, 0.935), size=10) + #label within plot C
  draw_plot_label(c("A", "AB", "B"), c(0.115,0.19, 0.285), c(0.445, 0.445, 0.445), size=10) + #label within plot D
  draw_plot_label(c("A", "A", "A"), c(0.43, 0.522, 0.607), c(0.445, 0.445, 0.445), size=10) + #label within plot E
  draw_plot_label(c("A", "AB", "B"), c(0.755, 0.83, 0.926), c(0.445, 0.445, 0.445), size=10) #label within plot F

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing sexes separately, use the following lines of code instead of the above 3 lines, then proceed.
ggdraw() +  
  draw_plot(p,0,0,1,1) +
  draw_plot_label(c("A", "B", "C", "D", "E", "F"), c(0.075,0.395,0.715,0.075,0.395,0.715), c(0.99,0.99,0.99,0.50,0.50,0.50), size=15) +
  draw_plot_label(c("A", "A", "A"), c(0.115, 0.20, 0.285), c(0.935, 0.935, 0.935), size=8) + #label within plot A
  draw_plot_label(c("A", "A", "A"), c(0.43, 0.522, 0.607), c(0.935, 0.935, 0.935), size=8) + #label within plot B
  draw_plot_label(c("A", "A", "A"), c(0.755, 0.84, 0.926), c(0.935, 0.935, 0.935), size=8) + #label within plot C
  draw_plot_label(c("A", "AB", "B"), c(0.115,0.20, 0.285), c(0.445, 0.445, 0.445), size=8) + #label within plot D
  draw_plot_label(c("A", "A", "A"), c(0.43, 0.522, 0.607), c(0.445, 0.445, 0.445), size=8) + #label within plot E
  draw_plot_label(c("A", "AB", "B"), c(0.755, 0.83, 0.926), c(0.445, 0.445, 0.445), size=8) #label within plot F

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###
  
dev.off()

### Make the box plots of variation in morphometrics
png(file="Phenotype Variation.png",width=7,height=5.5,units="in",res=500)

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing sexes separately, use the following lines of code instead of the above 3 lines, then proceed.
png(file="Phenotype Variation_Females.png",width=7,height=5.5,units="in",res=500)
### OR
#png(file="Phenotype Variation_Males.png",width=7,height=5.5,units="in",res=500)

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

otu_cv_avg<-subset(otu_cv_avg, select=-c(Tarsus.Average)) #Remove tarsus from the df to be able to call BL:TL below

col<-c("black","gray40","gray70","black","gray40","gray70","black","gray40","gray70","black","gray40","gray70","black","gray40","gray70","black","gray40","gray70")

responsevariable.labs<-c("Bill Length", "Bill Width", "Bill Depth", "Kipp's Distance", "Wing Length", "Tail Length")
names(responsevariable.labs)<- c("BL.Average", "BW.Average", "BD.Average", "Kipp.s.Distance", "WC.Average", "Tail")

p<-otu_cv_avg %>%
  pivot_longer(BL.Average:Tail, names_to="responsevariable", values_to = "variation") %>% 
  ggplot(aes(y=variation, x=Strategy, scale_fill_manual(values=c("gray70","gray40","black")))) +
  geom_boxplot(color=col) +
  facet_wrap(vars(responsevariable), ncol=3, labeller=labeller(responsevariable = responsevariable.labs)) +
  panel_border() +
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        strip.text=element_text(size=15)) +
  scale_y_continuous(expand=c(0,0.03)) +
  labs(x="Migration Strategy", y= "Variation")

draw_plot(plot, x = 0, y = 0, width = 1, height = 1)
ggdraw() +  
  draw_plot(p,0,0,1,1) +
  draw_plot_label(c("G", "H", "I", "J", "K", "L"), c(0.075,0.385,0.69,0.075,0.385,0.69), c(0.99,0.99,0.99,0.53,0.53,0.53), size=15) +
  draw_plot_label(c("A", "A", "A"), c(0.12, 0.22, 0.31), c(0.935, 0.935, 0.935), size=10) + #label within plot G
  draw_plot_label(c("A", "A", "A"), c(0.4262, 0.522, 0.614), c(0.935, 0.935, 0.935), size=10) + #label within plot H
  draw_plot_label(c("A", "A", "A"), c(0.73, 0.826, 0.924), c(0.935, 0.935, 0.935), size=10) + #label within plot I
  draw_plot_label(c("A", "A", "A"), c(0.12,0.22, 0.31), c(0.475, 0.475, 0.475), size=10) + #label within plot J
  draw_plot_label(c("A", "A", "A"), c(0.4262, 0.522, 0.614), c(0.475, 0.475, 0.475), size=10) + #label within plot K
  draw_plot_label(c("A", "A", "A"), c(0.73, 0.826, 0.924), c(0.475, 0.475, 0.475), size=10) #label within plot L

dev.off()

### Build a single figure from the two .pngs above
require(png)

png(file="Phenotype Residuals and Variation.png", width=7, height=5.5, units = "in", res=500)

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing sexes separately, use the following lines of code instead of the above 3 lines, then proceed.
png(file="Phenotype Residuals and Variation_Females.png", width=7, height=5.5, units = "in", res=500)
### OR
#png(file="Phenotype Residuals and Variation_Males.png", width=7, height=5.5, units = "in", res=500)

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

rl<-lapply(list("./Output Files/Phenotype Residuals.png", "./Output Files/Phenotype Variation.png"), png::readPNG)

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing sexes separately, use the following lines of code instead of the above 3 lines, then proceed.
rl<-lapply(list("./Output Files/Phenotype Residuals_Females.png", "./Output Files/Phenotype Variation_Females.png"), png::readPNG)
### OR
#rl<-lapply(list("./Output Files/Phenotype Residuals_Males.png", "./Output Files/Phenotype Variation_Males.png"), png::readPNG)

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

gl<-lapply(rl, grid::rasterGrob)
do.call(gridExtra::grid.arrange, gl)

dev.off()

### Build boxplots of PPCA scores
png(file="PPC scores boxplots.png",width=7,height=5.5,units="in",res=500)

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing sexes separately, use the following lines of code instead of the above 3 lines, then proceed.
png(file="PPC scores boxplots_Females.png",width=7,height=2.75,units="in",res=500)
### OR
#png(file="PPC scores boxplots_Males.png",width=7,height=5.5,units="in",res=500)

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

responsevariable.labs<-c("Bill PPC1", "Bill PPC2")
names(responsevariable.labs)<- c("BillPC1", "BillPC2")

col<-c("black","gray40","gray70","black","gray40","gray70")

p<-otu_avg %>%
  pivot_longer(BillPC1:BillPC2, names_to="responsevariable", values_to = "value") %>% 
  ggplot(aes(y=value, x=Strategy)) +
  geom_boxplot(color=col) +
  facet_wrap(vars(responsevariable), ncol=2, labeller=labeller(responsevariable = responsevariable.labs), scales="free_y") +
  panel_border() +
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        strip.text=element_text(size=15)) +
  scale_y_continuous(expand=c(0.05,8)) +
  labs(x="Migration Strategy", y= "PPC Scores")

draw_plot(plot, x = 0, y = 0, width = 1, height = 1)
ggdraw() +  
  draw_plot(p,0,0,1,1) +
  draw_plot_label(c("A", "B"), c(0.085,0.565), c(0.99,0.99), size=15) + #main plot identifiers 
  draw_plot_label(c("A", "A", "A"), c(0.146, 0.287, 0.420), c(0.87, 0.87, 0.87), size=10) + #letters for plot A
  draw_plot_label(c("A", "AB", "B"), c(0.625, 0.746, 0.9), c(0.87, 0.87, 0.87), size=10) #letters for plot B


###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing females separately, use the following lines of code instead of the above 3 lines, then proceed.
#ggdraw() +  
#  draw_plot(p,0,0,1,1) +
#  draw_plot_label(c("A", "B", "C", "D"), c(0.085,0.565,0.085,0.565), c(0.99,0.99,0.53,0.53), size=15) + #main plot identifiers 
#  draw_plot_label(c("A", "A", "A"), c(0.156, 0.287, 0.420), c(0.946, 0.946, 0.946), size=8) + #letters for plot A
#  draw_plot_label(c("A", "AB", "B"), c(0.635, 0.766, 0.9), c(0.946, 0.946, 0.946), size=8) + #letters for plot B
#  draw_plot_label(c("A", "A", "A"), c(0.156, 0.287, 0.420), c(0.487, 0.487, 0.487), size=8) + #letters for plot C
#  draw_plot_label(c("A", "A", "A"), c(0.635, 0.766, 0.9), c(0.487, 0.487, 0.487), size=8) #letters for plot D

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

dev.off()

### Build ancestral state reconstruction for wing chord length for supplementary material
require(phytools)
png(file="WC_ASR.png",width=7,height=5.5,units="in",res=500)
png(file="WC_ASR_Females.png",width=7,height=5.5,units="in",res=500)


### Read in phylogeny 
phy<-read.tree('./Output Files/Tyrannus_phylogeny.tre')

WC<-residuals$WC
rownames(residuals)<-residuals$Row.names
names(WC)<- rownames(residuals)

obj<-contMap(phy, WC)
n<-length(obj$cols) #get the length of the color ramp
obj$cols[1:n]<-grey(0:(n-1)/(n-1)) #change the color ramp to a grey scale

plot(obj, ftype="on", mar=c(5,5,5,5), legend=0.7*max(nodeHeights(phy), fsize=c(0.7,0.9)))

format_tl<-phy$tip.label
format_tl<-gsub("_"," ",format_tl)
format_tl<-gsub("Tyrannus","T.",format_tl)

text(x=rep(5.2,Ntip(phy)),y=1:Ntip(phy),label=format_tl,adj=c(0,0.5),font=3)

dev.off()

### Build ancestral state reconstruction for Kipp's distance for supplementary material
png(file="KD_ASR.png",width=7,height=5.5,units="in",res=500)

KD<-residuals$KD
names(KD)<- rownames(residuals)

obj2<-contMap(phy, KD)
n<-length(obj2$cols) #get the length of the color ramp
obj2$cols[1:n]<-grey(0:(n-1)/(n-1)) #change the color ramp to a grey scale

plot(obj2, ftype="on", mar=c(5,5,5,5), legend=0.7*max(nodeHeights(phy), fsize=c(0.7,0.9)))

format_tl<-phy$tip.label
format_tl<-gsub("_"," ",format_tl)
format_tl<-gsub("Tyrannus","T.",format_tl)

text(x=rep(5.2,Ntip(phy)),y=1:Ntip(phy),label=format_tl,adj=c(0,0.5),font=3)

dev.off()

### Build ancestral state reconstruction for Bill PC2 for supplementary material
png(file="BillPPC2_ASR.png",width=7,height=5.5,units="in",res=500)

BillPC2<-otu_avg$BillPC2
names(BillPC2)<- rownames(otu_avg)

obj2<-contMap(phy, BillPC2)
n<-length(obj2$cols) #get the length of the color ramp
obj2$cols[1:n]<-grey(0:(n-1)/(n-1)) #change the color ramp to a grey scale

plot(obj2, ftype="on", mar=c(5,5,5,5), legend=0.7*max(nodeHeights(phy), fsize=c(0.7,0.9)))

format_tl<-phy$tip.label
format_tl<-gsub("_"," ",format_tl)
format_tl<-gsub("Tyrannus","T.",format_tl)

text(x=rep(5.2,Ntip(phy)),y=1:Ntip(phy),label=format_tl,adj=c(0,0.5),font=3)

dev.off()

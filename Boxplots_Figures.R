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
  ggplot(aes(y=residuals, x=migration, scale_fill_manual(values=c("gray70","gray40","black")))) + #the greys are not appearing
  geom_boxplot() +
  facet_wrap(vars(responsevariable), ncol=3, labeller=labeller(responsevariable = responsevariable.labs)) +
  panel_border() +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y= "Residuals") 

draw_plot(plot, x = 0, y = 0, width = 1, height = 1)
ggdraw() +  
  draw_plot(p,0,0,1,1) +
  draw_plot_label(c("A", "B", "C", "D", "E", "F"), c(0.075,0.385,0.69,0.075,0.385,0.69), c(0.99,0.99,0.99,0.50,0.50,0.50), size=15) +
  draw_plot_label(c("A", "A", "A"), c(0.12, 0.22, 0.31), c(0.945, 0.945, 0.945), size=7) + #label within plot A
  draw_plot_label(c("A", "A", "A"), c(0.43, 0.522, 0.617), c(0.945, 0.945, 0.945), size=7) + #label within plot B
  draw_plot_label(c("A", "A", "A"), c(0.73, 0.826, 0.924), c(0.945, 0.945, 0.945), size=7) + #label within plot C
  draw_plot_label(c("A", "A", "A"), c(0.12,0.22, 0.31), c(0.455, 0.455, 0.455), size=7) + #label within plot D
  draw_plot_label(c("A", "A", "A"), c(0.43, 0.522, 0.617), c(0.455, 0.455, 0.455), size=7) + #label within plot E
  draw_plot_label(c("A", "A", "A"), c(0.73, 0.826, 0.924), c(0.455, 0.455, 0.455), size=7) #label within plot F
  
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
  draw_plot_label(c("G", "H", "I", "J", "K", "L"), c(0.075,0.385,0.69,0.075,0.385,0.69), c(0.99,0.99,0.99,0.53,0.53,0.53), size=15) +
  draw_plot_label(c("A", "A", "A"), c(0.12, 0.22, 0.31), c(0.945, 0.945, 0.945), size=7) + #label within plot G
  draw_plot_label(c("A", "A", "A"), c(0.4262, 0.522, 0.614), c(0.945, 0.945, 0.945), size=7) + #label within plot H
  draw_plot_label(c("A", "A", "A"), c(0.73, 0.826, 0.924), c(0.945, 0.945, 0.945), size=7) + #label within plot I
  draw_plot_label(c("A", "A", "A"), c(0.12,0.22, 0.31), c(0.487, 0.487, 0.487), size=7) + #label within plot J
  draw_plot_label(c("A", "A", "A"), c(0.4262, 0.522, 0.614), c(0.487, 0.487, 0.487), size=7) + #label within plot K
  draw_plot_label(c("A", "A", "A"), c(0.73, 0.826, 0.924), c(0.487, 0.487, 0.487), size=7) #label within plot L

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

responsevariable.labs<-c("Bill PPC1", "Bill PPC2", "Feather PPC1", "Feather PPC2")
names(responsevariable.labs)<- c("BillPC1", "BillPC2", "BodyPC1", "BodyPC2")

p<-otu_avg %>%
  pivot_longer(BillPC1:BodyPC2, names_to="responsevariable", values_to = "value") %>% 
  ggplot(aes(y=value, x=Strategy)) +
  geom_boxplot() +
  facet_wrap(vars(responsevariable), ncol=2, labeller=labeller(responsevariable = responsevariable.labs), scales="free") +
  background_grid(major="none", minor="none") +
  panel_border() +
  labs(x="Migration Strategy", y= "PPC Scores")

draw_plot(plot, x = 0, y = 0, width = 1, height = 1)
ggdraw() +  
  draw_plot(p,0,0,1,1) +
  draw_plot_label(c("A", "B", "C", "D"), c(0.085,0.565,0.085,0.565), c(0.99,0.99,0.51,0.51), size=15) + #main plot identifiers 
  draw_plot_label(c("A", "A", "A"), c(0.156, 0.29, 0.422), c(0.946, 0.946, 0.946), size=7) + #letters for plot A
  draw_plot_label(c("A", "AB", "B"), c(0.635, 0.766, 0.9), c(0.946, 0.946, 0.946), size=7) + #letters for plot B
  draw_plot_label(c("A", "A", "A"), c(0.156, 0.29, 0.422), c(0.4715, 0.4715, 0.4715), size=7) + #letters for plot C
  draw_plot_label(c("A", "A", "A"), c(0.635, 0.766, 0.9), c(0.4715, 0.4715, 0.4715), size=7) #letters for plot D

dev.off()

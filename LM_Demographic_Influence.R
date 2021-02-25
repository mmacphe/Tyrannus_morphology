rm(list=ls())
getwd()

>library(readx1)
Tyrannus_voucher_table<-read_excel("C:/Users/maggi/Downloads/Tyrannus voucher table (1).xlsx")
Morph_data<-Tyrannus_voucher_table_uk_removed

#What we're trying to do is to see if:
#1) morphologies differ between the sexes (in each taxon?)
#2) morphologies differ between age classes (in each taxon?)
#3) morphologies vary with body size (in each taxon)

#First step is to classify each taxon 
Morph_data$subspecies<-trimws(Morph_data$subspecies)
Morph_data$trinomial<-paste(Morph_data$Species,Morph_data$subspecies)
Morph_data$trinomial<-trimws(Morph_data$trinomial)
table(Morph_data$trinomial) #unknown subspecies were removed

#/\/\/\/\/# Q1) Do the morphologies differ between the sexes? #unknown sex were removed
Morph_data<-Morph_data[!Morph_data$Sex=="",]
Morph_data$Sex[Morph_data$Sex=="H"]<-"Female"
Morph_data$Sex[Morph_data$Sex=="M"]<-"Male"
Morph_data$Sex[Morph_data$Sex=="D"]<-"Female"
Morph_data$Sex[Morph_data$Sex=="F"]<-"Female"
Morph_data$Sex[Morph_data$Sex=="Female "]<-"Female"
Morph_data$Sex[Morph_data$Sex=="NA"]<-"Unknown"
Morph_data$Sex[Morph_data$Sex=="unknown"]<-"Unknown"
table(Morph_data$Sex)

lm1<-lm(Morph_data$`BL Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`BL Average`~ Morph_data$trinomial*Morph_data$Sex)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 4072.592
#lm2 57 4086.302

lm1<-lm(Morph_data$`BW Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`BW Average`~ Morph_data$trinomial*Morph_data$Sex)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 2409.433
#lm2 57 2415.325

lm1<-lm(Morph_data$`BD Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`BD Average`~ Morph_data$trinomial*Morph_data$Sex)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 1483.731
#lm2 57 1478.459

lm1<-lm(Morph_data$`Kipp's Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`Kipp's Average`~ Morph_data$trinomial*Morph_data$Sex)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 5903.397
#lm2 57 5528.015

lm1<-lm(Morph_data$`WC Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`WC Average`~ Morph_data$trinomial*Morph_data$Sex)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 7943.603
#lm2 57 7813.485

lm1<-lm(Morph_data$Tail ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$Tail ~ Morph_data$trinomial*Morph_data$Sex)

AIC(lm1, lm2)
#    df       AIC
#lm1 29 10195.159
#lm2 57  9657.571

#/\/\/\/\/#Q2) Do the morphologies differ between age classes?
Morph_data<-Morph_data[!Morph_data$Age=="",]
Morph_data$Age[Morph_data$Age=="Adult "]<-"Adult"
Morph_data$Age[Morph_data$Age=="AD"]<-"Adult"
Morph_data$Age[Morph_data$Age=="SY"]<-"Adult"
Morph_data$Age[Morph_data$Age=="nestling"]<-"Nestling"
Morph_data$Age[Morph_data$Age=="JV"]<-"Juvenile"
Morph_data$Age[Morph_data$Age=="HY"]<-"Juvenile"
Morph_data<-Morph_data[!Morph_data$Age=="Nestling",]
Morph_data<-Morph_data[!Morph_data$Age=="Fledgling",]
names(table(Morph_data$Age)) #unknown + Juvenile (n=267) removed
table(Morph_data$Age)
#There is only one age category (Adult) so LM1 should always be better than LM2

lm1<-lm(Morph_data$`BL Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`BL Average`~ Morph_data$trinomial*Morph_data$Age)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 4072.592
#lm2 57 9657.571

lm1<-lm(Morph_data$`BW Average`~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`BW Average`~ Morph_data$trinomial*Morph_data$Age)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 2409.433
#lm2 57 9657.571

lm1<-lm(Morph_data$`BD Average`~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`BD Average`~ Morph_data$trinomial*Morph_data$Age)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 1483.731
#lm2 57 9657.571

lm1<-lm(Morph_data$`Kipp's Average`~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`Kipp's Average`~ Morph_data$trinomial*Morph_data$Age)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 5903.397
#lm2 57 9657.571

lm1<-lm(Morph_data$`WC Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`WC Average`~ Morph_data$trinomial*Morph_data$Age)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 7943.603
#lm2 57 9657.571

lm1<-lm(Morph_data$Tail ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$Tail ~ Morph_data$trinomial*Morph_data$Age)

AIC(lm1, lm2)
#    df       AIC
#lm1 29 10195.159
#lm2 57  9657.571

#/\/\/\/\/#Q3) Do the morphologies differ with tarsus length?
lm1<-lm(Morph_data$`BL Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`BL Average`~ Morph_data$trinomial*Morph_data$`Tarsus Average`)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 4072.592
#lm2 57 3839.143

lm1<-lm(Morph_data$`BW Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`BW Average`~ Morph_data$trinomial*Morph_data$`Tarsus Average`)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 2409.433
#lm2 57 2218.043

lm1<-lm(Morph_data$`BD Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`BD Average`~ Morph_data$trinomial*Morph_data$`Tarsus Average`)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 1483.731
#lm2 57 1337.364

lm1<-lm(Morph_data$`Kipp's Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`Kipp's Average`~ Morph_data$trinomial*Morph_data$`Tarsus Average`)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 5903.397
#lm2 57 5896.600

lm1<-lm(Morph_data$`WC Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`WC Average`~ Morph_data$trinomial*Morph_data$`Tarsus Average`)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 7943.603
#lm2 57 7366.929

lm1<-lm(Morph_data$Tail ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$Tail ~ Morph_data$trinomial*Morph_data$`Tarsus Average`)

AIC(lm1, lm2)
#    df       AIC
#lm1 29 10195.159
#lm2 57  9557.213
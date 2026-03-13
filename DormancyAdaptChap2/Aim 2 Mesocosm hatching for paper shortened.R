#Aim 2 Mesocosm hatching


library(data.table)
library(ggplot2)
library(foreach)
library(lattice)
library(tidyr)
library(tidyverse)
library(gridExtra)
library(dplyr)
library(reshape2)
library(patchwork)
library(knitr)
library(lme4)
library(cowplot)
library(ggsignif)


#Data sheets uploaded as of 3/17/25

hatching_data <- fread("C:/Users/rjpor/Downloads/Aim1 Mesocosm Hatching.csv")
hatching_data <- subset(hatching_data, Plate != "")
Mesolist <- fread("C:/Users/rjpor/Downloads/Mesolist20240306.csv")
sites <- read.csv("C:\\Users\\rjpor\\Downloads\\Sites.csv")

#from Genoseqlocal.R
merged_species <- read.csv("C:/Users/rjpor/Downloads/filtered_species0.20.csv")

wellxclone <- read.csv("C:\\Users\\rjpor\\Downloads\\wellxclone.csv")
#
hatching_data$plate_clone <- paste(hatching_data$Pond, hatching_data$Clone)


hatching_data$plate_clone <- gsub(" ", "_", hatching_data$plate_clone)
unique(hatching_data$plate_clone)

hatching_data2 <- hatching_data

#this was for putative species

hatching_data2$Condition <- str_replace(hatching_data2$Condition, "Control", "Cont")
hatching_data2$Condition <- str_replace(hatching_data2$Condition, "Cont", "Control")
hatching_data2$Condition <- str_replace(hatching_data2$Condition, "Dry", "Warm")


unique(hatching_data$Condition)
Conditiondrieddry <- subset(hatching_data2, Condition == "Warm" & Dried == "Y")

cutoff_date_before <- as.Date("6/7", format = "%m/%d")
cutoff_date_after <- as.Date("6/7", format = "%m/%d")

# Update the dates column based on the conditions
hatching_data2$Date.Collected <- as.Date(hatching_data2$Date.Collected, tryFormats = c("%m/%d", "%m/%d/%Y", "%d-%m-%Y","%d-%b"))

NAS2 <- subset(hatching_data2, is.na(hatching_data2$Date.Collected))

hatching_data2$Date.Collected[hatching_data2$Date.Collected < cutoff_date_before] <- as.Date(format(hatching_data2$Date.Collected[hatching_data2$Date.Collected < cutoff_date_before], "%m/%d/2024"), "%m/%d/%Y")
hatching_data2$Date.Collected[hatching_data2$Date.Collected > cutoff_date_after] <- as.Date(format(hatching_data2$Date.Collected[hatching_data2$Date.Collected > cutoff_date_after], "%m/%d/2023"), "%m/%d/%Y")

hatching_data2$Hatch1 <- as.Date(hatching_data2$Hatch1, , tryFormats = c("%m/%d", "%m/%d/%Y", "%d-%m-%Y","%d-%b"))
hatching_data2$Hatch2 <- as.Date(hatching_data2$Hatch2, , tryFormats = c("%m/%d", "%m/%d/%Y", "%d-%m-%Y","%d-%b"))
hatching_data2$Hatch3 <- as.Date(hatching_data2$Hatch3, , tryFormats = c("%m/%d", "%m/%d/%Y", "%d-%m-%Y","%d-%b"))
hatching_data2$Hatch4 <- as.Date(hatching_data2$Hatch4, , tryFormats = c("%m/%d", "%m/%d/%Y", "%d-%m-%Y","%d-%b"))

hatching_data2$Hatch1[hatching_data2$Hatch1 < cutoff_date_before & !is.na(hatching_data2$Hatch1)] <- as.Date(format(hatching_data2$Hatch1[hatching_data2$Hatch1 < cutoff_date_before & !is.na(hatching_data2$Hatch1)], "%m/%d/2024"), "%m/%d/%Y")
hatching_data2$Hatch1[hatching_data2$Hatch1 > cutoff_date_after & !is.na(hatching_data2$Hatch1)] <- as.Date(format(hatching_data2$Hatch1[hatching_data2$Hatch1 > cutoff_date_after & !is.na(hatching_data2$Hatch1)], "%m/%d/2023"), "%m/%d/%Y")

hatching_data2$Hatch2[hatching_data2$Hatch2 < cutoff_date_before & !is.na(hatching_data2$Hatch2)] <- as.Date(format(hatching_data2$Hatch2[hatching_data2$Hatch2 < cutoff_date_before & !is.na(hatching_data2$Hatch2)], "%m/%d/2024"), "%m/%d/%Y")
hatching_data2$Hatch2[hatching_data2$Hatch2 > cutoff_date_after & !is.na(hatching_data2$Hatch2)] <- as.Date(format(hatching_data2$Hatch2[hatching_data2$Hatch2 > cutoff_date_after & !is.na(hatching_data2$Hatch2)], "%m/%d/2023"), "%m/%d/%Y")

hatching_data2$Hatch3[hatching_data2$Hatch3 < cutoff_date_before & !is.na(hatching_data2$Hatch3)] <- as.Date(format(hatching_data2$Hatch3[hatching_data2$Hatch3 < cutoff_date_before & !is.na(hatching_data2$Hatch3)], "%m/%d/2024"), "%m/%d/%Y")
hatching_data2$Hatch3[hatching_data2$Hatch3 > cutoff_date_after & !is.na(hatching_data2$Hatch3)] <- as.Date(format(hatching_data2$Hatch3[hatching_data2$Hatch3 > cutoff_date_after & !is.na(hatching_data2$Hatch3)], "%m/%d/2023"), "%m/%d/%Y")

hatching_data2$Hatch4[hatching_data2$Hatch4 < cutoff_date_before & !is.na(hatching_data2$Hatch4)] <- as.Date(format(hatching_data2$Hatch4[hatching_data2$Hatch4 < cutoff_date_before & !is.na(hatching_data2$Hatch4)], "%m/%d/2024"), "%m/%d/%Y")
hatching_data2$Hatch4[hatching_data2$Hatch4 > cutoff_date_after & !is.na(hatching_data2$Hatch4)] <- as.Date(format(hatching_data2$Hatch4[hatching_data2$Hatch4 > cutoff_date_after & !is.na(hatching_data2$Hatch4)], "%m/%d/2023"), "%m/%d/%Y")

unique(hatching_data2$Dried)

#total ephipia tracked
hatching_data2<- subset(hatching_data2, Dried == "Y"| Dried == "N")

Mesolist <- Mesolist[,c(1,2,3,8)]
colnames(Mesolist)[3] ="Replicate"
Mesolist$Clone <- as.character(Mesolist$Clone)
Mesolist$plate_clone <- paste(Mesolist$Pond, Mesolist$Clone)
Mesolist$plate_clone <- gsub(" ", "_", Mesolist$plate_clone)

setDT(hatching_data2)
setDT(Mesolist)
setkey(hatching_data2, Pond, Clone, Replicate, plate_clone)
setkey(Mesolist, Pond, Clone, Replicate, plate_clone)
hatching_data2plusnoepp <- full_join(hatching_data2, Mesolist, by=c("plate_clone","Replicate","Pond","Clone"))

unique(hatching_data2plusnoepp$Plate)
unique(hatching_data2plusnoepp$plate_clone)

hatching_data2plusnoepp$Plate[is.na(hatching_data2plusnoepp$Plate)] <- 1000
subsetnoepp <- subset(hatching_data2plusnoepp, Plate == 1000)
hatching_data3 <- subset(hatching_data2plusnoepp, Plate != 1000)

NumofEpp <- hatching_data3[, .N, by=list(Pond, Clone, Replicate, plate_clone, Species)]

#total number of ephippia
sum(NumofEpp$N)                         

#Add in the mesocosms that didnt produce any ephippia
noeppproduced <- subset(hatching_data2plusnoepp, Plate == 1000)
noeppproduced <- noeppproduced[,c(1,3,4,5,22,23)]

noeppproducedR1 <- noeppproduced
noeppproducedR1$Replicate <- "R1"
noeppproducedR2 <- noeppproduced
noeppproducedR2$Replicate <- "R2"

noeppproducedall <- rbind(noeppproducedR1, noeppproducedR2)

noeppproducedall$total <- 0
noeppproducedall$N <- 0
noeppproducedall$prop <- 0

NumofEpp2 <- rbind(NumofEpp,noeppproducedall, fill = TRUE)

NumofEpp2 <- as.data.table(NumofEpp2)
NumofEpp2 <- subset(NumofEpp2, Species!= "NA")

noeppproduced <- subset(hatching_data2plusnoepp, Plate == 1000)
noeppproduced <- noeppproduced[,c(1,3,4,5,22,23)]

noeppproducedwarm <- noeppproduced
noeppproducedwarm$Condition <- "Warm"
noeppproducedwarm$Dried <- "Y"

noeppproducedwet <- noeppproduced
noeppproducedwet$Condition <- "Cold"
noeppproducedwet$Dried <- "N"

noeppproducedcontrol <- noeppproduced
noeppproducedcontrol$Condition <- "Control"
noeppproducedcontrol$Dried <- "N"

noeppproducedwetdry <- noeppproduced
noeppproducedwetdry$Condition <- "Cold"
noeppproducedwetdry$Dried <- "Y"

noeppproducedall2 <- rbind(noeppproducedwarm,noeppproducedwet,noeppproducedcontrol,noeppproducedwetdry)
noeppproducedall2$total <- 0
noeppproducedall2$N <- 0
noeppproducedall2$prop <- 0

hatching_data2plusnoepp2 <- rbind(hatching_data3, noeppproducedall2, fill =TRUE)

hatching_data3h1 <- hatching_data3[,c(1,2,5,6,7,8,15,16,17,20,22,23)]
hatching_data3h1$sib <- "a"
hatching_data3h2 <- hatching_data3[,c(1,2,5,6,9,10,15,16,17,20,22,23)]
hatching_data3h2$sib <- "b"
hatching_data3h3 <- hatching_data3[,c(1,2,5,6,11,12,15,16,17,20,22,23)]
hatching_data3h3$sib <- "c"
hatching_data3h4 <- hatching_data3[,c(1,2,5,6,13,14,15,16,17,20,22,23)]
hatching_data3h4$sib <- "d"

colnames(hatching_data3h2)[5] ="Hatch1"
colnames(hatching_data3h2)[6] ="Hatch1HD"
colnames(hatching_data3h3)[5] ="Hatch1"
colnames(hatching_data3h3)[6] ="Hatch1HD"
colnames(hatching_data3h4)[5] ="Hatch1"
colnames(hatching_data3h4)[6] ="Hatch1HD"
hatching_data_long <- rbind(hatching_data3h1,hatching_data3h2)
#hatching_data_long <- rbind(hatching_data3h1,hatching_data3h2,hatching_data3h3,hatching_data3h4)
hatching_data_long$h1 <- hatching_data_long$Hatch1- hatching_data_long$Date.Collected 

hatching_data_long <- hatching_data_long %>% separate(h1, c('Days', 'Day'), sep = " ")
hatching_data_long$Days[is.na(hatching_data_long$Days)] <- -100

hatching_data_long <- as.data.table(hatching_data_long)


hatching_data_long$Days<- as.numeric(hatching_data_long$Days)
Early <- subset(hatching_data_long, Days <= 45 & Days >=0)
Early$time <- "early"
Early2 <- Early[, .N, by=list(plate_clone, Replicate, Species, time)]

Later<- subset(hatching_data_long, Days > 45)
Later$time <- "later"
Later2 <- Later[, .N, by=list(plate_clone, Replicate, Species, time)]

Nohatch<- subset(hatching_data_long, Days < 0)
Nohatch$time <- "Nohatch"

hatchings <- rbind(Early, Later, Nohatch)
hatchshort <- hatchings[, .N, by=list(plate_clone, time, Replicate)]

hatchshortsmol <- hatchings[, .N, by=list(plate_clone, Replicate)]
hatchshortsmol
colnames(hatchshortsmol)[3] ="total"

hatchshortpluscond <- hatchings[, .N, by=list(plate_clone, time, Condition, Replicate, Dried, Species)]

setDT(hatchshort)
setDT(hatchshortsmol)
setkey(hatchshort, plate_clone, Replicate)
setkey(hatchshortsmol, plate_clone, Replicate)
hatchshortmerge <- merge(hatchshort, hatchshortsmol, all.x=TRUE, all.y=TRUE)

hatchshortmerge$prop <- hatchshortmerge$N/hatchshortmerge$total

hatchshortmerge$se <- sqrt(hatchshortmerge$prop * (1-hatchshortmerge$prop)/ (hatchshortmerge$total))
hatchshortmerge$lci <- hatchshortmerge$prop - (1.96 * hatchshortmerge$se)
hatchshortmerge$uci <- hatchshortmerge$prop + (1.96 * hatchshortmerge$se)
hatchshortonlyhatch <- subset(hatchshort, time =="later"| time== "early")

hatchshort <- hatchings[, .N, by=list(plate_clone, time, Replicate, Condition, Dried)]
hatchshort2 <- subset(hatchshort, Condition =="Cold" & Dried == "N")
hatchshort2early <- subset(hatchshort2, time =="early")

hatchingscoldwet <- subset(hatchings, Condition =="Cold" & Dried== "N")

hatchshortsmol2 <- subset(hatchingscoldwet, time =="later"| time== "early")[, .N, by=list(plate_clone, Replicate,Species)]
hatchshortsmol2
colnames(hatchshortsmol2)[4] ="total"

setDT(hatchshortonlyhatch)
setDT(hatchshortsmol2)
setkey(hatchshortonlyhatch, plate_clone, Replicate)
setkey(hatchshortsmol2, plate_clone, Replicate)
hatchshortmergeshort <- merge(hatchshortonlyhatch, hatchshortsmol2, all.x=TRUE, all.y=TRUE)

hatchshortmergeshort$prop <- hatchshortmergeshort$N/hatchshortmergeshort$total

colnames(hatchshortmergeshort)[3] ="Hatching_time"

hatchshortmergeshort$se <- sqrt(hatchshortmergeshort$prop * (1-hatchshortmergeshort$prop)/ (hatchshortmergeshort$total))
hatchshortmergeshort$lci <- hatchshortmergeshort$prop - (1.96 * hatchshortmergeshort$se)
hatchshortmergeshort$uci <- hatchshortmergeshort$prop + (1.96 * hatchshortmergeshort$se)

#late hatching

ColdWet <- subset(Later, Condition == "Cold" & Dried == "N")[, .N, by=list(plate_clone, Replicate, Condition, Dried, Species)]
ColdDry <- subset(Later, Condition == "Cold" & Dried == "Y")[, .N, by=list(plate_clone, Replicate, Condition, Dried, Species)]
WarmDry  <- subset(Later, Condition == "Warm" & Dried == "Y")[, .N, by=list(plate_clone, Replicate, Condition, Dried, Species)]
Control <- subset(Later, Condition == "Control")[, .N, by=list(plate_clone, Replicate, Condition, Dried, Species)]


all_cond <- rbind(ColdWet,ColdDry,WarmDry,Control)

Dryepps <- subset(hatching_data3, Condition == "Warm" & Dried == "Y")[, .N, by=list(plate_clone, Replicate, Condition, Dried)]
Onlydrys <- subset(hatching_data3, Condition == "Warm" & Dried == "Y")

hatching_data2only120 <- subset(Onlydrys, Leftover_eggs == 0| Leftover_eggs == 1| Leftover_eggs == 2)

hatching_data2only120$Leftover_eggs <- as.numeric(hatching_data2only120$Leftover_eggs)

Dryepps_summary <- hatching_data2only120[Condition == "Warm" & Dried == "Y", 
                                  .(total_Leftover_eggs = sum(Leftover_eggs, na.rm = TRUE)), 
                                  by = .(plate_clone, Replicate, Condition, Dried, Species)]

Dryepps_summary2 <- hatching_data2only120[Condition == "Warm" & Dried == "Y", 
                                         .N, by = .(plate_clone, Replicate, Condition, Dried, Species)]
Dryepps_summary$putativetotal <- Dryepps_summary2$N*2

onlydryepp<- subset(Dryepps_summary2, Species == "D_obtusa" | Species == "D_pulex")
sum(onlydryepp$N)

Later3 <- Later[, .N, by=list(plate_clone, Replicate, Condition, Dried, Species, time)]

Later4 <- subset(Later3,Condition == "Warm" & Dried == "Y")

Early3 <- Early[, .N, by=list(plate_clone, Replicate, Condition, Dried, Species, time)]
Early4 <- subset(Early3,Condition == "Warm" & Dried == "Y")

Dryepps_summary3 <- full_join(Dryepps_summary, Later4, by=c("plate_clone","Replicate","Condition","Dried", "Species"))
colnames(Dryepps_summary3)[9] ="totalhatchlater"

Dryepps_summary4 <- full_join(Dryepps_summary3, Early4, by=c("plate_clone","Replicate","Condition","Dried", "Species"))
colnames(Dryepps_summary4)[11] ="totalhatchearly"
Dryepps_summary4$totalhatchearly[is.na(Dryepps_summary4$totalhatchearly)] <- 0
Dryepps_summary4$totalhatchlater[is.na(Dryepps_summary4$totalhatchlater)] <- 0

Dryepps_summary4$totalemb <- Dryepps_summary4$totalhatchearly + Dryepps_summary4$totalhatchlater + Dryepps_summary4$total_Leftover_eggs

Dryepps_summary4$fill_rate <- Dryepps_summary4$totalemb/Dryepps_summary4$putativetotal

Dryepps_summary4$se <- sqrt(Dryepps_summary4$fill_rate * (1-Dryepps_summary4$fill_rate)/ (Dryepps_summary4$putativetotal))
Dryepps_summary4$lci <- Dryepps_summary4$fill_rate - (1.96 * Dryepps_summary4$se)
Dryepps_summary4$uci <- Dryepps_summary4$fill_rate + (1.96 * Dryepps_summary4$se)

Dryepps_withspecies <- subset(Dryepps_summary4, Species == "D_obtusa" | Species == "D_pulex")
Dryepps_withspecies_noempty <- subset(Dryepps_withspecies, fill_rate != "0")

Dryepps_withspecies$filled <- round(Dryepps_withspecies$fill_rate * Dryepps_withspecies$totalemb)

fillrate <- Dryepps_summary4[,c(1,2,5,13)]




latertotal<- rbind(Nohatch, Later)

ColdWet2 <- subset(latertotal, Condition == "Cold" & Dried == "N")[, .N, by=list(plate_clone, Replicate, Condition, Dried, Species)]
colnames(ColdWet2)[6] ="total"
ColdWet2$Condition_Dry <- "Cold_Wet"

ColdDry2 <- subset(latertotal, Condition == "Cold" & Dried == "Y")[, .N, by=list(plate_clone, Replicate, Condition, Dried, Species)]
colnames(ColdDry2)[6] ="total"
ColdDry2$Condition_Dry <- "Cold_Dry"

WarmDry2  <- subset(latertotal, Condition == "Warm" & Dried == "Y")[, .N, by=list(plate_clone, Replicate, Condition, Dried, Species)]
colnames(WarmDry2)[6] ="total"
WarmDry2$Condition_Dry <- "Warm_Dry"

Control2  <- subset(latertotal, Condition == "Control")[, .N, by=list(plate_clone, Replicate, Condition, Dried, Species)]
colnames(Control2)[6] ="total"
Control2$Condition_Dry <- "Warm_Wet"

totals <- rbind(ColdWet2,ColdDry2,WarmDry2, Control2)

setDT(all_cond)
setDT(totals)
setkey(all_cond, plate_clone, Replicate, Condition, Dried, Species)
setkey(totals, plate_clone, Replicate, Condition, Dried, Species)
all_cond_totals <- merge(all_cond, totals, all.x=TRUE, all.y=TRUE)
totals$total

sum(subset(totals, Condition_Dry == "Cold_Wet" & Species == "D_pulex")$total)
sum(subset(totals, Condition_Dry == "Cold_Wet" & Species == "D_obtusa")$total)


all_cond_totals$N[is.na(all_cond_totals$N)] <- 0

#all_cond_totals$prop <- all_cond_totals$N/all_cond_totals$total


sum(subset(all_cond_totals, Condition_Dry == "Cold_Wet" & Species == "D_pulex")$N)
sum(subset(all_cond_totals, Condition_Dry == "Cold_Wet" & Species == "D_obtusa")$N)

sum(subset(all_cond_totals, Species == "D_pulex")$N)
sum(subset(all_cond_totals, Species == "D_obtusa")$N)

sum(subset(totals, Species == "D_pulex")$total)
sum(subset(totals, Species == "D_obtusa")$total)

#all_cond_totals$se <- sqrt(all_cond_totals$prop * (1-all_cond_totals$prop)/ (all_cond_totals$total))
#all_cond_totals$lci <- all_cond_totals$prop - (1.96 * all_cond_totals$se)
#all_cond_totals$uci <- all_cond_totals$prop + (1.96 * all_cond_totals$se)

#Male data
maledata <- read.csv("C:/Users/rjpor/Downloads/Meso Wk 5 Count - Sheet1.csv")

#Get rid of the old mesocosms with "other" algae
maledata <- subset(maledata, Clone!="222 Reed")
maledata <- subset(maledata, Clone!="222 Old")

maledata$pondclone <- paste(maledata$Pond, maledata$Clone)

maledata$percentmale <- maledata$Males.L/(maledata$Males.L+maledata$Females.L)

write.csv(maledata, "C:/Users/rjpor/Downloads/maledata.csv")

#Fill rate

fillrate <- Dryepps_summary4[,c(1,2,5,13)]

all_cond_fill <- left_join(all_cond_totals, fillrate, by=c("plate_clone","Replicate", "Species"))
subset(all_cond_fill, Species == "D_pulex")

all_cond_fill <- all_cond_fill[,c(1,2,5,6,7,8,9)]
#Total here would be the total of putative offspring if each ephippia were filled
all_cond_fill$realtotal <- all_cond_fill$total*all_cond_fill$fill_rate

all_cond_fill$prop_fixed <- all_cond_fill$N/(all_cond_fill$total*all_cond_fill$fill_rate)
all_cond_fill$se <- sqrt(all_cond_fill$prop_fixed * (1-all_cond_fill$prop_fixed)/ (all_cond_fill$total*all_cond_fill$fill_rate))
all_cond_fill$lci <- all_cond_fill$prop - (1.96 * all_cond_fill$se)
all_cond_fill$uci <- all_cond_fill$prop + (1.96 * all_cond_fill$se)


sum(subset(all_cond_fill, Species == "D_pulex")$realtotal)
sum(subset(all_cond_fill, Species == "D_obtusa")$realtotal)
subset(all_cond_fill, Species == "D_pulex")

sum(subset(all_cond_fill, Species == "D_pulex" & Condition_Dry == "Cold_Wet")$realtotal)
sum(subset(all_cond_fill, Species == "D_obtusa"& Condition_Dry == "Cold_Wet")$realtotal)


# D. obtusa
prop.test(x = 2232, n = 16745, correct = FALSE)

# D. pulex
prop.test(x = 152, n = 1296, correct = FALSE)

#Proportion early hatched

prop.test(x = 10, n = 1296, correct = FALSE)
prop.test(x = 51, n = 16745, correct = FALSE)

prop.test(x = 145, n = 539, correct = FALSE)
prop.test(x = 2161, n = 6849, correct = FALSE)


hatchshortsmol2time <- subset(hatchingscoldwet, time =="later"| time== "early")[, .N, by=list(plate_clone, Replicate,Species, time)]

lateronly <- subset(hatchshortsmol2time, time == "later")
earlyonly <- subset(hatchshortsmol2time, time == "early")

onlyhatch<- full_join(lateronly,earlyonly, by=c("plate_clone","Replicate","Species"))
onlyhatch$N.y[is.na(onlyhatch$N.y)] <- 0
onlyhatch$total <- onlyhatch$N.x + onlyhatch$N.y

onlyhatch$prop <- onlyhatch$N.y/onlyhatch$total


onlyhatch$se <- sqrt(onlyhatch$prop * (1-onlyhatch$prop)/ (onlyhatch$total))
onlyhatch$lci <- onlyhatch$prop - (1.96 * onlyhatch$se)
onlyhatch$uci <- onlyhatch$prop + (1.96 * onlyhatch$se)

all_cond_filltotals <- all_cond_fill[,c(1,2,3,6,8)]

dt_summed <- all_cond_filltotals[, .(realtotal = sum(realtotal)), by = .(plate_clone, Replicate, Species)]

onlyearlyhatched <- onlyhatch[,c(1,2,3,4,6,7)]

earlytotals <- left_join(onlyearlyhatched, dt_summed, by=c("plate_clone","Replicate", "Species"))

sum(subset(earlytotals, Species == "D_pulex")$realtotal)
sum(subset(earlytotals, Species == "D_obtusa")$realtotal)









setDT(hatchshort2early)
setDT(hatchshortsmol2)
setkey(hatchshort2early, plate_clone, Replicate)
setkey(hatchshortsmol2, plate_clone, Replicate)
hatchshortmergeshort <- merge(hatchshort2early, hatchshortsmol2, all.x=TRUE, all.y=TRUE)

hatchshortmergeshort$N[is.na(hatchshortmergeshort$N)] <- 0


hatchshortmergeshort$prop <- hatchshortmergeshort$N/hatchshortmergeshort$total


colnames(hatchshortmergeshort)[3] ="Hatching_time"



hatchshortmergeshort$se <- sqrt(hatchshortmergeshort$prop * (1-hatchshortmergeshort$prop)/ (hatchshortmergeshort$total))
hatchshortmergeshort$lci <- hatchshortmergeshort$prop - (1.96 * hatchshortmergeshort$se)
hatchshortmergeshort$uci <- hatchshortmergeshort$prop + (1.96 * hatchshortmergeshort$se)

#Diversified Hatching

Diverse <- fread("C:/Users/rjpor/Downloads/Diversified_Hatching2.csv")


Diverse <- subset(Diverse, Hatch3 =="")
Diverse <- Diverse[,c(1,2,3,4,5,6,7,8,9,10,11,12,17,18,19,20,23)]

unique(Diverse$Leftover_eggs)
dissected <- subset(Diverse, Leftover_eggs != "")
unique(dissected$Leftover_eggs)

sum(dissected$Leftover_eggs)/(2*192)

Diversea <- Diverse[,c(1,2,3,4,5,6,7,8,9,10,13,14,15,16)]
Diversea$sib <- "a"
Diverseb <- Diverse[,c(1,2,3,4,5,6,7,8,11,12,13,14,15,16)]
Diverseb$sib <- "b"
colnames(Diverseb)[9] ="Hatch1"
colnames(Diverseb)[10] ="Hatch1HD"

Diverselong <- rbind(Diversea,Diverseb)

Diverselong$Hatch1 <- str_replace(Diverselong$Hatch1, "X", "6/7")

Diverselong$Hatch1 <- as.Date(Diverselong$Hatch1, "%m/%d")
Diverselong$`Date Collected` <- as.Date(Diverselong$`Date Collected`, "%m/%d/24")

Diverselong$h1 <- Diverselong$Hatch1- Diverselong$`Date Collected` 
Diverselong$h1[is.na(Diverselong$h1)] <- -100



LaterDiverse<- subset(Diverselong, h1 > 45)
LaterDiverse$time <- "LaterDiverse"

EarlyDiverse <- subset(Diverselong, h1 < 45 & h1 > 0)
EarlyDiverse$time <- "EarlyDiverse"

nohatch <-  subset(Diverselong, h1 < 0)
nohatch$time <- "nohatch"

ColdWetLight <- subset(LaterDiverse, Temp == "Cold" & Dried == "Wet" & Light =="Light")[, .N, by=list(Temp, Dried, Light)]
ColdDryLight <- subset(LaterDiverse, Temp == "Cold" & Dried == "Dry" & Light =="Light")[, .N, by=list(Temp, Dried, Light)]
WarmDryLight  <- subset(LaterDiverse, Temp == "Warm" & Dried == "Dry" & Light =="Light")[, .N, by=list(Temp, Dried, Light)]
WarmWetLight <- subset(LaterDiverse, Temp == "Warm" & Dried == "Wet" & Light =="Light")[, .N, by=list(Temp, Dried, Light)]
ColdWetDark <- subset(LaterDiverse, Temp == "Cold" & Dried == "Wet" & Light =="Dark")[, .N, by=list(Temp, Dried, Light)]
ColdDryDark <- subset(LaterDiverse, Temp == "Cold" & Dried == "Dry" & Light =="Dark")[, .N, by=list(Temp, Dried, Light)]
WarmDryDark  <- subset(LaterDiverse, Temp == "Warm" & Dried == "Dry" & Light =="Dark")[, .N, by=list(Temp, Dried, Light)]
WarmWetDark <- subset(LaterDiverse, Temp == "Warm" & Dried == "Wet" & Light =="Dark")[, .N, by=list(Temp, Dried, Light)]

all_condRep <- rbind(ColdWetLight,ColdDryLight,WarmDryLight,WarmWetLight,ColdWetDark,ColdDryDark,WarmDryDark,WarmWetDark)

ColdWetLight2 <- subset(Diverselong, Temp == "Cold" & Dried == "Wet" & Light =="Light")[, .N, by=list(Temp, Dried, Light)]
ColdDryLight2 <- subset(Diverselong, Temp == "Cold" & Dried == "Dry" & Light =="Light")[, .N, by=list(Temp, Dried, Light)]
WarmDryLight2  <- subset(Diverselong, Temp == "Warm" & Dried == "Dry" & Light =="Light")[, .N, by=list(Temp, Dried, Light)]
WarmWetLight2 <- subset(Diverselong, Temp == "Warm" & Dried == "Wet" & Light =="Light")[, .N, by=list(Temp, Dried, Light)]
ColdWetDark2 <- subset(Diverselong, Temp == "Cold" & Dried == "Wet" & Light =="Dark")[, .N, by=list(Temp, Dried, Light)]
ColdDryDark2 <- subset(Diverselong, Temp == "Cold" & Dried == "Dry" & Light =="Dark")[, .N, by=list(Temp, Dried, Light)]
WarmDryDark2  <- subset(Diverselong, Temp == "Warm" & Dried == "Dry" & Light =="Dark")[, .N, by=list(Temp, Dried, Light)]
WarmWetDark2 <- subset(Diverselong, Temp == "Warm" & Dried == "Wet" & Light =="Dark")[, .N, by=list(Temp, Dried, Light)]

all_condReptotals <- rbind(ColdWetLight2,ColdDryLight2,WarmDryLight2,WarmWetLight2,ColdWetDark2,ColdDryDark2,WarmDryDark2,WarmWetDark2)

colnames(all_condReptotals)[4] ="total"

setDT(all_condRep)
setDT(all_condReptotals)
setkey(all_condRep, Temp, Dried, Light)
setkey(all_condReptotals, Temp, Dried, Light)
all_conds <- merge(all_condRep, all_condReptotals, all.x=TRUE, all.y=TRUE)

all_conds$N[is.na(all_conds$N)] <- 0
all_conds$prop <- all_conds$N/(all_conds$total*0.6588542)
all_conds$fixed <- all_conds$total*0.6588542
sum(all_conds$total*0.6588542)

all_conds$se <- sqrt(all_conds$prop * (1-all_conds$prop)/ (all_conds$total*0.6588542))
all_conds$lci <- all_conds$prop - (1.96 * all_conds$se)
all_conds$uci <- all_conds$prop + (1.96 * all_conds$se)

sum(subset(all_conds, N==12|N==14|N==22|N==2)$N)/sum(subset(all_conds, N==12|N==14|N==22|N==2)$fixed)

prop.test(x = sum(subset(all_conds, N==12|N==14|N==22|N==2)$N), n = sum(subset(all_conds, N==12|N==14|N==22|N==2)$fixed), correct = FALSE)

#map of sites

library(maps)

sites <- read.csv("C:\\Users\\rjpor\\Downloads\\Sites.csv")

sites2 <- subset(sites, Pond== "Elvis" | Pond == "Stav1" | Pond == "Birdhut")

#Data
worldmap = map_data('world')
#Option 2
b<- ggplot() + geom_polygon(data = worldmap, 
                            aes(x = long, 
                                y = lat, 
                                group = group,
                            ),
                            fill = 'white', 
                            color = 'black'
                            ) + 
  #  coord_fixed(xlim = c(-2.15,-2.05), 
  coord_fixed(xlim = c(-4,1), 
              #  ylim = c(50.6, 50.7))+
              ylim = c(50, 56))+
  theme_bw()+ geom_jitter(data=subset(sites, Site != ""), aes(x=Long, y=Lat, color=Site), pch=1,
                          size=2, width = 0.05, height = 0.05)+
  # geom_jitter(width = 1000, height = 1000)+ 
  xlab("Longitude")+ylab("Latitude")
b

setDT(sites)
numofspec <- sites[, .N, by=list(Site, Putative_species)]


ggsave("C:/Users/rjpor/Downloads/MapofEngland.png", plot = b, width = 6, height = 4, dpi = 300)

#pH graph
locations <- fread("C:/Users/rjpor/Downloads/PondLocationsMar2023.csv")

ggplot(locations,aes(x=ph1, y= ph2)) +
  geom_point()+
  theme_bw()+ ylab("pH2")+ xlab("pH1")+
  geom_abline(slope=1)+ xlim(6,8.5)+ylim(6,8.5)

aaa <- glm(ph2~ph1, locations, family=poisson())
anova(aaa, test = "Chisq")

lm(locations$ph2~locations$ph1)

PH <- ggplot(locations,aes(x=ph1, col=species)) +
  geom_histogram()+
  theme_bw()+ ylab("Number of Ponds")+ xlab("pH")

Depth <- ggplot(locations,aes(x=depth_ft, group=Pond)) +
  geom_histogram(bins=10)+
  theme_bw()+ xlab("Depth in feet")+ ylab("")


Genhard <- ggplot(locations,aes(x=GeneralHardness, group=Pond)) +
  geom_histogram(bins=15)+
  theme_bw()+ ylab("")+ xlab("General Hardness")

ggplot(subset(locations, species != "NA"),aes(x=depth_ft,y = ph1, group=Pond, color = species)) +
  geom_point(size=2)+ 
  theme_bw()+ ylab("pH")+ xlab("Depth (ft)")

ggplot(subset(locations, species != "NA"),aes(x=depth_ft,y = GeneralHardness, group=Pond, color = species)) +
  geom_point(size=2)+ 
  theme_bw()+ ylab("GeneralHardness")+ xlab("Depth (ft)")

ggplot(subset(locations, species != "NA"),aes(x=depth_ft,y = watertemp, group=Pond, color = species)) +
  geom_point(size=2)+ 
  theme_bw()+ ylab("TempC")+ xlab("Depth (ft)")

PH+Depth+Genhard

pdf("C:\\Users\\rjpor\\Desktop\\conditions.pdf", width=8, height=4)
PH+Depth+Genhard
dev.off()

locationsspecfixed <- locations
locationsspecfixed <- subset(locationsspecfixed, species != "NA")
locationsspecfixed$species <- gsub("D_pulex", "D. pulex", locationsspecfixed$species)
locationsspecfixed$species <- gsub("D_obtusa", "D. obtusa", locationsspecfixed$species)




phxDepth <- ggplot(subset(locationsspecfixed, species != "NA"),aes(x=depth_ft,y = ph1, group=Pond, color = species)) +
  geom_point(size=2)+ xlim(0,6)+
  theme_bw()+ ylab("pH")+ xlab("Depth (ft)")
phxDepth
ggsave("C:/Users/rjpor/Downloads/phxDepth.pdf", plot = phxDepth, width = 4, height = 3.5, dpi = 600)



phbox <- ggplot(data=locationsspecfixed, aes(x=species, y=ph1)) +
  geom_boxplot()+
  ylab("pH") +
  xlab("species") + theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_signif(comparisons = list(c("D. pulex", "D. obtusa")),
              map_signif_level = TRUE)

depthbox <- ggplot(data=locationsspecfixed, aes(x=species, y=depth_ft)) +
  geom_boxplot()+
  ylab("depth_ft") +
  xlab("species") + theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_signif(comparisons = list(c("D. pulex", "D. obtusa")),
              map_signif_level = TRUE)






pkgs <- unique(c(
  "stringr","tidyr","tibble","dplyr","ggplot2","data.table",
  "ggtree","ape","foreach","lattice","tidyverse","gridExtra",
  "reshape2","patchwork","knitr","lme4","cowplot","ggsignif"
))

get_apa_citation <- function(pkg) {
  cit <- citation(pkg)
  txt <- capture.output(print(cit, style = "text"))
  cat("\n\n---", pkg, "---\n")
  cat(paste(txt, collapse = "\n"))
}

for (p in pkgs) {
  get_apa_citation(p)
}

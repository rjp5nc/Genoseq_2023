
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
library(ggpubr)
library(ggpmisc)







Timeseries <- fread("C:/Users/rjpor/Downloads/Aim1 Mesocosm Hatching - Time Series Hatching (current to 9_3_24).csv")


timeseriesa <- Timeseries[,c(1,2,3,4,5,6,7,8,9,12)]
timeseriesa$sib <- "a"
timeseriesb <- Timeseries[,c(1,2,3,4,5,6,7,10,11,12)]
timeseriesb$sib <- "b"

colnames(timeseriesb)[8] ="Hatch1"
colnames(timeseriesb)[9] ="Hatch1HD"
#rename 9 and 10

timeserieslong <- rbind(timeseriesa,timeseriesb)

timeserieslong$Hatch1 <- as.Date(paste(timeserieslong$Hatch1, "2024"), format = "%d-%b %Y")
timeserieslong$DateCollected <- as.Date(timeserieslong$`Date Collected`, format = "%m/%d/%Y")

# Calculate the number of days between collected and hatch date
timeserieslong$h1 <- as.numeric(timeserieslong$Hatch1 - timeserieslong$DateCollected)


timeserieslongnona <- na.omit(timeserieslong)

timeseriesshortRep <- timeserieslongnona[, .N, by=list(Plate,Replicate,Time)]
timeseriesshort <- timeserieslongnona[, .N, by=list(Plate,Time)]

timeseriesshort$N <- as.numeric(timeseriesshort$N)

ggplot(data=timeseriesshort, aes(x=Time, y=N)) +
  geom_point(size=2)+ ggtitle("Number of Hatchlings Over Time")+
  ylab("Number of Hatchlings") + # facet_wrap("Replicate") +
  xlab("Time in Cold Treatment (Weeks)") + theme_bw() + ylim(0,45)


timeseriesshort$prop <- timeseriesshort$N/(2*72*0.6588542)

ggplot(data=timeseriesshort, aes(x=Time, y=prop)) +
  geom_point(size=2)+ ggtitle("Number of Hatchlings Over Time")+
  ylab("Number of Hatchlings") + # facet_wrap("Replicate") +
  xlab("Time in Cold Treatment (Weeks)") + theme_bw() + ylim(0,1)+ 
  stat_smooth(method = "lm", formula = y ~ x+ I(x^2), se= F)+
  stat_fit_glance(method = 'lm', label.x="middle", label.y="bottom",
                  method.args = list(formula =  y ~ x+ I(x^2)),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                  size = 3)

pdf("C:\\Users\\rjpor\\Downloads\\time_series.pdf", width=4, height=3)

ggplot(data=timeseriesshort, aes(x=Time, y=N)) +
  geom_point(size=2)+ ggtitle("Number of Hatchlings Over Time")+
  ylab("Number of Hatchlings") + 
  xlab("Time in Cold Treatment (Weeks)") + theme_bw() + ylim(0,45)+ 
  stat_smooth(method = "lm", formula = y ~ x+ I(x^2), se= F)+
  stat_fit_glance(method = 'lm', label.x="middle", label.y="bottom",
                method.args = list(formula =  y ~ x+ I(x^2)),
                geom = 'text',
                aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                 size = 3)

dev.off()



ggplot(data=timeseriesshortRep, aes(x=Time, y=N, color= Replicate)) +
  #geom_point(size=2)+
 # geom_line()+
  geom_point()+
  ylab("Number of Hatchlings") + # facet_wrap("Replicate") +
  xlab("Time in Cold (weeks)") + theme_bw()















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



Later<- subset(Diverselong, h1 > 45)
Later$time <- "later"

early <- subset(Diverselong, h1 < 45 & h1 > 0)
early$time <- "early"

nohatch <-  subset(Diverselong, h1 < 0)
nohatch$time <- "nohatch"


ColdWetLight <- subset(Later, Temp == "Cold" & Dried == "Wet" & Light =="Light")[, .N, by=list(Replicate, Temp, Dried, Light)]
ColdDryLight <- subset(Later, Temp == "Cold" & Dried == "Dry" & Light =="Light")[, .N, by=list(Replicate, Temp, Dried, Light)]
WarmDryLight  <- subset(Later, Temp == "Warm" & Dried == "Dry" & Light =="Light")[, .N, by=list(Replicate, Temp, Dried, Light)]
WarmWetLight <- subset(Later, Temp == "Warm" & Dried == "Wet" & Light =="Light")[, .N, by=list(Replicate, Temp, Dried, Light)]
ColdWetDark <- subset(Later, Temp == "Cold" & Dried == "Wet" & Light =="Dark")[, .N, by=list(Replicate, Temp, Dried, Light)]
ColdDryDark <- subset(Later, Temp == "Cold" & Dried == "Dry" & Light =="Dark")[, .N, by=list(Replicate, Temp, Dried, Light)]
WarmDryDark  <- subset(Later, Temp == "Warm" & Dried == "Dry" & Light =="Dark")[, .N, by=list(Replicate, Temp, Dried, Light)]
WarmWetDark <- subset(Later, Temp == "Warm" & Dried == "Wet" & Light =="Dark")[, .N, by=list(Replicate, Temp, Dried, Light)]


all_condRep <- rbind(ColdWetLight,ColdDryLight,WarmDryLight,WarmWetLight,ColdWetDark,ColdDryDark,WarmDryDark,WarmWetDark)


ColdWetLight2 <- subset(Diverselong, Temp == "Cold" & Dried == "Wet" & Light =="Light")[, .N, by=list(Replicate, Temp, Dried, Light)]
ColdDryLight2 <- subset(Diverselong, Temp == "Cold" & Dried == "Dry" & Light =="Light")[, .N, by=list(Replicate, Temp, Dried, Light)]
WarmDryLight2  <- subset(Diverselong, Temp == "Warm" & Dried == "Dry" & Light =="Light")[, .N, by=list(Replicate, Temp, Dried, Light)]
WarmWetLight2 <- subset(Diverselong, Temp == "Warm" & Dried == "Wet" & Light =="Light")[, .N, by=list(Replicate, Temp, Dried, Light)]
ColdWetDark2 <- subset(Diverselong, Temp == "Cold" & Dried == "Wet" & Light =="Dark")[, .N, by=list(Replicate, Temp, Dried, Light)]
ColdDryDark2 <- subset(Diverselong, Temp == "Cold" & Dried == "Dry" & Light =="Dark")[, .N, by=list(Replicate, Temp, Dried, Light)]
WarmDryDark2  <- subset(Diverselong, Temp == "Warm" & Dried == "Dry" & Light =="Dark")[, .N, by=list(Replicate, Temp, Dried, Light)]
WarmWetDark2 <- subset(Diverselong, Temp == "Warm" & Dried == "Wet" & Light =="Dark")[, .N, by=list(Replicate, Temp, Dried, Light)]


all_condReptotals <- rbind(ColdWetLight2,ColdDryLight2,WarmDryLight2,WarmWetLight2,ColdWetDark2,ColdDryDark2,WarmDryDark2,WarmWetDark2)

colnames(all_condReptotals)[5] ="total"


setDT(all_condRep)
setDT(all_condReptotals)
setkey(all_condRep, Replicate, Temp, Dried, Light)
setkey(all_condReptotals, Replicate, Temp, Dried, Light)
all_conds <- merge(all_condRep, all_condReptotals, all.x=TRUE, all.y=TRUE)

all_conds$N[is.na(all_conds$N)] <- 0


all_conds$prop <- all_conds$N/(all_conds$total*0.6588542)



all_conds$se <- sqrt(all_conds$prop * (1-all_conds$prop)/ (all_conds$total*0.6588542))
all_conds$lci <- all_conds$prop - (1.96 * all_conds$se)
all_conds$uci <- all_conds$prop + (1.96 * all_conds$se)


conditionglm <- glm(prop~Light*Temp*Dried, all_conds, family=binomial(), weights = total)
anova(conditionglm, test = "Chisq")

allcond_plot<- ggplot(data=all_conds, aes(x=Dried, y=prop, col=Temp, shape= Light)) +
  geom_point(size=2, position = position_dodge(width = 0.5))+
  ylab("Proportion hatched") + facet_wrap("Replicate") +
  geom_errorbar(aes(ymin=uci, ymax=lci), width=.2, position = position_dodge(width = 0.5)) +
  xlab("Condition") + theme_bw() + ylim(-.1,1)
allcond_plot


allcond_plot2 <- ggplot(data=all_cond_totals, aes(x=Condition_Dry, y=total, col=Lineage)) +
  geom_point(size=2, position = position_dodge(width = 0.5))+
  ylab("total ephippia in condition") + facet_wrap("Replicate") +
  xlab("Condition") + theme_bw()

allcond_plot3 <- ggplot(data=all_cond_totals, aes(x=Condition_Dry, y=N, col=Lineage)) +
  geom_point(size=2, position = position_dodge(width = 0.5))+
  ylab("total hatchlings from condition") + facet_wrap("Replicate") +
  xlab("Condition") + theme_bw()



pdf("C:\\Users\\rjpor\\Desktop\\EarlyHatchingAim1\\Allcond.pdf", width=8, height=6)

allcond_plot + ggtitle("Proportion of Late Hatchers")

dev.off()






























ColdWetLight <- subset(Later, Temp == "Cold" & Dried == "Wet" & Light =="Light")[, .N, by=list(Temp, Dried, Light)]
ColdDryLight <- subset(Later, Temp == "Cold" & Dried == "Dry" & Light =="Light")[, .N, by=list(Temp, Dried, Light)]
WarmDryLight  <- subset(Later, Temp == "Warm" & Dried == "Dry" & Light =="Light")[, .N, by=list(Temp, Dried, Light)]
WarmWetLight <- subset(Later, Temp == "Warm" & Dried == "Wet" & Light =="Light")[, .N, by=list(Temp, Dried, Light)]
ColdWetDark <- subset(Later, Temp == "Cold" & Dried == "Wet" & Light =="Dark")[, .N, by=list(Temp, Dried, Light)]
ColdDryDark <- subset(Later, Temp == "Cold" & Dried == "Dry" & Light =="Dark")[, .N, by=list(Temp, Dried, Light)]
WarmDryDark  <- subset(Later, Temp == "Warm" & Dried == "Dry" & Light =="Dark")[, .N, by=list(Temp, Dried, Light)]
WarmWetDark <- subset(Later, Temp == "Warm" & Dried == "Wet" & Light =="Dark")[, .N, by=list(Temp, Dried, Light)]


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



all_conds$se <- sqrt(all_conds$prop * (1-all_conds$prop)/ (all_conds$total*0.6588542))
all_conds$lci <- all_conds$prop - (1.96 * all_conds$se)
all_conds$uci <- all_conds$prop + (1.96 * all_conds$se)


allcond_plot<- ggplot(data=all_conds, aes(x=Light, y=prop, col=Temp)) +
  geom_point(size=2, position = position_dodge(width = 0.5))+
  ylab("Proportion hatched") + facet_wrap("Dried") +
  geom_errorbar(aes(ymin=uci, ymax=lci), width=.2, position = position_dodge(width = 0.5)) +
  xlab("Condition") + theme_bw() + ylim(-.1,1)
allcond_plot

conditionglm <- glm(prop~Light*Temp*Dried, all_conds, family=binomial(), weights = total)
anova(conditionglm, test = "Chisq")


pdf("C:\\Users\\rjpor\\Downloads\\AllcondDiversehatching.pdf", width=6, height=4)

allcond_plot + ggtitle("Diversified Hatching")

dev.off()



ggsave("C:/Users/rjpor/Downloads/AllcondDiversehatching.png", plot = allcond_plot, width = 6, height = 4, dpi = 300)

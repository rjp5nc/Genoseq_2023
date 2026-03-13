#Aim 2 figures



popbox <- ggplot(data=maledata, aes(x=Species, y=Indiv.Meso)) +
  geom_boxplot()+
  ylab("Number in Meso") + theme_bw()+ ylim(0,4000)+
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
geom_signif(comparisons = list(c("D. pulex", "D. obtusa")),
              test = "t.test",
              map_signif_level = TRUE)


save(maledata, file = "/Users/rjpor/Downloads/maledata.RData")



eppbox <- ggplot(data=NumofEpp2, aes(x=Species, y=N)) +
  geom_boxplot()+
  ylab("Number of Ephippia")+ 
  xlab("Species") + theme_bw()+ ylim(0,1100)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),   
        axis.text.y  = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  geom_signif(comparisons = list(c("D_pulex", "D_obtusa")),
              map_signif_level = TRUE)


mergedNumofEpp2 <- merge(NumofEpp2, Mesolist, by = c("Pond", "Clone", "Replicate", "Species", "plate_clone"), all = TRUE)


numofeppplot <- ggplot(data=mergedNumofEpp2, aes(x=plate_clone, y=N, group=plate_clone, col=Replicate)) +
  geom_point(size=2, position = position_dodge(width = 0.5))+
  ylab("")+ facet_wrap("Species", scales = "free_x")+
  xlab("Clone") + theme_bw()+ ylim(0,1100)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),   
        axis.text.y  = element_text(size = 14),
        axis.title.y = element_text(size = 14))



maledata2 <- maledata

maledata2[29,1] <- "PottCarr1"
maledata2[29, 2] <- "2"
maledata2[30, 1] <- "PottCarr1"
maledata2[30, 2] <- "2"
maledata2[31, 1] <- "PottCarr3"
maledata2[31, 2] <- "50.2"
maledata2[32, 1] <- "PottCarr3"
maledata2[32, 2] <- "50.2"

maledata2$Pond <- gsub(" ","",maledata2$Pond)

mergedmaledata <- merge(maledata2, Mesolist, by.x = c("Pond", "Clone", "Rep"), by.y = c("Pond", "Clone", "Replicate"), all = TRUE)

popplot <- ggplot(data=mergedmaledata, aes(x=pondclone, y=Indiv.Meso, col=Rep)) +
  geom_point(size=2)+
  ylab("")+ facet_wrap("Species.x", scales = "free_x")+  ylab("Number in Meso")+ 
  xlab("Clone") + theme_bw()+ ylim(0,4000)+ ylab("")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),   
        axis.text.y  = element_text(size = 14),
        axis.title.y = element_text(size = 14))


mergedmaledata
mergedmaledata$se <- sqrt(mergedmaledata$percentmale * (1-mergedmaledata$percentmale)/ (mergedmaledata$Indiv.L))
mergedmaledata$lci <- mergedmaledata$percentmale - (1.96 * mergedmaledata$se)
mergedmaledata$uci <- mergedmaledata$percentmale + (1.96 * mergedmaledata$se)


maledataplot <- ggplot(mergedmaledata, aes(x = pondclone, y= percentmale, col=Rep)) +
  geom_point(size=2, position = position_dodge(width = 0.5))+
  ylab("")+ facet_wrap("Species.x", scales = "free_x")+
  xlab("Clone") + theme_bw()+ ylim(0,1)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),   
        axis.text.y  = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  geom_errorbar(aes(ymin=uci, ymax=lci), width=.2,position = position_dodge(width = 0.5))
  

malebox <- ggplot(data=maledata, aes(x=Species, y=percentmale)) +
  geom_boxplot()+
  ylab("PercentMale")+ ylim(0,1)+
  xlab("Species") + theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),   
        axis.text.y  = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  geom_signif(comparisons = list(c("D. pulex", "D. obtusa")),
              map_signif_level = TRUE)



Dryepps_summary4_separated <- Dryepps_summary4 %>%
  separate(plate_clone, into = c("Pond", "Clone", "oth"), sep = "_")

mergedDryepps_summary4 <- merge(Dryepps_summary4_separated, Mesolist, by = c("Pond", "Clone", "Replicate", "Species"), all = TRUE)

fillpoint <- ggplot(subset(mergedDryepps_summary4, Species == "D_obtusa" | Species == "D_pulex"), aes(x=plate_clone, y=fill_rate, col = Replicate)) +
  geom_point(size=2, position = position_dodge(width = 0.5))+
  ylab("") + ylim(0,1.2)+
  geom_errorbar(aes(ymin=uci, ymax=lci), width=.2, position = position_dodge(width = 0.5))+
  xlab("Clone") + theme_bw()+ facet_wrap("Species", scales = "free_x")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),   
        axis.text.y  = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

fillbox <- ggplot(subset(Dryepps_summary4, Species == "D_obtusa" | Species == "D_pulex"), aes(x=Species, y=fill_rate)) +
  geom_boxplot()+
  ylab("Fill Rate") +
  xlab("Species") + theme_bw()+ ylim(0,1.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y  = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  geom_signif(comparisons = list(c("D_pulex", "D_obtusa")),
              test = "t.test",
              map_signif_level = TRUE)





onlyhatch_separated <- earlytotals %>%
  separate(plate_clone, into = c("Pond", "Clone", "oth"), sep = "_")

fills <- mergedDryepps_summary4[,c(1,2,3,4,9,15)]

totals2222 <- left_join(Mesolist, fills, by = c("Pond", "Clone", "Replicate", "Species"))

mergedonlyhatch <- merge(onlyhatch_separated, totals2222, by = c("Pond", "Clone", "Replicate", "Species"), all = TRUE)

mergedonlyhatch$realtotal <- mergedonlyhatch$putativetotal*mergedonlyhatch$fill_rate

sum(subset(mergedonlyhatch, Species == "D_pulex")$realtotal, na.rm=TRUE)

mergedonlyhatch$N.y[is.na(mergedonlyhatch$N)] <- 0

mergedonlyhatch$prop <- mergedonlyhatch$N.y/mergedonlyhatch$realtotal
mergedonlyhatch$se <- sqrt(mergedonlyhatch$prop * (1-mergedonlyhatch$prop)/ (mergedonlyhatch$realtotal))
mergedonlyhatch$lci <- mergedonlyhatch$prop - (1.96 * mergedonlyhatch$se)
mergedonlyhatch$uci <- mergedonlyhatch$prop + (1.96 * mergedonlyhatch$se)


earlyprop2 <- ggplot(data=mergedonlyhatch, aes(x=plate_clone, y=prop, col=Replicate)) +
  geom_point(size=2, position = position_dodge(width = 0.5))+
  labs(y = "", x = "Clone") +
  geom_errorbar(aes(ymin=uci, ymax=lci), width=.2, position = position_dodge(width = 0.5)) + theme_bw()+  facet_wrap("Species", scales = "free_x")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),   
        axis.text.y  = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

onlyearlyhatchpulex <- subset(mergedonlyhatch, Species == "D_pulex")
onlyearlyhatchobtusa <- subset(mergedonlyhatch, Species == "D_obtusa")
sum(onlyearlyhatchpulex$total, na.rm = TRUE)
sum(onlyearlyhatchpulex$N.y, na.rm = TRUE)
sum(onlyearlyhatchobtusa$total, na.rm = TRUE)
sum(onlyearlyhatchobtusa$N.y, na.rm = TRUE)

sum(mergedonlyhatch$N.y, na.rm = TRUE)

sum(onlyearlyhatchpulex$N.x, na.rm = TRUE)
sum(onlyearlyhatchobtusa$N.x, na.rm = TRUE)


earlybox2 <- ggplot(data=mergedonlyhatch, aes(x=Species, y=prop)) +
  geom_boxplot()+
  ylab("Proportion of Early Hatchers\nout of total embryos") +
  xlab("Species") + theme_bw()+ ylim(0,0.8)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),   
        axis.text.y  = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  geom_signif(comparisons = list(c("D_pulex", "D_obtusa")),
              map_signif_level = TRUE)


allcond_plot2<- ggplot(data=subset(all_cond_fill, Species == "D_obtusa"|Species=="D_pulex"), aes(x=Condition_Dry, y=prop_fixed, group=plate_clone)) +
  geom_point(size=2, position = position_dodge(width = 0.5))+
  ylab("") + facet_wrap("Species") +
  geom_errorbar(aes(ymin=uci, ymax=lci), width=.2, position = position_dodge(width = 0.5)) +
  xlab("Condition") + theme_bw() + ylim(-.0,1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),   
        axis.text.y  = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))


testplot2 <-ggplot(data=subset(all_cond_fill, Species == "D_obtusa"|Species=="D_pulex"), aes(x=Condition_Dry, y=prop_fixed, col=Species, group=Replicate)) +
  geom_point(size=2)+
  geom_line()+
  ylab("Proportion hatched") + facet_wrap("plate_clone") +
  geom_errorbar(aes(ymin=uci, ymax=lci), width=.2) +
  xlab("Condition") + theme_bw() + ylim(-.0,1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

test2 <- subset(all_cond_fill, Condition_Dry == "Cold_Dry" | Condition_Dry == "Cold_Wet")
testplot2 <-ggplot(data=subset(test2, Species == "D_obtusa"|Species=="D_pulex" ), aes(x=Condition_Dry, y=prop_fixed, col=Species, group=Replicate)) +
  geom_point(size=2)+
  geom_line()+
  ylab("Proportion hatched") + facet_wrap("plate_clone") +
  geom_errorbar(aes(ymin=uci, ymax=lci), width=.2) +
  xlab("Condition") + theme_bw() + ylim(-.0,1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ReactionNormsallplot <-ggplot(data=subset(all_cond_fill, Species == "D_obtusa"|Species=="D_pulex" ), aes(x=Condition_Dry, y=prop_fixed, col=Species, group=Replicate)) +
  geom_point(size=2)+
  geom_line()+
  ylab("Proportion hatched during late time out of total embryos") + facet_wrap("plate_clone") +
  geom_errorbar(aes(ymin=uci, ymax=lci), width=.2) +
  xlab("Condition") + theme_bw() + ylim(-.0,1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("C:/Users/rjpor/Downloads/ReactionNormsallplot.png", plot = ReactionNormsallplot, width = 6, height = 8, dpi = 300)



all_cond_fill77 <- all_cond_fill
all_cond_fill77$clonerep <- paste(all_cond_fill77$plate_clone, all_cond_fill77$Rep)
ReactionNormsallplotbyspec <-ggplot(data=subset(all_cond_fill77, Species == "D_obtusa"|Species=="D_pulex" ), aes(x=Condition_Dry, y=prop_fixed, col=plate_clone, group=clonerep)) +
  geom_point(size=2)+
  geom_line()+
  ylab("Proportion hatched") + facet_wrap("Species") +
  geom_errorbar(aes(ymin=uci, ymax=lci), width=.2) +
  xlab("Condition") + theme_bw() + ylim(-.0,1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y  = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))


allcond_box2<- ggplot(data=subset(all_cond_fill, Species == "D_obtusa"|Species=="D_pulex"), aes(x=Condition_Dry, y=prop_fixed)) +
  geom_boxplot()+
  ylab("Proportion hatched") + facet_wrap("Species") +
  xlab("Condition") + theme_bw() + ylim(-.0,1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),   
        axis.text.y  = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  geom_signif(comparisons = list(c("Cold_Dry", "Cold_Wet"), 
                                 c("Cold_Wet", "Warm_Dry"), 
                                 c("Cold_Wet", "Warm_Wet")),
              test = "wilcox.test", 
              test.args = list(exact = FALSE),
              map_signif_level = TRUE,
              y_position = c(0.8, 0.85, 0.9))
allcond_box2

allplots <- (
  popbox + popplot +
    eppbox + numofeppplot +
    malebox + maledataplot +
    fillbox + fillpoint +
    earlybox2 + earlyprop2 +
    allcond_box2 + allcond_plot2
) + plot_layout(ncol = 2, widths = c(1, 2))

ggsave("C:/Users/rjpor/Downloads/allplots.png", plot = allplots, width = 12, height = 16, dpi = 300)
ggsave("C:/Users/rjpor/Downloads/allplots.pdf", plot = allplots, width = 12, height = 16, dpi = 300)






beforeplots <- (
  popbox + popplot +
    malebox + maledataplot +
    eppbox + numofeppplot +
    fillbox + fillpoint
) + plot_layout(ncol = 2, widths = c(1, 2))

duringplots <- (  earlybox2 + earlyprop2)+ plot_layout(ncol = 2, widths = c(1, 2))

hatchplots <- (allcond_box2 + allcond_plot2)+ plot_layout(ncol = 2, widths = c(1, 2))





# Add subplot labels: A, B, C...
beforeplots <- (
  popbox + popplot +
    malebox + maledataplot +
    eppbox + numofeppplot +
    fillbox + fillpoint
) + 
  plot_layout(ncol = 2, widths = c(1, 5)) +
  plot_annotation(title = "", tag_levels = "A")

duringplots <- (
  earlybox2 + earlyprop2
) + 
  plot_layout(ncol = 2, widths = c(1, 5)) +
  plot_annotation(title = "", tag_levels = "A")

hatchplots <- (
  allcond_box2 + allcond_plot2
) + 
  plot_layout(ncol = 2, widths = c(1, 5)) +
  plot_annotation(title = "", tag_levels = "A")


ggsave("C:/Users/rjpor/Downloads/beforeplots.pdf", plot = beforeplots, width = 12, height = 12, dpi = 300)
ggsave("C:/Users/rjpor/Downloads/duringplots.pdf", plot = duringplots, width = 12, height = 4, dpi = 300)
ggsave("C:/Users/rjpor/Downloads/hatchplots.pdf", plot = hatchplots, width = 12, height = 4, dpi = 300)





pdf("C:\\Users\\rjpor\\Downloads\\allplots.pdf", width=15, height=20)

allplots + ggtitle("allplots")

dev.off()



#Anovas

#population size

popglm <- glmer(Indiv.L~Species+(1|pondclone ) + (1|pondclone:Rep), maledata, family=poisson())
popglm2 <- glmer(Indiv.L~1+(1|pondclone ) + (1|pondclone:Rep), maledata, family=poisson())
summary(popglm)
anova(popglm,popglm2)

maledatapulex <- subset(maledata, Species == "D. pulex")
maledataobtusa <- subset(maledata, Species == "D. obtusa")

poppulexglm <- glmer(Indiv.L~1 + (1|pondclone) + (1|pondclone:Rep), maledatapulex, family=poisson())
popglm_noclonepulex <- glmer(Indiv.L ~ 1 + (1 | pondclone:Rep),
                        data = maledatapulex,
                        family = poisson())
anova(poppulexglm,popglm_noclonepulex)


popobtusaglm <- glmer(Indiv.L~1 + (1|pondclone) + (1|pondclone:Rep), maledataobtusa, family=poisson())
popglm_nocloneobtusa <- glmer(Indiv.L ~ 1 + (1 | pondclone:Rep),
                             data = maledataobtusa,
                             family = poisson())

anova(popobtusaglm,popglm_nocloneobtusa)








#number of ephippia

#heritability within species as the variance effect of clone

numofeppglm <- glmer(N~Species+(1|plate_clone) + (1|plate_clone:Replicate), NumofEpp2, family=poisson())
numofeppglm2 <- glmer(N~1+(1|plate_clone) + (1|plate_clone:Replicate), NumofEpp2, family=poisson())
summary(numofeppglm)
anova(numofeppglm,numofeppglm2)


NumofEpp2pulex <- subset(NumofEpp2, Species == "D_pulex")
NumofEpp2obtusa <- subset(NumofEpp2, Species == "D_obtusa")

NumofEpp2pulexglm <- glmer(N~1+(1|plate_clone) + (1|plate_clone:Replicate), NumofEpp2pulex, family=poisson())
NumofEpp2pulexglm2 <- glmer(N~1 + (1|plate_clone:Replicate), NumofEpp2pulex, family=poisson())
summary(NumofEpp2pulexglm)
anova(NumofEpp2pulexglm,NumofEpp2pulexglm2)

NumofEpp2obtusaglm <- glmer(N~1+(1|plate_clone) + (1|plate_clone:Replicate), NumofEpp2obtusa, family=poisson())
NumofEpp2obtusaglm2 <- glmer(N~1 + (1|plate_clone:Replicate), NumofEpp2obtusa, family=poisson())
summary(NumofEpp2obtusaglm)
anova(NumofEpp2obtusaglm,NumofEpp2obtusaglm2)



overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp  <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  ratio <- Pearson.chisq / rdf
  pval  <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  
  data.frame(
    chisq = Pearson.chisq,
    ratio = ratio,
    rdf = rdf,
    p = pval
  )
}

overdisp_fun(numofeppglm)
overdisp_fun(numofeppglm2)

overdisp_fun(NumofEpp2pulexglm)
overdisp_fun(NumofEpp2pulexglm2)

overdisp_fun(NumofEpp2obtusaglm)
overdisp_fun(NumofEpp2obtusaglm2)

testZeroInflation(sim)



library(DHARMa)

sim <- simulateResiduals(numofeppglm)
testDispersion(sim)
plot(sim)

#chi2 DF Pr in text
#report "23.309  1   1.38e-06 ***" 



#percent male



malepercentglm <- glmer(Males.L/Indiv.L~Species+(1|pondclone ) + (1|pondclone:Rep), maledata, family=binomial(), weights=Indiv.L)
malepercentglm2 <- glmer(Males.L/Indiv.L~1+(1|pondclone ) + (1|pondclone:Rep), maledata, family=binomial(), weights=Indiv.L)
summary(malepercentglm)
anova(malepercentglm,malepercentglm2)




malepercentglmpulex <- glmer(Males.L/Indiv.L~1 + (1|pondclone) + (1|pondclone:Rep), maledatapulex, family=binomial(), weights=Indiv.L)
malepercentglm_noclonepulex <- glmer(Males.L/Indiv.L ~ 1 + (1 | pondclone:Rep),
                             data = maledatapulex,
                             family=binomial(), weights=Indiv.L)
anova(malepercentglmpulex,malepercentglm_noclonepulex)


malepercentglmobtusa <- glmer(Males.L/Indiv.L~1 + (1|pondclone) + (1|pondclone:Rep), maledataobtusa, family=binomial(), weights=Indiv.L)
malepercent_nocloneobtusa <- glmer(Males.L/Indiv.L ~ 1 + (1 | pondclone:Rep),
                              data = maledataobtusa,
                              family=binomial(), weights=Indiv.L)

anova(malepercentglmobtusa,malepercent_nocloneobtusa)
summary(malepercentglmobtusa)




#fill rate

fillrateglm <- glmer(totalemb/putativetotal~Species+(1|plate_clone ) + (1|plate_clone:Replicate), Dryepps_withspecies, family=binomial(), weights=putativetotal)
fillrateglm2 <- glmer(totalemb/putativetotal~1+(1|plate_clone ) + (1|plate_clone:Replicate), Dryepps_withspecies, family=binomial(), weights=putativetotal)
summary(fillrateglm)
anova(fillrateglm,fillrateglm2)




Dryepps_withspeciespulex <- subset(Dryepps_withspecies, Species == "D_pulex")
Dryepps_withspeciesobtusa <- subset(Dryepps_withspecies, Species == "D_obtusa")


Dryepps_withspeciespulexglm <- glmer(totalemb/putativetotal~1+(1|plate_clone ) + (1|plate_clone:Replicate), Dryepps_withspeciespulex, family=binomial(), weights=putativetotal)
Dryepps_withspeciespulexglm2 <- glmer(totalemb/putativetotal~1 + (1|plate_clone:Replicate), Dryepps_withspeciespulex, family=binomial(), weights=putativetotal)
summary(Dryepps_withspeciespulexglm)
anova(Dryepps_withspeciespulexglm,Dryepps_withspeciespulexglm2)

Dryepps_withspeciesobtusaglm <- glmer(totalemb/putativetotal~1+(1|plate_clone ) + (1|plate_clone:Replicate), Dryepps_withspeciesobtusa, family=binomial(), weights=putativetotal)
Dryepps_withspeciesobtusaglm2 <- glmer(totalemb/putativetotal~1 + (1|plate_clone:Replicate), Dryepps_withspeciesobtusa, family=binomial(), weights=putativetotal)
summary(Dryepps_withspeciesobtusaglm)
anova(Dryepps_withspeciesobtusaglm,Dryepps_withspeciesobtusaglm2)

#early rate


earlyprop2glm <- glmer(N.y/realtotal~Species+(1|plate_clone ) + (1|plate_clone:Replicate), mergedonlyhatch, family=binomial(), weights=realtotal)
earlyprop2glm2 <- glmer(N.y/realtotal~1+(1|plate_clone ) + (1|plate_clone:Replicate), mergedonlyhatch, family=binomial(), weights=realtotal)
summary(earlyprop2glm)
anova(earlyprop2glm,earlyprop2glm2)

earlytotalspulex <- subset(mergedonlyhatch, Species == "D_pulex")
earlytotalsobtusa <- subset(mergedonlyhatch, Species == "D_obtusa")

earlyproppulexglm <- glmer(N.y/realtotal~1+(1|plate_clone ) + (1|plate_clone:Replicate), earlytotalspulex, family=binomial(), weights=realtotal)
earlyproppulexglm2 <- glmer(N.y/realtotal~1 + (1|plate_clone:Replicate), earlytotalspulex, family=binomial(), weights=realtotal)
summary(earlyproppulexglm)
anova(earlyproppulexglm,earlyproppulexglm2)

earlypropobtusaglm <- glmer(N.y/realtotal~1+(1|plate_clone ) + (1|plate_clone:Replicate), earlytotalsobtusa, family=binomial(), weights=realtotal)
earlypropobtusaglm2 <- glmer(N.y/realtotal~1 + (1|plate_clone:Replicate), earlytotalsobtusa, family=binomial(), weights=realtotal)
summary(earlypropobtusaglm)
anova(earlypropobtusaglm,earlypropobtusaglm2)

#differences between conditions

allcond_complete_data <- all_cond_fill %>%
  dplyr::filter(!is.na(prop_fixed), 
                !is.na(Species), 
                !is.na(plate_clone), 
                !is.na(Replicate), 
                !is.na(total))

allcondglm <- glmer(prop_fixed~Species+(1|plate_clone ) + (1|plate_clone:Replicate), allcond_complete_data, family=binomial(), weights=realtotal)
allcondglm2 <- glmer(prop_fixed~1+(1|plate_clone ) + (1|plate_clone:Replicate), allcond_complete_data, family=binomial(), weights=realtotal)
summary(allcondglm)
anova(allcondglm,allcondglm2)








all_cond_fillsep <- all_cond_fill %>%
  separate(Condition_Dry, into = c("Condition", "Dry"), sep = "_")

allcondglm <- glmer(prop_fixed~Species*Condition*Dry+(1|plate_clone) + (1|plate_clone:Replicate), all_cond_fillsep, family=binomial(), weights=realtotal)


allcondglm <- glmer(prop_fixed~Species*Dry+(1|plate_clone) + (1|plate_clone:Replicate), all_cond_fillsep, family=binomial(), weights=realtotal)
allcondglm2 <- glmer(prop_fixed~Species*1+(1|plate_clone) + (1|plate_clone:Replicate), all_cond_fillsep, family=binomial(), weights=realtotal)
anova(allcondglm,allcondglm2)






all_cond_fillsepobtusa <- subset(all_cond_fillsep, Species == "D_obtusa")
allcondglm <- glmer(prop_fixed~Condition*Dry+(1|plate_clone) + (1|plate_clone:Replicate), all_cond_fillsepobtusa, family=binomial(), weights=realtotal)
allcondglm2 <- glmer(prop_fixed~1+(1|plate_clone) + (1|plate_clone:Replicate), all_cond_fillsepobtusa, family=binomial(), weights=realtotal)
anova(allcondglm,allcondglm2)

allcondglm <- glmer(prop_fixed~Condition*Dry+(1|plate_clone) + (1|plate_clone:Replicate), all_cond_fillsepobtusa, family=binomial(), weights=realtotal)
allcondglm2 <- glmer(prop_fixed~Condition+(1|plate_clone) + (1|plate_clone:Replicate), all_cond_fillsepobtusa, family=binomial(), weights=realtotal)
anova(allcondglm,allcondglm2)


allcondglm <- glmer(prop_fixed~Condition+(1|plate_clone) + (1|plate_clone:Replicate), all_cond_fillsepobtusa, family=binomial(), weights=realtotal)
allcondglm2 <- glmer(prop_fixed~1+(1|plate_clone) + (1|plate_clone:Replicate), all_cond_fillsepobtusa, family=binomial(), weights=realtotal)
anova(allcondglm,allcondglm2)





allcondglm <- glmer(prop_fixed~Dry+(1|plate_clone) + (1|plate_clone:Replicate), all_cond_fillsepobtusa, family=binomial(), weights=realtotal)
allcondglm2 <- glmer(prop_fixed~1+(1|plate_clone) + (1|plate_clone:Replicate), all_cond_fillsepobtusa, family=binomial(), weights=realtotal)
anova(allcondglm,allcondglm2)


all_cond_fillseppulex <- subset(all_cond_fillsep, Species == "D_pulex")
allcondglm <- glmer(prop_fixed~Condition*Dry+(1|plate_clone) + (1|plate_clone:Replicate), all_cond_fillseppulex, family=binomial(), weights=realtotal)
allcondglm2 <- glmer(prop_fixed~1+(1|plate_clone) + (1|plate_clone:Replicate), all_cond_fillseppulex, family=binomial(), weights=realtotal)
anova(allcondglm,allcondglm2)

allcondglm <- glmer(prop_fixed~Dry+(1|plate_clone) + (1|plate_clone:Replicate), all_cond_fillseppulex, family=binomial(), weights=realtotal)
allcondglm2 <- glmer(prop_fixed~1+(1|plate_clone) + (1|plate_clone:Replicate), all_cond_fillseppulex, family=binomial(), weights=realtotal)
anova(allcondglm,allcondglm2)








Coldwetsep <- subset(all_cond_fillsep, Condition=="Cold" & Dry=="Wet" & Species == "D_obtusa")
Coldwetsepglm <- glmer(prop_fixed~1+(1|plate_clone ) + (1|plate_clone:Replicate), Coldwetsep, family=binomial(), weights=realtotal)
Coldwetsepglm2 <- glmer(prop_fixed~1+(1|plate_clone:Replicate), Coldwetsep, family=binomial(), weights=realtotal)
summary(Coldwetsepglm)
anova(Coldwetsepglm,Coldwetsepglm2)


Coldwetsep <- subset(all_cond_fillsep, Condition=="Cold" & Dry=="Wet" & Species == "D_pulex")
Coldwetsepglm <- glmer(prop_fixed~1+(1|plate_clone ) + (1|plate_clone:Replicate), Coldwetsep, family=binomial(), weights=realtotal)
Coldwetsepglm2 <- glmer(prop_fixed~1+(1|plate_clone:Replicate), Coldwetsep, family=binomial(), weights=realtotal)
summary(Coldwetsepglm)
anova(Coldwetsepglm,Coldwetsepglm2)


#diversified hatching

allcond_plot<- ggplot(data=all_conds, aes(x=Light, y=prop, col=Temp)) +
  geom_point(size=2, position = position_dodge(width = 0.5))+
  ylab("Proportion hatched") + facet_wrap("Dried") +
  geom_errorbar(aes(ymin=uci, ymax=lci), width=.2, position = position_dodge(width = 0.5)) +
  xlab("Condition") + theme_bw() + ylim(-.1,1) +
  theme(axis.text.x  = element_text(size = 14),   
        axis.text.y  = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
allcond_plot

ggsave("C:/Users/rjpor/Downloads/allcond_plot.png", plot = allcond_plot, width = 6, height = 6, dpi = 300)
ggsave("C:/Users/rjpor/Downloads/allcond_plot.pdf", plot = allcond_plot, width = 6, height = 6, dpi = 300)


pairwise_results <- list()

for (i in 1:(nrow(all_conds)-1)) {
  for (j in (i+1):nrow(all_conds)) {
    
    x1 <- all_conds$N[i]
    n1 <- all_conds$fixed[i]
    
    x2 <- all_conds$N[j]
    n2 <- all_conds$fixed[j]
    
    test <- prop.test(
      x = c(x1, x2),
      n = c(n1, n2),
      correct = FALSE
    )
    
    pairwise_results[[paste(i, j, sep = "_")]] <- list(
      cond1 = all_conds[i, .(Temp, Dried, Light)],
      cond2 = all_conds[j, .(Temp, Dried, Light)],
      Xsq   = unname(test$statistic),
      pval  = test$p.value
    )
  }
}

pairwise_results

pairtab <- rbindlist(lapply(names(pairwise_results), function(k) {
  x <- pairwise_results[[k]]
  
  data.table(
    comparison = k,
    
    # Condition 1
    Temp1  = x$cond1$Temp,
    Dried1 = x$cond1$Dried,
    Light1 = x$cond1$Light,
    
    # Condition 2
    Temp2  = x$cond2$Temp,
    Dried2 = x$cond2$Dried,
    Light2 = x$cond2$Light,
    
    # Stats
    Xsq  = x$Xsq,
    pval = x$pval
  )
}))

pairtab


target <- all_conds[Temp=="Cold" & Dried=="Dry" & Light=="Light"]

other_conds <- all_conds[!(Temp=="Cold" & Dried=="Dry" & Light=="Light")]

results <- other_conds[, {
  test <- prop.test(
    x = c(target$N, N),
    n = c(target$fixed, fixed),
    correct = FALSE
  )
  list(
    compare_to = paste(Temp, Dried, Light, sep="-"),
    Xsq = as.numeric(test$statistic),
    pval = test$p.value
  )
}, by = 1:nrow(other_conds)]

results





focal <- all_conds[Temp=="Cold" & Dried=="Wet" & Light=="Light"]

# Pool all other treatments
others <- all_conds[!(Temp=="Cold" & Dried=="Wet" & Light=="Light")]

# Sum counts for others
x_other <- sum(others$N)       # total hatched in all other treatments
n_other <- sum(others$fixed)   # total effective denominator in all other treatments

# Values for focal treatment
x_focal <- focal$N
n_focal <- focal$fixed

# Two-sample proportions test
test_result <- prop.test(
  x = c(x_focal, x_other),
  n = c(n_focal, n_other),
  correct = FALSE
)

test_result


prop.test(x = 190, n = 252, correct = FALSE)





prop.test(x = 190, n = 251.6823, correct = FALSE)
prop.test(x = 108, n = 250.3646 , correct = FALSE)



allcondglm <- glm(N/fixed~Temp*Dried*Light, all_conds, family=binomial(), weights=fixed)

anova(allcondglm, test = "Chisq")

#Time series


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


timeseriesshort$se <- sqrt(timeseriesshort$prop * (1-timeseriesshort$prop)/ (2*72*0.6588542))
timeseriesshort$lci <- timeseriesshort$prop - (1.96 * timeseriesshort$se)
timeseriesshort$uci <- timeseriesshort$prop + (1.96 * timeseriesshort$se)



#try weighting by fill rate per replicate, rather than overall?
pdf("C:\\Users\\rjpor\\Downloads\\time_series_prop.pdf", width=4, height=3)

ggplot(data=timeseriesshort, aes(x=Time, y=prop)) +
  geom_point(size=2)+ ggtitle("Proportion of Hatchlings Over Time")+
  ylab("Proportion of Hatchlings") + 
  xlab("Time in Cold Treatment (Weeks)") + theme_bw() + ylim(0,0.6)+
  geom_errorbar(aes(ymin=uci, ymax=lci), width=.2)

dev.off()



proptime <- ggplot(data=timeseriesshort, aes(x=Time, y=prop)) +
  geom_point(size=2)+ ggtitle("Proportion of Hatchlings Over Time")+
  ylab("Proportion of Hatchlings") + 
  xlab("Time in Cold/Dark Treatment (Weeks)") + theme_bw() + ylim(0,0.6)+
  geom_errorbar(aes(ymin=uci, ymax=lci), width=.2)

ggsave("C:\\Users\\rjpor\\Downloads\\time_series_prop.png", proptime, width = 4, height = 4, dpi = 600 )




 # stat_smooth(method = "lm", formula = y ~ x+ I(x^2), se= F)+
#  stat_fit_glance(method = 'lm', label.x="middle", label.y="bottom",
#                  method.args = list(formula =  y ~ x+ I(x^2)),
#                  geom = 'text',
#                  aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
#                  size = 3)

pdf("C:\\Users\\rjpor\\Downloads\\time_series.pdf", width=4, height=3)

ggplot(data=timeseriesshort, aes(x=Time, y=N)) +
  geom_point(size=2)+ ggtitle("Number of Hatchlings Over Time")+
  ylab("Number of Hatchlings (out of total embryos)") + 
  xlab("Time in Cold Treatment (Weeks)") + theme_bw() + ylim(0,45)+ 
  stat_smooth(method = "lm", formula = y ~ x+ I(x^2), se= F)+
  stat_fit_glance(method = 'lm', label.x="middle", label.y="bottom",
                  method.args = list(formula =  y ~ x+ I(x^2)),
                  geom = 'text',
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                  size = 3)

dev.off()




library(ggplot2)
library(maps)

us_states <- map_data("state")

p <- ggplot(us_states, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "white", color = "black", linewidth = 0.3) +
  coord_fixed(
    ratio = 1.3,              # correct lat/long scaling for US maps
    xlim = c(-95, -66),       # crop region
    ylim = c(25, 48)
  ) +
  theme_void()

# save as square
ggsave("C:\\Users\\rjpor\\Downloads\\east_midwest_map.png", p,
       width = 6, height = 6, dpi = 300)
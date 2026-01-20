### This script contains the data analyses of water holding traits and hydrophobicity ###

### Packages ###
library("lme4")
library("ggeffects")
require("ggplot2")
library("gridExtra")
library("MASS")
library("glmmTMB")
library("car")
library("sjPlot")
library("emmeans")
library("multcomp")
library("Rcpp")
library("broom.mixed")
library("dplyr")
library("kableExtra")

############# Data reading ##########

setwd("C:/R/tesis/trasplante") #set working directory, this is different for each user

DFt <- read.table("CRH_trasplante.csv",
                  sep = ";",
                  dec = ".",
                  header = T) 

DF <- read.table("hidrofobicidad.csv",
                 sep = ";",
                 dec = ".",
                 header = T) 


################## GLMMs ##############
head(DFt)

DFt$tmt <- factor(DFt$tmt, levels = c("control_baja", "baja_media", "control_sur",
                                      "sur_norte","control_norte","norte_sur","media_alta"), 
                  ordered = F) # ordering treatments


fithidrof=glmmTMB(tinicio ~ tmt + status_initial +(1|roca) + (1|spp),
                  family = nbinom1,
                  data = DF) # model for hydrophobicity in which Tini was modeled as a response variable, tmt is a fixed effect, 
# while the status of the droplet is a covariable. Each rock and lichen species are modeled as random effects
summary(fit1)
Anova(fit1)

fit.WHCtot <- lmer(log(CRHtot) ~ tmto  + (1 | roca) + (1 | spp), data = DFt)
summary(fit.WHCtot)
Anova(fit.WHCtot)

fit.WHCint <- lmer(CRHint ~ tmto  + (1 | roca) + (1 | spp), data = DFt)
summary(fit.WHCint)
Anova(fit.WHCint)

fit.WHCext <- lmer(CRHext ~ tmto  + (1 | roca) + (1 | spp), data = DFt)
summary(fit.WHCext)
Anova(fit.WHCext)

## Model descriptions ##
# A model per water retention trait was carried in which each parameter was modeled as response variable, 
# and treatments was modeled as a fixed effect. Species and Rocks were treated as independent random effects  

### Estimating means from the carried models ###

A=ggemmeans(fithidrof, terms = c("tmt"), ci.lvl= 0.95) # estimated means for hydrophobicity

gA=plot(A, facets = F, colors = "bw", ci.style = "ribbon",connect.lines = F,
        show.legend = T, jitter=T)+
  theme_bw() + labs(x="Tratamiento",title = "A")+
  ylab(bquote("T"[ini] ~("s")))+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

B=ggemmeans(fit.WHCtot, terms = c("tmto"), ci.lvl=0.95) # estimated means for total water holding capacity

gB=plot(B, ci.style="ribbon", colors = "bw") + 
  theme_bw()+labs(x="Tratamiento", title = "B")+
  ylab(bquote(WHCtot ~(mg[H2O] ~cm^-2)))+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

C=ggemmeans(fit.WHCint, terms = c("tmto"), ci.lvl=0.95) # estimated means for internal water holding capacity
gC=plot(C, ci.style="ribbon", colors = "bw") + 
  theme_bw()+labs(x="Tratamiento", title = "C")+
  ylab(bquote(WHCint ~(mg[H2O] ~cm^-2)))+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

D=ggemmeans(fit.WHCext, terms = c("tmto"), ci.lvl=0.95) # estimated means for external water holding capacity
gD=plot(D, ci.style="ribbon", colors = "bw") + 
  theme_bw()+labs(x="Tratamiento", title = "D")+
  ylab(bquote(WHCext ~(mg[H2O] ~cm^-2)))+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

windows()
grid.arrange(gA, gB, gC, gD)

## Posthoc/a posteriori tests ##

ph.hidrof <- glht(fithidrof, emm(pairwise ~ tmt))
summary(ph.hidrof, test=adjusted(type="none"))

ph.CRHtot <- glht(fit.WHCtot, emm(pairwise ~ tmto))
summary(ph.CRHtot, test=adjusted(type="none"))

ph.CRHint <- glht(fit.WHCint, emm(pairwise ~ tmto))
summary(ph.CRHint, test=adjusted(type="none"))

ph.CRHext <- glht(fit.WHCext, emm(pairwise ~ tmto))
summary(ph.CRHext, test=adjusted(type="none"))

### Water retention and hydrophobicity traits analyses per species  ###

CRHu<-subset(DFt,spp=="Usnea_amblyoclada") # a subset for water retention traits of Usnea amblyoclada data
CRHp<-subset(DFt,spp=="Parmotrema_reticulatum") # subset for water retention traits of Parmotrema reticulatum data

head(DF)

hidro.u<-subset(DF,spp=="U.amblyoclada") # a subset for hydrophobicity of Usnea amblyoclada data
hidro.p<-subset(DF,spp=="P.reticulatum") # a subset for hydrophobicity of Parmotrema reticulatum data

CRHu$tmto <- factor(CRHu$tmto, levels = c("control_baja", 
                                          "baja_media", 
                                          "control_sur",
                                          "sur_norte",
                                          "control_norte",
                                          "norte_norte",
                                          "norte_sur",
                                          "media_alta"), 
                    ordered = F) #ordering treatment levels

# GLMMs for each species for Usnea amblyoclada data #

fit.hidr.u=glmmTMB(tinicio ~ tmt + status_initial +(1|roca),
                   family = nbinom1,
                   data = hidro.u)
summary(fit.hidr.u)
Anova(fit.hidr.u)

fit.WHCtot <- lmer(log(CRHtot) ~ tmto  + (1 | roca), data = CRHu)
summary(fit.WHCtot)
Anova(fit.WHCtot)

fit.WHCint <- lmer(CRHint ~ tmto  + (1 | roca), data = CRHu)
summary(fit.WHCint)
Anova(fit.WHCint)

fit.WHCext <- lmer(CRHext ~ tmto  + (1 | roca), data = CRHu)
summary(fit.WHCext)
Anova(fit.WHCext)

# Estimated means for Usnea amblyoclada #

A=ggemmeans(fit.hidr.u, terms = c("tmt"), ci.lvl= 0.95)
gA=plot(A, facets = F, colors = "bw", ci.style = "ribbon",connect.lines = F,
        show.legend = T, jitter=T)+
  theme_bw() + labs(x="Tratamiento",title = "A")+
  ylab(bquote("T"[ini] ~("s")))+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

B=ggemmeans(fit.WHCtot, terms = c("tmto"),
            ci.lvl=0.95)
summary(B)
gB=plot(B, ci.style="ribbon", colors = "bw") + 
  theme_bw()+labs(x="Tratamiento", title = "B")+
  ylab(bquote(WHCtot ~(mg[H2O] ~cm^-2)))+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

C=ggemmeans(fit.WHCint, terms = c("tmto"),
            ci.lvl=0.95)
gC=plot(C, ci.style="ribbon", colors = "bw") + 
  theme_bw()+labs(x="Tratamiento", title = "C")+
  ylab(bquote(WHCint ~(mg[H2O] ~cm^-2)))+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

D=ggemmeans(fit.WHCext, terms = c("tmto"),
            ci.lvl=0.95)
gD=plot(D, ci.style="ribbon", colors = "bw") + 
  theme_bw()+labs(x="Tratamiento", title = "D")+
  ylab(bquote(WHCext ~(mg[H2O] ~cm^-2)))+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

windows()
grid.arrange(gA, gB, gC, gD) # final plots

## Posthoc tests ##

ph.hidrof <- glht(fit.hidr.u, emm(pairwise ~ tmt))
summary(ph.hidrof, test=adjusted(type="none"))

ph.CRHtot <- glht(fit.WHCtot, emm(pairwise ~ tmto))
summary(ph.CRHtot, test=adjusted(type="none"))

ph.CRHint <- glht(fit.WHCint, emm(pairwise ~ tmto))
summary(ph.CRHint, test=adjusted(type="none"))

ph.CRHext <- glht(fit.WHCext, emm(pairwise ~ tmto))
summary(ph.CRHext, test=adjusted(type="none"))

####### GLMs for Parmotrema reticulatum data ############

CRHp$tmto <- factor(CRHp$tmto, levels = c("control_baja", 
                                          "baja_media", 
                                          "control_sur",
                                          "sur_norte",
                                          "control_norte",
                                          "norte_norte",
                                          "norte_sur",
                                          "media_alta"), 
                    ordered = F)

fit.hidr.p=glmmTMB(tinicio ~ tmt + status_initial +(1|roca),
                   family = nbinom1,
                   data = hidro.p) # Parmotrema reticulatum's Hydrophobicity model 
summary(fit.hidr.p)
Anova(fit.hidr.p)
irr_hidrof.p <- tidy(fit.hidr.p, conf.int = TRUE, exponentiate = TRUE, effects = "fixed")
irr_hidrof.p %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  rename(IRR = estimate, CI_lower = conf.low, CI_upper = conf.high)

fit.WHCtot <- lmer(log(CRHtot) ~ tmto  + (1 | roca), data = CRHp) # Parmotrema reticulatum's total water holding capacity model
summary(fit.WHCtot)
Anova(fit.WHCtot)

fit.WHCint <- lmer(CRHint ~ tmto  + (1 | roca), data = CRHp) # Parmotrema reticulatum's internal water holding capacity model
summary(fit.WHCint)
Anova(fit.WHCint)

fit.WHCext <- lmer(CRHext ~ tmto  + (1 | roca), data = CRHp) # Parmotrema reticulatum's external water holding capacity model
summary(fit.WHCext)
Anova(fit.WHCext)

# Estimated means for Parmotrema reticulatum's models #

A=ggemmeans(fit.hidr.p, terms = c("tmt"), ci.lvl= 0.95)
gA=plot(A, facets = F, colors = "bw", ci.style = "ribbon",connect.lines = F,
        show.legend = T, jitter=T)+
  theme_bw() + labs(x="Tratamiento",title = "A")+
  ylab(bquote("T"[ini] ~("s")))+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

B=ggemmeans(fit.WHCtot, terms = c("tmto"),
            ci.lvl=0.95)
gB=plot(B, ci.style="ribbon", colors = "bw") + 
  theme_bw()+labs(x="Tratamiento", title = "B")+
  ylab(bquote(CRHtot ~(mg[H2O] ~cm^-2)))+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

C=ggemmeans(fit.WHCint, terms = c("tmto"),
            ci.lvl=0.95)
gC=plot(C, ci.style="ribbon", colors = "bw") + 
  theme_bw()+labs(x="Tratamiento", title = "C")+
  ylab(bquote(CRHint ~(mg[H2O] ~cm^-2)))+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

D=ggemmeans(fit.WHCext, terms = c("tmto"),
            ci.lvl=0.95)
gD=plot(D, ci.style="ribbon", colors = "bw") + 
  theme_bw()+labs(x="Tratamiento", title = "D")+
  ylab(bquote(CRHext ~(mg[H2O] ~cm^-2)))+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

windows()
grid.arrange(gA, gB, gC, gD) # Final plots

## posthoc ##

ph.hidrof <- glht(fit.hidr.p, emm(pairwise ~ tmt))
summary(ph.hidrof, test=adjusted(type="none"))

ph.CRHtot <- glht(fit.WHCtot, emm(pairwise ~ tmto))
summary(ph.CRHtot, test=adjusted(type="none"))

ph.CRHint <- glht(fit.WHCint, emm(pairwise ~ tmto))
summary(ph.CRHint, test=adjusted(type="none"))

ph.CRHext <- glht(fit.WHCext, emm(pairwise ~ tmto))
summary(ph.CRHext, test=adjusted(type="none"))

##### Hydrophobicity curves #

setwd("C:/R/tesis/trasplante") #optional
require("GGally")
library("survival")
library("gridExtra")

DF <- read.table("hidrofobicidad.csv",
                 sep = ";",
                 dec = ".",
                 header = T)

head(DF)

altura <- DF[DF$tmt %in% c("control_baja", "baja_media", "control_sur", "media_alta"), ] #creating subsets for elevation data
micrositio <- DF[DF$tmt %in% c("control_norte", "control_sur", "norte_sur", "sur_norte"), ] #subset for microsite data

# Subsets for species
usnea <- DF[DF$spp %in% c("U.amblyoclada"), ]
usnea.alt <- usnea[usnea$tmt %in% c("control_baja", "baja_media", "control_sur", "media_alta"), ]
usnea.micro <- usnea[usnea$tmt %in% c("control_norte", "control_sur", "norte_sur", "sur_norte"), ]

parmotrema <- DF[DF$spp %in% c("P.reticulatum"), ]
parmotrema.alt <- parmotrema[parmotrema$tmt %in% c("control_baja", "baja_media", "control_sur", "media_alta"), ]
parmotrema.micro <- parmotrema[parmotrema$tmt %in% c("control_norte", "control_sur", "norte_sur", "sur_norte"), ]

# Fit survival model
fit1 <- survfit(Surv(tini, status_initial) ~ tmt, data = DF)
fit2 <- survfit(Surv(tini, status_initial) ~ tmt, data = altura)
fit3 <- survfit(Surv(tini, status_initial) ~ tmt, data = micrositio)

g1 <- ggsurv(fit1, back.white = T, CI=F, 
             xlab="Tiempo de absorción (s)", ylab="Estado de absorción",
             main="A")

g2 <- ggsurv(fit2, back.white = T, CI=F, 
             xlab="Tiempo de absorción (s)", ylab="Estado de absorción",
             main="A")

g3 <- ggsurv(fit3, back.white = T, CI=F, 
             xlab="Tiempo de absorción (s)", ylab="Estado de absorción",
             main="B")

fit4 <- survfit(Surv(tini, status_initial) ~ tmt, data = usnea.alt)
fit5 <- survfit(Surv(tini, status_initial) ~ tmt, data = usnea.micro)

fit6 <- survfit(Surv(tini, status_initial) ~ tmt, data = parmotrema.alt)
fit7 <- survfit(Surv(tini, status_initial) ~ tmt, data = parmotrema.micro)

g4 <- ggsurv(fit4, back.white = T, CI=F, 
             xlab="Tiempo de absorción (s)", ylab="Estado de absorción",
             main="A")

g5 <- ggsurv(fit5, back.white = T, CI=F, 
             xlab="Tiempo de absorción (s)", ylab="Estado de absorción",
             main="B")

g6 <- ggsurv(fit6, back.white = T, CI=F, 
             xlab="Tiempo de absorción (s)", ylab="Estado de absorción",
             main="C")

g7 <- ggsurv(fit7, back.white = T, CI=F, 
             xlab="Tiempo de absorción (s)", ylab="Estado de absorción",
             main="D")

windows()
grid.arrange(g4,g6,g5,g7) # final plots

########################################################################################################################
################### Model for usnic acid concentration of Usnea amblyoclada samples ####################################

setwd("C:/R/tesis/trasplante") # optional

DFus <- read.table("usnico_nov.csv",
                   sep = ";",
                   dec = ".",
                   header = T) #reading data

# GLMs for Usnic acid concentration #

m.usnico <- glm(formula = usnico ~ tmto, family = gaussian(link = "log"), 
                data = DFus)
summary(m.usnico)
Anova(m.usnico)

head(DFus)

# GLM for extract yield EY #

m.yield <- glmmTMB(rendimiento ~ tmto, 
                   DFus, family = "beta_family")
summary(m.yield)
Anova(m.yield)

# Plots #

DFus$tmto <- factor(DFus$tmto, levels = c("control_baja", 
                                          "baja_media", 
                                          "control_sur",
                                          "sur_norte",
                                          "control_norte",
                                          "norte_norte",
                                          "norte_sur",
                                          "media_alta"), 
                    ordered = F)
head(DFus)

gA <- ggplot(DFus , aes(x=tmto, y=usnico))+ 
  geom_boxplot()+  
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(x="Tratamiento", y= "AU (ppm)", title = "A")

gB <- ggplot(DFus , aes(x=tmto, y=rendimiento))+ 
  geom_boxplot()+  
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(x="Tratamiento", y=" rendimiento (%)", title = "B")

windows()
grid.arrange(gA, gB) # final plots

## Posthoc for Usnic acid concentration and EY models ##

ph.AU <- glht(m.usnico, emm(pairwise ~ tmto))
summary(ph.AU, test=adjusted(type="none"))

ph.yield <- glht(m.yield, emm(pairwise ~ tmto))
summary(ph.yield, test=adjusted(type="none"))

### Photosynthetic yield related parameters ###
###############################################################################################3
############## FVFM ##########################################################################

rm(list=ls()) # optional

library("lattice")
library("glmmTMB")
library("lme4")
library("car")
library("multcomp")
library("emmeans")
library("gridExtra")
library("ggeffects")
library("ggplot2")

setwd("C:/R/tesis/trasplante")
datos <- read.table("fvfm.csv",
                    sep = ";",
                    dec = ".",
                    header = T)

# subset per species #

usnea <- subset(datos,spp=="Usnea_amblyoclada")
parmotrema <- hidro.p<-subset(datos,spp=="Parmotrema_reticulatum")

m1 = lmer(FvFm ~ tmto + (1|roca), data=usnea)
summary(m1)
Anova(m1)

#Fo#
m2 = lmer(Fo ~ tmto + (1|roca), data=usnea)
summary(m2)
Anova(m2)

m3 = lmer(Fm ~ tmto + (1|roca), data=usnea)
summary(m3)
Anova(m3)

#reticulatum
m1.2 = lmer(FvFm ~ tmto+(1|roca), data=parmotrema)
summary(m1.2)
Anova(m1.2)

#########
m2 = lmer(Fo ~ tmto + (1|roca), data=parmotrema)
summary(m2)
Anova(m2)

#############
m3 = lmer(Fm ~ tmto + (1|muestra), data=parmotrema)
summary(m3)
Anova(m3)

# none sifnificantive differences were founded #

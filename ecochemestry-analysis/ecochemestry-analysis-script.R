# Loading packages

library("lattice")
library("ggeffects")
require("glmmTMB")
require("car")
library("fitdistrplus")
library("multcompView")
library("postHoc")
require("ggplot2")

# set wd (different for each user)
setwd("C:/....")
DFus <- read.table("cc_usnico_2.csv",
                  sep = ";",
                  dec = ".",
                  header = T)

# ensure that alt and muestra are factor variables
DFus$alt <- as.factor(DFus$alt) 
DFus$muestra <- as.factor(DFus$muestra)

# Descriptive statistics #

head(DFus)

se <- function(x) sqrt(var(x)/length(x)) # a function to calculate standard error

summary(DFus$usnico)
se(DFus$usnico)

# Max and min values of usnic acid concentration
DFus[DFus$usnico == max(DFus$usnico),]
DFus[which.min(DFus$usnico),]

summary(DFus$rendimiento)
se(DFus$rendimiento)
summary(DFus$peso)
se(DFus$peso)

### average, max and min values per elevation, orientation ### 
library("tidyr")
library("dplyr")
se <- function(x) sqrt(var(x)/length(x))

head(DFus)

cc <- DFus %>%
  group_by(alt)%>%
  summarise(med= mean(usnico),
            desv=se(usnico),
            mim=min(usnico),
            mam=max(usnico))
22.6/16.4
22.6/19

16.4*1.38

rr <- DFus %>%
  group_by(alt, orientacion)%>%
  summarise(med=mean(rendimiento),
            desv=se(rendimiento))

pp <- DFus %>%
  group_by(alt, orientacion)%>%
  summarise(med=mean(peso_liquen),
            desv=se(peso_liquen))

amb <- DFus %>%
  group_by(alt)%>%
  summarise(med= mean(exp),
            desv=se(exp),
            mim=min(exp),
            mam=max(exp))

################################################
# Usnic acid concentration per elevation plots #

head(DFus)
windows()
ggplot(DFus, aes(x=alt, y=usnico))+
  geom_boxplot()+  
  labs(x="Elevation (m.a.s.l.)",y="Usnic acid concentration (ppm)")+ 
  theme_bw()

windows()
ggplot(DFus, aes(x=orientacion, y=usnico))+
  geom_boxplot()+  
  labs(x="Elevation (m.a.s.l.)",y="Usnic acid concentration (ppm)")+ 
  theme_bw()

ggplot(DFus, aes(x=alt, y=usnico))+
  geom_boxplot()+  
  labs(x="Elevation (m.a.s.l.)",y="Usnic acid concentration (ppm)")+ 
  theme_bw()+
facet_wrap(~orientacion, scale="free")

windows()
ggplot(DFus, aes(x=orientacion, y=usnico))+
  geom_boxplot()+  
  labs(x="Elevation (m.a.s.l.)",y="Usnic acid concentration (ppm)")+ 
  theme_bw()+
  facet_wrap(~alt, scale="free")

windows()
ggplot(DFus, aes(x=alt, y=rendimiento))+
  geom_boxplot()+  
  labs(x="Elevation (m.a.s.l.)",y="Extract yield (%)")+ 
  theme_bw()

head(DFus)
windows()
ggplot(DFus, aes(x=orientacion, y=rendimiento))+
  geom_boxplot()+  
  labs(x="Rock aspect",y="Extract yield (%)")+ 
  theme_bw()

ggplot(DFus, aes(x=alt, y=rendimiento))+
  geom_boxplot()+  
  labs(x="Elevation (m.a.s.l.)",y="Extract yield (%)")+ 
  theme_bw()+facet_wrap(~orientacion)

windows()
ggplot(DFus, aes(x=orientacion, y=rendimiento))+
  geom_boxplot()+  
  labs(x="Microsite condition x Elevation (m.a.s.l.)",y="Extract yield (%)")+ 
  theme_bw()+facet_wrap(~alt)

head(DFus)
ggplot(DFus, aes(x=slope, y=usnico))+
  geom_boxplot()+
  labs(x="Microsite slope",y="Usnic acid concentration (ppm)")+ 
  theme_bw()

ggplot(DFus, aes(x=slope, y=rendimiento))+
  geom_boxplot()+
  labs(x="Microsite slope",y="Extract yield (%)")+ 
  theme_bw()

ggplot(DFus, aes(x=exp, y=usnico))+
  geom_point()+
  labs(x="Exposition",y="Extract yield (%)")+ 
  theme_bw()+geom_smooth(method=lm , color="red", se=T)

head(DFus)
ggplot(DFus, aes(x=orientacion, y=exp))+
  geom_point()+
  labs(x="Orientation",y="Exposition (°)")+ 
  theme_bw()+geom_smooth(method=lm , color="red", se=T)

head(DFus)
xyplot(exp ~ alt|orientacion, data = DFus)

### POST HOC ###
library("postHoc")
fit1 <- lm(usnico ~ alt + 0, data = DFus)
summary(fit1)
anova(fit1)
coef(fit1)

fitnull <- lm(usnico ~ 1, data = DFus)
anova(fit1, fitnull)

TT <- posthoc(Model = fit1, 
              EffectLabels = c("900","1800","2100"),
              digits = 1)
summary(TT)

print(TT)

windows()
barplot(TT, ylim = c(0,30),
        xlab="Elevation (m.a.s.l.)", ylab = "Usnic acid concentration (ppm)")
abline(h = 0)

head(DFus)
fit1 <- lm(usnico ~ orientacion + 0, data = DFus)
summary(fit1)
anova(fit1)
coef(fit1)

fitnull <- lm(usnico ~ 1, data = DFus)
anova(fit1, fitnull)

TT <- posthoc(Model = fit1, 
              EffectLabels = c("900","1800","2100"),
              digits = 1)
summary(TT)

print(TT)

windows()
barplot(TT, ylim = c(0,30),
        xlab="Pendiente (categoria)", ylab = "riqueza por roca")
abline(h = 0)

# tuckey en boxplot #
library(multcompView)
model=lm( usnico ~ alt, data = DFus )
ANOVA=aov(model)

TUKEY <- TukeyHSD(x=ANOVA, 'alt', conf.level=0.95)
plot(TUKEY , las=2 , col="brown")

generate_label_df <- function(TUKEY, variable){
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

LABELS <- generate_label_df(TUKEY , "alt")
LABELS$Letters <- factor(c("a","ab","b"))

windows()
a <- boxplot(DFus$usnico ~ DFus$alt, 
             ylim=c(min(DFus$usnico), 
                    1.1*max(DFus$usnico)), 
             ylab="Usnic acid concentration (ppm)", 
             xlab="Elevation (m.a.s.l.)", main="")


#text(3,31,expression("p > 0.01 ; R"^2*"= 99%"))

over <- 0.1*max(a$stats[nrow(a$stats),])

text(c(1:nlevels(DFus$alt)), 
     a$stats[nrow(a$stats),]+over, 
     LABELS[,1])

#ajustes con GLMs#

#cc usnico
fit1 <- lm(usnico ~ alt, data = DFus)
summary(fit1)
anova(fit1)
plot(fit1)

fit1.1 <- lm(usnico ~ orientacion + alt, data = DFus)
summary(fit1.1)
anova(fit1.1)

fit2 <- lm(usnico ~ alt*orientacion, data = DFus)
summary(fit2)
anova(fit2)

fit2.1 <- lm(usnico ~ orientacion*alt, data = DFus)
summary(fit2.1)
anova(fit2.1)

head(DFus)
windows()
ggplot(DFus, aes(x=alt, y=usnico))+
  geom_boxplot()+
  labs(x="Altitud (m.s.n.m.)",y="Concentración de ácido úsnico (ppm)", title = "")+ 
  theme_bw()

#######################################################################
################# Extract yield ###################
fit3 <- lm(rendimiento ~ alt + orientacion, data = DFus)
summary(fit3)
anova(fit3)

fit3.1 <- lm(rendimiento ~ orientacion + alt, data = DFus)
summary(fit3.1)
anova(fit3.1)

fit4<- lm(rendimiento ~ alt * orientacion, data = DFus)
summary(fit4)
anova(fit4)

#DFw <- na.omit(DFw)
head(DFus)
ggplot(DFus, aes(x=usnico))+
  geom_histogram(binwidth=0.9, fill="blue", colour="black")+  
  labs(x="cc Usnico total",y="cc Usnico acumulado")+ 
  theme_bw()+scale_y_discrete()

head(DFus)
require(fitdistrplus)
normal=fitdist(DFus$usnico,"norm")
lognormal=fitdist(DFus$usnico,"lnorm")
gama=fitdist(DFus$usnico,"gamma")

windows()
par(mfrow=c(1,2), 
    mai=c(1,1,0.25,0.25), 
    cex.lab=1.5, cex.axis=1.5, cex=1)

cdfcomp(list(normal,lognormal,gama),
        horizontals=F, ylogscale = T,
        xlogscale = F, 
        #lwd=2,
        addlegend=T,
        legendtext=c("normal","lognormal","gama"))
qqcomp(list(normal,lognormal,gama),
       addlegend=T,
       legendtext=c("normal","lognormal","gama"))
par(mfrow=c(1,1))

#modelo con gamma o normal#

s1.g <- glm(usnico ~ alt * orientacion, 
          data=DFus, 
          family=Gamma(link="log"))
summary(s1.g)

s1.n <- glm(usnico ~ alt * orientacion, 
            data=DFus, 
            family=gaussian(link="log"))
summary(s1.n)

windows()
par(mfrow=c(2,2))
resid.m1effS1.G=resid(s1.g,type="pearson")
resid.m1effS1.N=resid(s1.n,type="pearson")
plot(fitted(s1.g),resid.m1effS1.G); abline(h=0)
plot(fitted(s1.n),resid.m1effS1.N); abline(h=0)
qqnorm(resid(s1.g));qqline(resid(s1.g))
qqnorm(resid(s1.n));qqline(resid(s1.n))
par(mfrow=c(1,1))

anova(s1.n, test = "Chisq")

s2.1 <- update(s1.n, ~. -orientacion)
anova(s1.n, s2.1, test="Chisq")

anova(s2.1, test="Chisq")

s2.2 <- update(s2.1, ~. -alt:orientacion)
anova(s2.1, s2.2, test="Chisq")

summary(s2.2)
anova(s2.2, test="Chisq")

m.usnico <- glm(formula = usnico ~ alt, family = gaussian(link = "log"), 
                data = DFus)
summary(m.usnico)
Anova(m.usnico)

#########################################################################
windows()
# Obtain model residuals
residuals <- resid(m.usnico)

# Plot diagnostic plots
par(mfrow = c(2, 2))  # Set up a 2x2 grid for plots

# 1. Residuals vs Fitted
plot(fitted(m.usnico), residuals,
     main = "Residuals vs Fitted",
     xlab = "Fitted values",
     ylab = "Residuals")

# 2. Normal Q-Q Plot
qqnorm(residuals, main = "Normal Q-Q Plot")
qqline(residuals)

# 3. Scale-Location Plot (Square Root of Standardized Deviance Residuals vs Fitted)
plot(fitted(m.usnico), sqrt(abs(residuals)),
     main = "Scale-Location Plot",
     xlab = "Fitted values",
     ylab = "Square Root of Standardized Deviance Residuals")

# 4. Residuals vs Leverage
plot(hatvalues(m.usnico), residuals,
     main = "Residuals vs Leverage",
     xlab = "Leverage",
     ylab = "d")

#library(sjPlot)
#tab_model(m.usnico, show.loglik = T)

# modelo GLM for extract yield #

library(fitdistrplus)

DFus$rend <- DFus$rendimiento/100

beta=fitdist(DFus$rend,"beta")
normal=fitdist(DFus$rend,"norm")
gamma=fitdist(DFus$rend,"gamma")

windows()
par(mfrow=c(1,2), 
    mai=c(1,1,0.25,0.25), 
    cex.lab=1.5, cex.axis=1.5, cex=1)

cdfcomp(list(beta,normal,gamma),
        horizontals=F, ylogscale = T,
        xlogscale = F, 
        #lwd=2,
        addlegend=T,
        legendtext=c("beta","normal","gamma"))
qqcomp(list(beta,normal,gamma),
       addlegend=T,
       legendtext=c("beta","normal","gamma"))
par(mfrow=c(1,1))
gofstat(list(beta,normal,gamma),
        fitnames=c("beta", "normal", "gamma"))$bic
#beta oikoveta

DFus$alt <- as.factor(DFus$alt)

with(DFus, plot(density, rend,cex=2, pch=19, cex.lab=1.5, cex.axis=1.5))

require(glmmTMB)
require(car)

p.1 <- glmmTMB(rend ~ alt * orientacion + exp, DFus, family = "beta_family")
summary(p.1)
Anova(p.1)

p.2 <- update(p.1, ~. -exp)
anova(p.1, p.2)

summary(p.2)
Anova(p.2)
confint(p.2)

p.3 <- update(p.2, ~. -alt)
summary(p.3)
Anova(p.3)

head(DFus)

m.yield <- glmmTMB(rend ~ orientacion + alt*orientacion-alt, 
                   DFus, family = "beta_family")
summary(m.yield)
Anova(m.yield)

m10Pielou1<- m.yield

windows()
par(mfrow=c(2,2)) 
resid.m10Pielou=resid(m.yield,type="pearson")
plot(fitted(m.yield),resid.m10Pielou, 
     ylab="Residuals", 
     xlab="Predicted",main=""); abline(h=0)

plot(m.yield$frame$alt, resid.m10Pielou, 
     xlab="Elevation", ylab="Residuales",main="", las=3); abline(h=0)

qqnorm(resid.m10Pielou, ylab="Theorical quantiles", xlab="Theorical residuals quantiles",main="")
qqline(resid.m10Pielou,lwd=1.5, col="red")
par(mfrow=c(1,1),cex.lab=1, cex.axis=1)

b.m10=as.vector(unlist(fixef(m.yield)$cond)) # coefs con TODOS los datos
mcov=solve(vcov(m.yield)$cond) # inverse of matriz de var covarianza de fixefs
cook.m10=data.frame(dato=1:nrow(DFus), dist=NA) # dataframe para guardar las dist de Cook
for (i in 1:nrow(DFus)){
  m=glmmTMB(rend~orientacion+orientacion*alt-alt,data=DFus[-i,],family=beta_family)
  bi=as.vector(unlist(fixef(m)$cond)) # coefs sin fila i 
  cook.m10$dist[i]= t(b.m10-bi)%*% mcov %*% (b.m10-bi)
}

# Cook distances for model fitness diagnostics #
ggplot(data=cook.m10, aes(x=dato, y=dist))+ geom_point(size=1.5)+ 
  geom_linerange(aes(x=dato, ymin=0, ymax=dist))+ theme_bw()+
  labs(y="Distancia de Cook")
which(cook.m10$dist>0.5)

# m10Pielou output visualization:  curves
library(ggeffects)
# Curva predicha por modelo m10Pielou
pred.m.yield1=ggpredict(m.yield, terms = c("orientacion", "alt"))
pred.m.yield2=ggpredict(m.yield, terms = c("orientacion"))

# Final model
g1=plot(pred.m.yield1, facets = F)+ 
  theme_bw()+labs(x="Rock aspect", y="Predicted Extract yield (%)")+
  theme(legend.position = c(0.25,0.35), plot.title=element_blank())
g2=plot(pred.m.yield2)+ theme_bw()+labs(x="Orientacion", y="")+
  theme( plot.title=element_blank())

require(cowplot)
windows()
plot_grid(g1, g2, labels = c('A', 'B'), label_size = 12)

# Real data plots #

library("lattice")
head(DFus)
bwplot(rendimiento ~ alt|orientacion, data=DFus)

ggplot( DFus , aes(x=orientacion, y=rendimiento))+ 
  geom_boxplot()+  
  facet_wrap(~alt)+
  theme_bw()+
  theme(legend.position="none")

summary(m.yield)
Anova(m.yield)

require("emmeans")
require("multcomp")
means <- emmeans(m.yield, ~orientacion*alt)

pigs.emm <- emmeans(m.yield, "orientacion", type = "response")
posthoc <- multcomp::cld(means, alpha = 0.5, Letters = LETTERS)

posthoc$.group
mult.letters <- posthoc$.group

library(dplyr)
rend.summarized = DFus %>% 
  group_by(orientacion,alt) %>% 
  summarize(Max.rendimiento=max(rendimiento))

g1 <- ggplot(DFus, aes(x=orientacion, y=rendimiento))+ 
  geom_boxplot()+  
  facet_wrap(~alt)+
  theme_bw()+
  xlab("Microsite aspect*Elevation (m.a.s.l.)")+
  ylab("Extract yield (%)")

g2 <- g1 +
  geom_text(data=rend.summarized,
            aes(x=orientacion,
                y=0.2+Max.rendimiento,
                label=posthoc$.group),
            vjust=0)
windows()
plot(g2)

license
version
cite

ggplot( DFus , aes(x=orientacion, y=rendimiento))+ 
  geom_boxplot()+  
  #facet_wrap(~alt)+
  theme_bw()+
  theme(legend.position="none")

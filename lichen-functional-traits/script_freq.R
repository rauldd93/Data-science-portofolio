rm(list=ls())

setwd("C:/R/Data-science-portofolio/lichen-functional-traits")

library("ggeffects")
library("glmmTMB")
library("readr")
library("fitdistrplus")
library("lme4")
library("nlme")
library("DHARMa")
library("MuMIn")
library("car")
library("sjPlot")
library("ggplot2")
library("multcompView")
library("multcomp")
library("emmeans")
library("MASS")
library("pscl")
library("viridis")

# Tareas: agregar leyendas We tested for differences 
#dividing elevation a priori into three groups: low elevation (900–1900 m s.a.l), 
#medium elevation (1901–2400 m s.a.l.) and high elevation (2400–2680 m s.a.l.).

datos <- read.table("dat_freq.csv", 
                    sep = ";", 
                    dec = ",", 
                    header = T)

datos$Elevation<-as.factor(datos$Elevation)
#levels(datos$Elevation) <- c("Low (900 - 1900 mm s.a.l.)", "Medium (1901–2400 m s.a.l.)", "High (2400–2680 m s.a.l.)")
levels(datos$Elevation) <- c("Low", "Medium", "High")
datos$Slope<-as.factor(datos$Slope)
levels(datos$Slope) <- c("Low", "High")

#####################################################
#############folioso##############################
#####################################################
#####################################################

m2 <- glm(foliose ~ Elevation + Slope + Orientation, data = datos, family = quasipoisson(link = "log"))
summary(m2) #modelo para folioso
Anova(m2)

#Elev tukey and plot
comps <- glht(m2, linfct = mcp(Elevation = "Tukey"))
#summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Elevation = "Tukey")), Letters = letters)
y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Elevation=x,Pred = y)

g1 <- ggplot(df, aes(x = Elevation, y = Pred, group=Elevation)) +
  geom_boxplot() +
  labs(x = "Elevation", y = "Predicted frequency of foliose species", title = "A") +
  theme_bw()

tukey_df <- data.frame(
  Elevation = names(tukey_letters$mcletters$Letters),  # The Elevation levels ("0", "1", "2")
  Letters = tukey_letters$mcletters$Letters)            # The Tukey letters ("a", "a", "a")

plot.parameters <- scale_color_manual(values = c("Low" = "#40B0A6", "Medium" ="#E1BE6A", "High" = "#D55E00"))

#windows()
g.folioso.elev <- g1 + geom_text(data = tukey_df, 
                                 aes(x = Elevation, y = max(tukey_letters$lp) + 0.01, label = Letters), 
                                 vjust = -0.5, size = 5, color = "black")+
  geom_jitter(data= datos, aes(x=Elevation, y = foliose, color = Slope), width = 0.2, alpha=0.9) + 
  plot.parameters

g.folioso.elev

### Slope tukey and plot ###
comps <- glht(m2, linfct = mcp(Slope = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Slope = "Tukey")), Letters = letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Slope=x,Pred = y)

g1 <- ggplot(df, aes(x = Slope, y = Pred, group=Slope)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Slope", y = "Predicted richness of foliose species", title = "B") +
  theme_bw()

tukey_df <- data.frame(
  Slope = names(tukey_letters$mcletters$Letters),  
  Letters = tukey_letters$mcletters$Letters)   

g.folioso.slope <- g1 + geom_text(data = tukey_df, 
                                  aes(x = Slope, y = max(tukey_letters$lp) + 0.01, label = Letters), 
                                 vjust = -0.5, size = 5, color = "black")+ 
  geom_jitter(data=datos,aes(x=Slope, y = foliose, color = Elevation), 
              width = 0.2, alpha=0.9) +
  plot.parameters
  
#windows()
gridExtra::grid.arrange(g.folioso.elev, g.folioso.slope)

#####################################################
#############FructiculosoCWM#########################
#####################################################
#####################################################

m2 <- glm(fructicose ~ Elevation + Slope + Orientation, data = datos, family = quasipoisson(link = "log"))

summary(m2)
Anova(m2,type = "II")

# Elevation plot
comps <- glht(m2, linfct = mcp(Elevation = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Elevation = "Tukey")), Letters = letters)
#plot(tukey_letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Elevation=x,Pred = y)
#df$fructicose <- datos$fructicose

g1 <- ggplot(df, aes(x = Elevation, y = Pred, group=Elevation)) +
  geom_boxplot(outliers = F) +
  theme_minimal() +
  labs(x = "Elevation", y = "Fructicose", title = "A") +
  theme_bw()

tukey_df <- data.frame(
  Elevation = names(tukey_letters$mcletters$Letters),  # The Elevation levels ("0", "1", "2")
  Letters = tukey_letters$mcletters$Letters)            # The Tukey letters ("a", "a", "a")

g2 <- g1 + geom_text(data = tukey_df, 
                                    aes(x = Elevation, y = max(tukey_letters$lp) + 
                                          0.01, label = Letters), vjust = -0.5, size = 5, color = "black")

g.fructicose.elev <- g2 + geom_jitter(data= datos, aes(x=Elevation, y = fructicose, color = Slope), 
                                      width = 0.2, alpha=0.9) + plot.parameters

g.fructicose.elev

### Slope plot ###

comps <- glht(m2, linfct = mcp(Slope = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Slope = "Tukey")), Letters = letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Slope=x,Pred = y)

g1 <- ggplot(df, aes(x = Slope, y = Pred, group=Slope)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Slope", y = "Predicted richness of fructicose species", title = "B") +
  theme_bw()

tukey_df <- data.frame(
  Slope = names(tukey_letters$mcletters$Letters),  
  Letters = tukey_letters$mcletters$Letters)   

g2 <- g1 + geom_text(data = tukey_df, aes(x = Slope, y = max(tukey_letters$lp) + 0.01, label = Letters), 
                                  vjust = -0.5, size = 5, color = "black")

g.fructicose.slope <- g2 + geom_jitter(data=datos,aes(x=Slope, y = fructicose, color = Elevation), 
                                       width = 0.2, alpha=0.9) + plot.parameters

g.fructicose.slope

### Aspect plot ###

or.prediction <- ggpredict(m2, terms = "Orientation")
g1 <- plot(or.prediction) + theme_bw() + labs(x="Aspect", y="Fructicose", title = "C")

or.plot <- g1 + geom_jitter(data= datos, aes(x=Orientation, y = fructicose, color = Elevation), 
                            width = 0.2, alpha=0.9) + plot.parameters 
or.plot

gridExtra::grid.arrange(g.fructicose.elev, g.fructicose.slope, or.plot)

#####################################################
############Crustoso##############################
#####################################################
#####################################################

#m2 <- glm(crustose ~ Elevation + Slope + Orientation, data = datos,family = gaussian(link = "identity"))
m2 <- glm(crustose ~ Elevation + Slope + Orientation, data = datos,family = quasipoisson(link = "log"))
summary(m2)
Anova(m2)

### Elevation plot ###
comps <- glht(m2, linfct = mcp(Elevation = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Elevation = "Tukey")), Letters = letters)
plot(tukey_letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Elevation=x,Pred = y)
df$crustose <- datos$crustose

g1 <- ggplot(df, aes(x = Elevation, y = Pred, group=Elevation)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Elevation", y = "Crustose", title = "A") +
  theme_bw()

tukey_df <- data.frame(
  Elevation = names(tukey_letters$mcletters$Letters),  # The Elevation levels ("0", "1", "2")
  Letters = tukey_letters$mcletters$Letters)            # The Tukey letters ("a", "a", "a")

g2 <- g1 + geom_text(data = tukey_df, 
                                    aes(x = Elevation, y = max(tukey_letters$lp) + 
                                          0.01, label = Letters), vjust = -0.5, size = 5, color = "black")

g.crustose.elev <- g2 + geom_jitter(data= datos, aes(x=Elevation, y = crustose, color = Slope), 
              width = 0.2, alpha=0.9) + plot.parameters

g.crustose.elev

### Slope plot ###

comps <- glht(m2, linfct = mcp(Slope = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Slope = "Tukey")), Letters = letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Slope=x,Pred = y)
df$crustose <- datos$crustose

g1 <- ggplot(df, aes(x = Slope, y = Pred, group=Slope)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Slope", y = "Predicted richness of crustose species", title = "B") +
  theme_bw()

tukey_df <- data.frame(
  Slope = names(tukey_letters$mcletters$Letters),  
  Letters = tukey_letters$mcletters$Letters)   

g2 <- g1 + geom_text(data = tukey_df, aes(x = Slope, y = max(tukey_letters$lp) + 0.01, label = Letters), 
                                     vjust = -0.5, size = 5, color = "black")

g.crustose.slope <- g2 + geom_jitter(data=datos,aes(x=Slope, y = crustose, color = Elevation), 
                                     width = 0.2, alpha=0.9) + plot.parameters

g.crustose.slope

### Aspect plot ###

or.prediction <- ggpredict(m2, terms = "Orientation")
g1 <- plot(or.prediction) + theme_bw() + labs(x="Aspect", y="Crustose", title = "C")

or.plot <- g1 + geom_jitter(data=datos,aes(x=Orientation, y = crustose, color = Elevation), 
                             width = 0.2, alpha=0.9) + plot.parameters
or.plot

gridExtra::grid.arrange(g.crustose.elev, g.crustose.slope, or.plot)

#####################################################
############Apotecios##############################
#####################################################
#####################################################

#m2 <- glm(apothecia ~ Elevation + Slope + Orientation, 
 #         data = datos,family = gaussian(link = "identity"))
#m2 <- glm(apothecia ~ Elevation + Slope + Orientation, 
  #        data = datos,family = quasipoisson(link = "log"))
m2 <- glm(apothecia ~ Elevation + Slope + Orientation, 
          data = datos,family = poisson(link = "log"))
Anova(m2)

### Elevation plot ###
comps <- glht(m2, linfct = mcp(Elevation = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Elevation = "Tukey")), Letters = letters)
plot(tukey_letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Elevation=x,Pred = y)
df$apotecios <- datos$apothecia

g1 <- ggplot(df, aes(x = Elevation, y = Pred, group=Elevation)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Elevation", y = "Apothecia", title = "A") +
  theme_bw()

tukey_df <- data.frame(
  Elevation = names(tukey_letters$mcletters$Letters),  # The Elevation levels ("0", "1", "2")
  Letters = tukey_letters$mcletters$Letters)            # The Tukey letters ("a", "a", "a")

g2 <- g1 + geom_text(data = tukey_df, 
                                    aes(x = Elevation, y = max(tukey_letters$lp) + 
                                          0.01, label = Letters), vjust = -0.5, size = 5, color = "black")

g.apothecia.elev <- g2 + geom_jitter(data= datos, aes(x=Elevation, y = apothecia, color = Slope), 
                                     width = 0.2, alpha=0.9) + plot.parameters

g.apothecia.elev

### Slope plot ###

comps <- glht(m2, linfct = mcp(Slope = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Slope = "Tukey")), Letters = letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Slope=x,Pred = y)

g1 <- ggplot(df, aes(x = Slope, y = Pred, group=Slope)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Slope", y = "Predicted frequency of apotheciate lichens", title = "B") +
  theme_bw()

tukey_df <- data.frame(
  Slope = names(tukey_letters$mcletters$Letters),  
  Letters = tukey_letters$mcletters$Letters)   

g2 <- g1 + geom_text(data = tukey_df, aes(x = Slope, y = max(tukey_letters$lp) + 0.01, label = Letters), 
                                     vjust = -0.5, size = 5, color = "black")

g.apothecia.slope <- g2 + geom_jitter(data=datos,aes(x=Slope, y = apothecia, color = Elevation), 
                                      width = 0.2, alpha=0.9) + plot.parameters

g.apothecia.slope

### Aspect plot ###

or.prediction <- ggpredict(m2, terms = "Orientation")
g1 <- plot(or.prediction) + theme_bw() + labs(x="Aspect", y="Apothecia", title = "C")

or.plot <- g1 + geom_jitter(data=datos,aes(x=Orientation, y = apothecia, color = Elevation), 
                             width = 0.2, alpha=0.9) + plot.parameters
or.plot

gridExtra::grid.arrange(g.apothecia.elev, g.apothecia.slope, or.plot)

#####################################################
############   Soredios   ###########################
#####################################################
#####################################################

#m2 <- glm(soredia ~ Elevation + Slope + Orientation, data = datos,family = gaussian(link = "identity"))
m2 <- glm(soredia ~ Elevation + Slope + Orientation, data = datos,family = quasipoisson(link = "log"))

summary(m2)
Anova(m2,type = "II")

### Elevation plot ###
comps <- glht(m2, linfct = mcp(Elevation = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Elevation = "Tukey")), Letters = letters)
plot(tukey_letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Elevation=x,Pred = y)
df$soredios <- datos$soredia

g1 <- ggplot(df, aes(x = Elevation, y = Pred, group=Elevation)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Elevation", y = "Soredia", title = "A") +
  theme_bw()

tukey_df <- data.frame(
  Elevation = names(tukey_letters$mcletters$Letters),  # The Elevation levels ("0", "1", "2")
  Letters = tukey_letters$mcletters$Letters)            # The Tukey letters ("a", "a", "a")

g2 <- g1 + geom_text(data = tukey_df, 
                                    aes(x = Elevation, y = max(tukey_letters$lp) + 
                                          0.01, label = Letters), vjust = -0.5, size = 5, color = "black")

g.soredia.elev <- g2 + geom_jitter(data= datos, aes(x=Elevation, y = soredia, color = Slope), 
                                   width = 0.2, alpha=0.9) + plot.parameters

g.soredia.elev

### Slope plot ###

comps <- glht(m2, linfct = mcp(Slope = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Slope = "Tukey")), Letters = letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Slope=x,Pred = y)

g1 <- ggplot(df, aes(x = Slope, y = Pred, group=Slope)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Slope", y = "Predicted richness of sorediate lichens", title = "B") +
  theme_bw()

tukey_df <- data.frame(
  Slope = names(tukey_letters$mcletters$Letters),  
  Letters = tukey_letters$mcletters$Letters)   

g2 <- g1 + geom_text(data = tukey_df, aes(x = Slope, y = max(tukey_letters$lp) + 0.01, label = Letters), 
                     vjust = -0.5, size = 5, color = "black")

g.soredia.slope <- g2 + geom_jitter(data=datos,aes(x=Slope, y = soredia, color = Elevation), 
                                    width = 0.2, alpha=0.9) + plot.parameters

g.soredia.slope

### Aspect plot ###

or.prediction <- ggpredict(m2, terms = "Orientation")

or.plot <- plot(or.prediction) + theme_bw() + labs(x="Aspect", y="Soredia", title = "C")+
  geom_jitter(data=datos,aes(x=Orientation, y = soredia, color = Elevation), 
               width = 0.2, alpha=0.9) + plot.parameters

gridExtra::grid.arrange(g.soredia.elev, g.soredia.slope, or.plot)

#####################################################
############Isidios##############################
#####################################################
#####################################################

#m2 <- glm(isidia ~ Elevation + Slope + Orientation, data = datos,family = gaussian(link = "identity"))
m2 <- glm(isidia ~ Elevation + Slope + Orientation, data = datos,family = quasipoisson(link = "log"))

summary(m2)
Anova(m2,type = "II")

### Elevation plot ###
comps <- glht(m2, linfct = mcp(Elevation = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Elevation = "Tukey")), Letters = letters)
plot(tukey_letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Elevation=x,Pred = y)

g1 <- ggplot(df, aes(x = Elevation, y = Pred, group=Elevation)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Elevation", y = "Isidia", title = "A") +
  theme_bw()

tukey_df <- data.frame(
  Elevation = names(tukey_letters$mcletters$Letters),  # The Elevation levels ("0", "1", "2")
  Letters = tukey_letters$mcletters$Letters)            # The Tukey letters ("a", "a", "a")

g.isidia.elev <- g1 + geom_text(data = tukey_df, 
                                    aes(x = Elevation, y = max(tukey_letters$lp) + 
                                          0.01, label = Letters), vjust = -0.5, size = 5, color = "black")+
  geom_jitter(data= datos, aes(x=Elevation, y = isidia, color = Slope), 
              width = 0.2, alpha=0.9) + plot.parameters

### Slope plot ###

comps <- glht(m2, linfct = mcp(Slope = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Slope = "Tukey")), Letters = letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Slope=x,Pred = y)

g1 <- ggplot(df, aes(x = Slope, y = Pred, group=Slope)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Slope", y = "Predicted Richness of isidiate species", title = "B") +
  theme_bw()

tukey_df <- data.frame(
  Slope = names(tukey_letters$mcletters$Letters),  
  Letters = tukey_letters$mcletters$Letters)   

g.isidia.slope <- g1 + geom_text(data = tukey_df, aes(x = Slope, y = max(tukey_letters$lp) + 0.01, label = Letters), 
                                     vjust = -0.5, size = 5, color = "black")+
  geom_jitter(data=datos,aes(x=Slope, y = isidia, color = Elevation), 
              width = 0.2, alpha=0.9) + plot.parameters

gridExtra::grid.arrange(g.isidia.elev, g.isidia.slope)

#####################################################
###########Cilia##############################
#####################################################
#####################################################

#m2 <- glm(cilia ~ Elevation + Slope + Orientation, data = datos,family = gaussian(link = "identity"))
m2 <- glm(cilia ~ Elevation + Slope + Orientation, data = datos,family = poisson(link = "log"))
summary(m2)
Anova(m2)

### Elevation plot ###
comps <- glht(m2, linfct = mcp(Elevation = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Elevation = "Tukey")), Letters = letters)
plot(tukey_letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Elevation=x,Pred = y)

g1 <- ggplot(df, aes(x = Elevation, y = Pred, group=Elevation)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Elevation", y = "Cilia", title = "A") +
  theme_bw()

tukey_df <- data.frame(
  Elevation = names(tukey_letters$mcletters$Letters),  # The Elevation levels ("0", "1", "2")
  Letters = tukey_letters$mcletters$Letters)            # The Tukey letters ("a", "a", "a")

g.cilia.elev <- g1 + geom_text(data = tukey_df, 
                                    aes(x = Elevation, y = max(tukey_letters$lp) + 
                                          0.01, label = Letters), vjust = -0.5, size = 5, color = "black")+
  geom_jitter(data= datos, aes(x=Elevation, y = cilia, color = Slope), 
              width = 0.2, alpha=0.9) + plot.parameters


### Slope plot ###

comps <- glht(m2, linfct = mcp(Slope = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Slope = "Tukey")), Letters = letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Slope=x,Pred = y)

g1 <- ggplot(df, aes(x = Slope, y = Pred, group=Slope)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Slope", y = "Predicted richness of ciliate species", title = "B") +
  theme_bw()

tukey_df <- data.frame(
  Slope = names(tukey_letters$mcletters$Letters),  
  Letters = tukey_letters$mcletters$Letters)   

g.cilia.slope <- g1 + geom_text(data = tukey_df, aes(x = Slope, y = max(tukey_letters$lp) + 0.01, label = Letters), 
                                     vjust = -0.5, size = 5, color = "black")+
  geom_jitter(data=datos,aes(x=Slope, y = cilia, color = Elevation), 
              width = 0.2, alpha=0.9) + plot.parameters


### Aspect plot ###

or.prediction <- ggpredict(m2, terms = "Orientation")
or.plot <- plot(or.prediction) + theme_bw() + labs(x="Aspect", y="Cilia", title = "C")+
  geom_jitter(data=datos,aes(x=Orientation, y = cilia, color = Elevation), 
               width = 0.2, alpha=0.9) + plot.parameters

#windows()
gridExtra::grid.arrange(g.cilia.elev, g.cilia.slope, or.plot)

#####################################################
###########UV protection##############################
#####################################################
#####################################################

#m2 <- glm(uv_protection ~ Elevation + Slope + Orientation, data = datos,family = gaussian(link = "identity"))
m2 <- glm(uv_protection ~ Elevation + Slope + Orientation, data = datos,family = poisson(link = "identity"))

summary(m2)
Anova(m2,type = "II")

### Slope plot ###
comps <- glht(m2, linfct = mcp(Slope = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Slope = "Tukey")), Letters = letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Slope=x,Pred = y)

g1 <- ggplot(df, aes(x = Slope, y = Pred, group=Slope)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Slope", y = "UV protection", title = "B") +
  theme_bw()

tukey_df <- data.frame(
  Slope = names(tukey_letters$mcletters$Letters),  
  Letters = tukey_letters$mcletters$Letters)   

g.uvprotection.slope <- g1 + geom_text(data = tukey_df, aes(x = Slope, y = max(tukey_letters$lp) + 0.01, label = Letters), 
                                     vjust = -0.5, size = 5, color = "black")+
  geom_jitter(data=datos,aes(x=Slope, y = uv_protection, color = Elevation), 
              width = 0.2, alpha=0.9) + plot.parameters

g.uvprotection.slope

#####################################################
########### Microbiological_protec ##################
#####################################################
#####################################################

#m2 <- glm(Microbiological_protec ~ Elevation + Slope + Orientation, data = datos,family = gaussian(link = "identity"))
m2 <- glm(Microbiological_protec ~ Elevation + Slope + Orientation, data = datos,family = quasipoisson(link = "log"))

summary(m2)
Anova(m2,type = "II")

### Elevation plot ###
comps <- glht(m2, linfct = mcp(Elevation = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Elevation = "Tukey")), Letters = letters)
plot(tukey_letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Elevation=x,Pred = y)

g1 <- ggplot(df, aes(x = Elevation, y = Pred, group=Elevation)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Elevation", y = "Microbiological", title = "A") +
  theme_bw()

tukey_df <- data.frame(
  Elevation = names(tukey_letters$mcletters$Letters),  # The Elevation levels ("0", "1", "2")
  Letters = tukey_letters$mcletters$Letters)            # The Tukey letters ("a", "a", "a")

g.micro.elev <- g1 + geom_text(data = tukey_df, 
                                    aes(x = Elevation, y = max(tukey_letters$lp) + 
                                          0.01, label = Letters), vjust = -0.5, size = 5, color = "black")+
  geom_jitter(data= datos, aes(x=Elevation, y = Microbiological_protec, color = Slope), 
              width = 0.2, alpha=0.9) + plot.parameters

### Slope plot ###

comps <- glht(m2, linfct = mcp(Slope = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Slope = "Tukey")), Letters = letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Slope=x,Pred = y)

g1 <- ggplot(df, aes(x = Slope, y = Pred, group=Slope)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Slope", y = "Microbiological protection", title = "B") +
  theme_bw()

tukey_df <- data.frame(
  Slope = names(tukey_letters$mcletters$Letters),  
  Letters = tukey_letters$mcletters$Letters)   

g.micro.slope <- g1 + geom_text(data = tukey_df, aes(x = Slope, y = max(tukey_letters$lp) + 0.01, label = Letters), 
                                     vjust = -0.5, size = 5, color = "black")+
  geom_jitter(data=datos,aes(x=Slope, y = Microbiological_protec, color = Elevation), 
              width = 0.2, alpha=0.9) + plot.parameters

g.micro.slope

### Aspect plot ###

or.prediction <- ggpredict(m2, terms = "Orientation")
or.micro.plot <- plot(or.prediction) + theme_bw() + labs(x="Aspect", y="microbiological", title = "C")+
  geom_jitter(data=datos,aes(x=Orientation, y = Microbiological_protec, color = Elevation), 
               width = 0.2, alpha=0.9) + plot.parameters

gridExtra::grid.arrange(g.micro.elev, g.micro.slope, or.micro.plot)

#####################################################
###########herbivory_protec##############################
#####################################################
#####################################################

m2 <- glm(herbivory_protec ~ Elevation + Slope + Orientation, data = datos,family = poisson(link = "log"))

#summary(m2)
Anova(m2,type = "II")

### Slope plot ###

comps <- glht(m2, linfct = mcp(Slope = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Slope = "Tukey")), Letters = letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Slope=x,Pred = y)

g1 <- ggplot(df, aes(x = Slope, y = Pred, group=Slope)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Slope", y = "Herbivory protection", title = "B") +
  theme_bw()

tukey_df <- data.frame(
  Slope = names(tukey_letters$mcletters$Letters),  
  Letters = tukey_letters$mcletters$Letters)   

g.herb.slope <- g1 + geom_text(data = tukey_df, aes(x = Slope, y = max(tukey_letters$lp) + 0.01, label = Letters), 
                                     vjust = -0.5, size = 5, color = "black")+
  geom_jitter(data=datos,aes(x=Slope, y = herbivory_protec, color = Elevation), 
              width = 0.2, alpha=0.9) + plot.parameters

g.herb.slope

### Aspect plot ###

or.prediction <- ggpredict(m2, terms = "Orientation")
or.herb <- plot(or.prediction) + theme_bw() + labs(x="Aspect", y="Herbivory protection", title = "C")+
  geom_jitter(data=datos,aes(x=Orientation, y = herbivory_protec, color = Elevation), 
              width = 0.2, alpha=0.9) + plot.parameters


#gridExtra::grid.arrange(g.herb.slope, or.herb)

windows()
gridExtra::grid.arrange(g.micro.elev, g.micro.slope, or.micro.plot,
                        g.herb.slope, or.herb,g.uvprotection.slope, ncol=3)


###########################################################################################################3
############################################################################################################
##################### N M S #############################################################################
##########################################################################################################
setwd()

library("vegan")
library("dplyr")
library("ggplot2")

bio <- read.table("main_reducida.csv",
                  header = T,
                  sep = ";",
                  dec = ".",
                  row.names = 1) 
head(bio)
env <- read.table("secundaria_reducida.csv",
                  header = T,
                  sep = ";",
                  dec = ".",
                  row.names = 1)
head(env)

nms <- metaMDS(comm = bio,
               autotransform = FALSE,
               distance = "bray",
               engine = "monoMDS",
               k = 3,
               weakties = TRUE,
               model = "global",
               maxit = 300,
               try = 40,
               trymax = 100)

envnms <- envfit(nms, env)

windows()
plot(nms)
plot(envnms, p.max = 0.05, col = "red", cex = 1)

species_scores <- scores(nms, display = "species")

site_data <- as.data.frame( scores(nms, display="sites") )
site_data <- cbind(rownames(site_data), site_data)
rownames(site_data) <- NULL
colnames(site_data) <- c("ID","NMDS1","NMDS2", "NMDS3")
site_data$ID

sitesID <- as.data.frame(env)
sitesID <- cbind(rownames(sitesID), sitesID)
rownames(sitesID) <- NULL
colnames(sitesID) <- c("ID","or_R", "tipo_roc", "or_L", "exp", "pend", "piso")
sitesID

combined_data <- inner_join(site_data, sitesID, by="ID")

df_environ  <- scores(envnms, display = 'bp')
df_environ <- as.data.frame(df_environ)

num_elevations <- length(unique(combined_data$piso))
shape_values <- 1:num_elevations
shape_values <- c(19,17,19,17,15,17,19)
elevation_colors <- c("black", "darkblue", "darkgreen", "darkred", "brown1", "darkgoldenrod", "deepskyblue3")

names <- as.data.frame(species_scores) 

names <- cbind(rownames(names), names)
rownames(names) <- NULL
colnames(names) <- c("species","NMDS1","NMDS2", "NMDS3")
names$species
names.names <- names$species

species_scores <- as.data.frame(species_scores)

windows()

g1 <- ggplot(combined_data, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(shape = as.factor(piso), color = as.factor(piso)),size=2) +  # Set shape within aes() and specify size
  geom_text(data = names, aes(x = NMDS1, y = NMDS2, label = species)) +
  labs(x = "NMDS1", y = "NMDS2", title = "") +
  scale_shape_manual(values = shape_values) + # Set shape values
  scale_color_manual(values = elevation_colors, name = "Elevation") +
  guides(shape = guide_legend(title = "Elevation"))+  # Add legend for elevation 
  theme_bw()
g1

g2 <- g1 +  geom_segment(data=df_environ, 
                         aes(x=0, #Starting coordinate in CCA1 = 0 
                             y=0, #Start in CCA2 = 0
                             xend=NMDS1,#Ending coordinate in CCA1  
                             yend=NMDS2), #Ending coordinate in CCA2 
                         color="darkred", #set color
                         linewidth=1,
                         arrow=arrow(length=unit(0.01,"npc")))#Set the size of the lines that form the tip of the arrow

g3 <- g2 +  geom_text(data=df_environ,
                      aes(x=NMDS1, 
                          y=NMDS2,
                          label=rownames(df_environ),
                          hjust=0.2*(1-sign(NMDS1)),#Add the text of each environmental var at the end of the arrow
                          vjust=0.5*(1-sign(NMDS2))),#Add the text of each environmental var at the end of the arrow 
                      color="darkred")  #Set x and y axis titles
windows()
g3


######################################################
################ NMS PLOT ############################
######################################################

rm(list = ls())
nms <- read.table("nms.csv", header = T, sep = ";", dec = ".")
head(nms)

###### CLAUDE

# Packages
library(ggplot2)
library(vegan)
library(ggrepel)

# Make sure grouping variables are factors
nms$altura3 <- factor(nms$altura3,
                      levels = c(1, 2, 3),
                      labels = c("1", "2", "3"))

nms$Pend2 <- factor(nms$Pend2,
                    levels = c(1, 2),
                    labels = c("1", "2"))

# Environmental variables to use as arrows
env_vars <- c(
  "FoliosoCWM", "FructiculosoCWM", "CrustosoCWM", "EscuamulosoCWM",
  "DimorficoCWM", "ApoteciosCWM", "SorediosCWM", "IsidiosCWM",
  "P.HerbCWM", "CiliadosCWM", "P.UVCWM", "Altura", "Pendiente"
)

# NMS axes
ord_axes <- nms[, c("nms1", "nms2")]

# Environmental matrix
env_data <- nms[, env_vars]

# Fit environmental vectors
env_fit <- envfit(ord_axes, env_data, permutations = 999)

# Extract arrow coordinates
arrows <- as.data.frame(scores(env_fit, display = "vectors"))
arrows$variable <- rownames(arrows)

# Rename columns for clarity
colnames(arrows)[1:2] <- c("nms1", "nms2")

# Scale arrows so they fit nicely in the plot
arrow_multiplier <- 1.2 * min(
  diff(range(nms$nms1)) / diff(range(arrows$nms1)),
  diff(range(nms$nms2)) / diff(range(arrows$nms2))
)

arrows$nms1 <- arrows$nms1 * arrow_multiplier
arrows$nms2 <- arrows$nms2 * arrow_multiplier

# Biplot
windows()
ggplot(nms, aes(x = nms1, y = nms2)) +
  
  # Points
  geom_point(
    aes(color = altura3, shape = Pend2),
    size = 3
  ) +
  
  # Environmental arrows
  geom_segment(
    data = arrows,
    aes(x = 0, y = 0, xend = nms1, yend = nms2),
    arrow = arrow(length = unit(0.25, "cm")),
    color = "black",
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  
  # Arrow labels
  geom_text_repel(
    data = arrows,
    aes(x = nms1, y = nms2, label = variable),
    size = 3.5,
    color = "black",
    inherit.aes = FALSE
  ) +
  
  # Manual colors for altura3
  scale_color_manual(
    values = c(
      "1" = "violet",
      "2" = "lightgrey",
      "3" = "lightblue"
    ),
    name = "Altura3"
  ) +
  
  # Manual shapes for Pend2
  scale_shape_manual(
    values = c(
      "1" = 16,  # circle
      "2" = 17   # triangle
    ),
    name = "Pend2"
  ) +
  
  labs(
    x = "NMS1",
    y = "NMS2",
    title = "NMS biplot with environmental vectors"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

################### BIPLOTS WITH BOXPLOTS OP1 ##############

library(ggplot2)
library(ggrepel)
library(vegan)
library(dplyr)
library(patchwork)

windows()

# Make sure factors are correctly defined
nms$altura3 <- factor(nms$altura3,
                      levels = c(1, 2, 3),
                      labels = c("1", "2", "3"))

nms$Pend2 <- factor(nms$Pend2,
                    levels = c(1, 2),
                    labels = c("1", "2"))

# Environmental variables
env_vars <- c(
  "FoliosoCWM", "FructiculosoCWM", "CrustosoCWM", "EscuamulosoCWM",
  "DimorficoCWM", "ApoteciosCWM", "SorediosCWM", "IsidiosCWM",
  "P.HerbCWM", "CiliadosCWM", "P.UVCWM", "Altura", "Pendiente"
)

ord_axes <- nms[, c("nms1", "nms2")]
env_data <- nms[, env_vars]

env_fit <- envfit(ord_axes, env_data, permutations = 999)

arrows <- as.data.frame(scores(env_fit, display = "vectors"))
arrows$variable <- rownames(arrows)
colnames(arrows)[1:2] <- c("nms1", "nms2")

# Scale arrows
arrow_multiplier <- 1.2 * min(
  diff(range(nms$nms1)) / diff(range(arrows$nms1)),
  diff(range(nms$nms2)) / diff(range(arrows$nms2))
)

arrows$nms1 <- arrows$nms1 * arrow_multiplier
arrows$nms2 <- arrows$nms2 * arrow_multiplier

# Main NMS biplot
p_nms <- ggplot(nms, aes(x = nms1, y = nms2)) +
  
  geom_point(
    aes(color = altura3, shape = Pend2),
    size = 3
  ) +
  
  geom_segment(
    data = arrows,
    aes(x = 0, y = 0, xend = nms1, yend = nms2),
    arrow = arrow(length = unit(0.25, "cm")),
    color = "black",
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  
  geom_text_repel(
    data = arrows,
    aes(x = nms1, y = nms2, label = variable),
    size = 3.5,
    color = "black",
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "1" = "darkviolet",
      "2" = "darkgreen",
      "3" = "blue"
    ),
    name = "Altura3"
  ) +
  
  scale_shape_manual(
    values = c(
      "1" = 16,
      "2" = 17
    ),
    name = "Pend2"
  ) +
  
  labs(
    x = "NMS1",
    y = "NMS2",
    title = "NMS biplot with environmental vectors"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

# Create grouping variable for boxplots
nms <- nms %>%
  mutate(group_box = interaction(altura3, Pend2, sep = " / Pend2 "))

# Horizontal boxplot of NMS1 scores
p_box_nms1 <- ggplot(nms, aes(x = nms1, y = group_box)) +
  
  geom_boxplot(
    aes(fill = altura3),
    coef = Inf,
    alpha = 0.7,
    color = "black"
  ) +
  
  geom_point(
    aes(shape = Pend2),
    size = 2,
    alpha = 0.7,
    position = position_jitter(height = 0.08, width = 0)
  ) +
  
  scale_fill_manual(
    values = c(
      "1" = "darkviolet",
      "2" = "darkgreen",
      "3" = "blue"
    ),
    name = "Altura3"
  ) +
  
  scale_shape_manual(
    values = c(
      "1" = 16,
      "2" = 17
    ),
    name = "Pend2"
  ) +
  
  labs(
    x = "NMS1 score",
    y = "Altura3 / Pend2",
    title = "Horizontal distribution of NMS1 scores by Altura3 and Pend2"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# Combine biplot and horizontal boxplot
p_nms / p_box_nms1 +
  plot_layout(heights = c(3, 1.2))

############ OPTION 2 ##########

library(tidyr)

nms_long <- nms %>%
  pivot_longer(
    cols = c(nms1, nms2),
    names_to = "axis",
    values_to = "score"
  ) %>%
  mutate(group_box = interaction(altura3, Pend2, sep = " / Pend2 "))

ggplot(nms_long, aes(x = score, y = group_box)) +
  
  geom_boxplot(
    aes(fill = altura3),
    coef = Inf,
    alpha = 0.7,
    color = "black"
  ) +
  
  geom_point(
    aes(shape = Pend2),
    size = 2,
    alpha = 0.7,
    position = position_jitter(height = 0.08, width = 0)
  ) +
  
  facet_wrap(~ axis, scales = "free_x") +
  
  scale_fill_manual(
    values = c(
      "1" = "violet",
      "2" = "lightgrey",
      "3" = "lightblue"
    ),
    name = "Altura3"
  ) +
  
  scale_shape_manual(
    values = c(
      "1" = 16,
      "2" = 17
    ),
    name = "Pend2"
  ) +
  
  labs(
    x = "NMS score",
    y = "Altura3 / Pend2",
    title = "Distribution of NMS scores by elevation and slope group"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

############ OPTION 3 ###############

library(dplyr)

# Summary of min, median, max and quartiles for NMS1
box_summary <- nms %>%
  group_by(altura3, Pend2) %>%
  summarise(
    xmin = min(nms1, na.rm = TRUE),
    q1 = quantile(nms1, 0.25, na.rm = TRUE),
    median = median(nms1, na.rm = TRUE),
    q3 = quantile(nms1, 0.75, na.rm = TRUE),
    xmax = max(nms1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(altura3, Pend2)

# Define y positions below the ordination cloud
yrange <- range(c(nms$nms2, arrows$nms2), na.rm = TRUE)
xrange <- range(c(nms$nms1, arrows$nms1), na.rm = TRUE)

y_step <- diff(yrange) * 0.08
box_summary$ypos <- yrange[1] - y_step * seq_len(nrow(box_summary))

box_summary$ymin <- box_summary$ypos - y_step * 0.20
box_summary$ymax <- box_summary$ypos + y_step * 0.20

box_summary$label <- paste0("A", box_summary$altura3, " / P", box_summary$Pend2)

# Biplot with embedded horizontal range summaries
ggplot(nms, aes(x = nms1, y = nms2)) +
  
  geom_point(
    aes(color = altura3, shape = Pend2),
    size = 3
  ) +
  
  geom_segment(
    data = arrows,
    aes(x = 0, y = 0, xend = nms1, yend = nms2),
    arrow = arrow(length = unit(0.25, "cm")),
    color = "black",
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  
  geom_text_repel(
    data = arrows,
    aes(x = nms1, y = nms2, label = variable),
    size = 3.5,
    color = "black",
    inherit.aes = FALSE
  ) +
  
  # Min-max line
  geom_segment(
    data = box_summary,
    aes(x = xmin, xend = xmax, y = ypos, yend = ypos, color = altura3),
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  
  # Interquartile rectangle
  geom_rect(
    data = box_summary,
    aes(
      xmin = q1,
      xmax = q3,
      ymin = ymin,
      ymax = ymax,
      fill = altura3
    ),
    color = "black",
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  
  # Median line
  geom_segment(
    data = box_summary,
    aes(
      x = median,
      xend = median,
      y = ymin,
      yend = ymax
    ),
    color = "black",
    linewidth = 0.9,
    inherit.aes = FALSE
  ) +
  
  # Group labels
  geom_text(
    data = box_summary,
    aes(
      x = xrange[1],
      y = ypos,
      label = label
    ),
    hjust = 1.1,
    size = 3,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "1" = "violet",
      "2" = "lightgrey",
      "3" = "lightblue"
    ),
    name = "Altura3"
  ) +
  
  scale_fill_manual(
    values = c(
      "1" = "violet",
      "2" = "lightgrey",
      "3" = "lightblue"
    ),
    name = "Altura3"
  ) +
  
  scale_shape_manual(
    values = c(
      "1" = 16,
      "2" = 17
    ),
    name = "Pend2"
  ) +
  
  coord_cartesian(
    ylim = c(min(box_summary$ymin) - y_step, yrange[2]),
    clip = "off"
  ) +
  
  labs(
    x = "NMS1",
    y = "NMS2",
    title = "NMS biplot with embedded horizontal group summaries"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 60)
  )

########## OTHER OPTIONS FOR A ACCURATE FILTER EFFECT VISUALIZATION #########

library(vegan)

# Use NMS coordinates
ord_dist <- dist(nms[, c("nms1", "nms2")])

# Test dispersion by elevation class
bd_altura <- betadisper(ord_dist, nms$altura3)

permutest(bd_altura, permutations = 999)

TukeyHSD(bd_altura)







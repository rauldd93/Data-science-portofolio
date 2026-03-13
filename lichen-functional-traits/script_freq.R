rm(list=ls())

setwd("C:/R/ramiro")

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
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Elevation = "Tukey")), Letters = letters)
y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Elevation=x,Pred = y)

g1 <- ggplot(df, aes(x = Elevation, y = Pred, group=Elevation)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Elevation", y = "Predicted frequency of foliose species", title = "A") +
  theme_bw()

tukey_df <- data.frame(
  Elevation = names(tukey_letters$mcletters$Letters),  # The Elevation levels ("0", "1", "2")
  Letters = tukey_letters$mcletters$Letters)            # The Tukey letters ("a", "a", "a")

windows()
g.folioso.elev <- g1 + geom_text(data = tukey_df, 
                                 aes(x = Elevation, y = max(tukey_letters$lp) + 0.01, label = Letters), 
                                 vjust = -0.5, size = 5, color = "black")+
  geom_jitter(data= datos, aes(x=Elevation, y = foliose, color = Slope), width = 0.2, alpha=0.7)+
  #scale_color_viridis_d()
  scale_color_manual(values = c("Low" = "#E1BE6A", "High" = "#40B0A6"))

g.folioso.elev

#+theme(axis.text.x = element_text(angle = 20, hjust = 1,vjust=1))
#g.folioso.elev + 
 # scale_fill_viridis(discrete = TRUE, alpha=0.6, option = "A") +
  #geom_jitter(color="black", size=1.5, alpha=0.5)


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
  geom_jitter(data=datos,aes(x=Slope, y = foliose, color = Elevation), width = 0.2, alpha=0.6)

#windows()
gridExtra::grid.arrange(g.folioso.elev, g.folioso.slope)

#####################################################
#############FructiculosoCWM##############################
#####################################################
#####################################################

## Histograma variable respuesta ####
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
df$fructicose <- datos$fructicose

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
g2$data$datos <- df$fructicose

g.fructicose.elev <- g2 + geom_jitter(aes(x=Elevation, y = datos), width = 0.2)
g.fructicose.elev

### Slope plot ###

comps <- glht(m2, linfct = mcp(Slope = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Slope = "Tukey")), Letters = letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Slope=x,Pred = y)
df$fructicose <- datos$fructicose

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
g2$data$datos <- df$fructicose

g.fructicose.slope <- g2 + geom_jitter(aes(x=Slope, y = datos), width = 0.2)
g.fructicose.slope

### Aspect plot ###

or.prediction <- ggpredict(m2, terms = "Orientation")
g1 <- plot(or.prediction) + theme_bw() + labs(x="Aspect", y="Fructicose", title = "C")

or.plot <- g1 + geom_jitter(data=datos, aes(x=Orientation, y = fructicose))
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

g2$data$datos <- df$crustose
g.crustose.elev <- g2 + geom_jitter(aes(x=Elevation, y = datos), width = 0.2)
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

g2$data$datos <- df$crustose

g.crustose.slope <- g2 + geom_jitter(aes(x=Slope, y = datos), width = 0.2)
g.crustose.slope

### Aspect plot ###

or.prediction <- ggpredict(m2, terms = "Orientation")
g1 <- plot(or.prediction) + theme_bw() + labs(x="Slope", y="Crustose", title = "C")
or.plot <- g1 + geom_jitter(data=datos, aes(x=Orientation, y = crustose))
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
g2$data$datos <- df$apotecios

g.apothecia.elev <- g2 + geom_jitter(aes(x=Elevation, y = datos), width = 0.2)
g.apothecia.elev

### Slope plot ###

comps <- glht(m2, linfct = mcp(Slope = "Tukey"))
summary(comps)
tukey_letters <- cld(glht(m2, linfct = mcp(Slope = "Tukey")), Letters = letters)

y <- tukey_letters$lp 
x <- tukey_letters$x
df <- data.frame(Slope=x,Pred = y)
df$apotecios <- datos$apothecia

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

g2$data$datos <- df$apotecios

g.apothecia.slope <- g2 + geom_jitter(aes(x=Slope, y = datos), width = 0.2)
g.apothecia.slope

### Aspect plot ###

or.prediction <- ggpredict(m2, terms = "Orientation")
g1 <- plot(or.prediction) + theme_bw() + labs(x="Aspect", y="Apothecia", title = "C")
or.plot <- g1 + geom_jitter(data=datos, aes(x=Orientation, y = apothecia))
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

g.soredia.elev <- g2 + geom_jitter(data = datos, aes(x=Elevation, y = soredia), width = 0.2)
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

g.soredia.slope <- g1 + geom_text(data = tukey_df, aes(x = Slope, y = max(tukey_letters$lp) + 0.01, label = Letters), 
                                     vjust = -0.5, size = 5, color = "black") +
  geom_jitter(data=datos, aes(x=Slope, y=soredia), width = 0.2)

g.soredia.slope

### Aspect plot ###

or.prediction <- ggpredict(m2, terms = "Orientation")

or.plot <- plot(or.prediction) + theme_bw() + labs(x="Aspect", y="Soredia", title = "C")+
  geom_jitter(data=datos, aes(x=Orientation, y = soredia))

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
  geom_jitter(data=datos, aes(x=Elevation, y=isidia), width = 0.2)

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
  geom_jitter(data=datos, aes(x=Slope, y=isidia), width = 0.2)

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
  geom_jitter(data=datos, aes(x=Elevation, y=cilia), width = 0.2)

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
  geom_jitter(data=datos, aes(x=Slope, y=isidia), width = 0.2)

### Aspect plot ###

or.prediction <- ggpredict(m2, terms = "Orientation")
or.plot <- plot(or.prediction) + theme_bw() + labs(x="Aspect", y="Cilia", title = "C")+
  geom_jitter(data=datos, aes(x=Orientation, y = cilia))

windows()
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
  geom_jitter(data=datos, aes(x=Slope, y=uv_protection), width = 0.2)

g.uvprotection.slope

#####################################################
###########Microbiological_protec##############################
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
  geom_jitter(data=datos, aes(x=Elevation, y=Microbiological_protec), width = 0.2)

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
  geom_jitter(data=datos, aes(x=Slope, y=Microbiological_protec), width = 0.2)
g.micro.slope

### Aspect plot ###

or.prediction <- ggpredict(m2, terms = "Orientation")
or.micro.plot <- plot(or.prediction) + theme_bw() + labs(x="Aspect", y="microbiological", title = "C")+
  geom_jitter(data=datos, aes(x=Orientation, y = Microbiological_protec))

#gridExtra::grid.arrange(g.micro.elev, g.micro.slope, or.plot)

#####################################################
###########herbivory_protec##############################
#####################################################
#####################################################

## Histograma variable respuesta ####

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
  geom_jitter(data=datos, aes(x=Slope, y=herbivory_protection), width = 0.2)
g.herb.slope

### Aspect plot ###

or.prediction <- ggpredict(m2, terms = "Orientation")
or.herb <- plot(or.prediction) + theme_bw() + labs(x="Aspect", y="Herbivory protection", title = "C")+
  geom_jitter(data=datos, aes(x=Orientation, y = herbivory_protection))

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











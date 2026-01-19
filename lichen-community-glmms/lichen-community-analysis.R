# calling neccesary libreries
library("lattice")
library("smatr")
library("data.table")
library("emmeans")
library("multcomp")
library("postHoc")
library("ggplot2")
library("multcompView")
library("GGally")
library("car")
library("gridExtra")
library("lme4")
library("glmmTMB")
library("fitdistrplus")
library("ggeffects")
library("dplyr")

setwd("C:/R/...")

### Descriptive/explorative analysis ###
datos <- read.table("dat_comunidades.csv", sep = ";", dec = ".", header = T) # specify working directory

# inclination and explosition plot #
datos$pend_c <- factor(datos$pend_c, levels = c("Flat", "Low_inclination", "High_inclination"), ordered = F) # order the factor levels, this is important for the plots

windows() # optional
datos$exp[datos$exp < 0] <- datos$exp[datos$exp < 0] * -1 # transform negative to positive values
bwplot(exp ~ pend_c|piso, data = datos)
summary(lm(exp ~ pend_c+as.factor(piso), data = datos))
anova(lm(exp ~ pend_c+as.factor(piso), data = datos))

head(datos) # optional

datos_env <- datos %>%
  group_by(ud_muestral) %>%
  summarise(
    exp = mean(exp, na.rm = TRUE),  # average values
    pend_c = first(pend_c),         
    piso = first(piso),
    or_R = first(or_R)
  ) %>%
  ungroup()

m1 <- lm(exp ~ pend_c*as.factor(piso), data = datos_env)
anova(m1)
summary(m1)

ggplot(datos_env, aes(x = factor(pend_c), y = exp)) +  
  geom_boxplot() +
  facet_wrap(~ piso) +
  labs(x = "Rock inclination (°)", y = "Microsite exposition (°)", title = "") +
  theme_bw()

datos$piso = as.factor(datos$piso)

# Data cleaning...
datos <- droplevels(datos[!datos$morfotipo == "Liquen",])
datos <- droplevels(datos[!datos$morfotipo == "musgo",])
datos <- droplevels(datos[!datos$morfotipo == "roca",])
datos <- droplevels(datos[!datos$morfotipo == "Musgo",])
datos <- droplevels(datos[!datos$morfotipo == "Roca",])

### Richness calculation ###
library("data.table")
setDT(datos)
datos2 <- datos[  , .
                  (riqueza = length(spp_liquenes)), 
                  by = .(grad, ud_muestral,piso,
                         or_roca, ORNS_roca, or_R, exp, sup, pend, pend_c, tipo_roca,
                         or_ladera, ORNS_lad, or_L, paj, roc, arb, so,
                         mg,rast,ces,alt)]

datos3 <- datos[  , .
                  (riqueza = length(spp_liquenes)), 
                  by = .(grad, ud_muestral,piso,morfotipo,
                         or_roca, ORNS_roca, or_R, exp, sup, pend, pend_c, tipo_roca,
                         or_ladera, ORNS_lad, or_L, paj, roc, arb, so,
                         mg,rast,ces,alt)]

# Cleaning data...
datos3 <- droplevels(datos3[!datos3$morfotipo == "Liquen",])
datos3 <- droplevels(datos3[!datos3$morfotipo == "musgo",])
datos3 <- droplevels(datos3[!datos3$morfotipo == "roca",])
datos3 <- droplevels(datos3[!datos3$morfotipo == "Musgo",])
datos3 <- droplevels(datos3[!datos3$morfotipo == "Roca",])

# Richness of each grow form (crustose, foliose and fruticose)
DF4 <- datos3 %>%
mutate(cr = ifelse(morfotipo=="cr",datos3$riqueza,0),
       fl = ifelse(morfotipo=="fl",datos3$riqueza,0),
       fr = ifelse(morfotipo=="fr",datos3$riqueza,0))

windows()
ggplot( aes(x=piso, y=cr), data = DF4) +
  geom_violin(linetype="dashed") +
  geom_boxplot()+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + theme_bw()+
  ggtitle("") +
  xlab("")

ggplot( aes(x=pend_c, y=riqueza), data = datos2) +
  geom_violin(linetype="dashed") +
  geom_boxplot()+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + theme_bw()+
  ggtitle("") +
  xlab("")

### Richness models ###
subset(datos2, or_R == "Flat")

datos2$pend_c <- factor(datos2$pend_c, levels = c("Flat", "Low_inclination", "High_inclination"), 
                        ordered = F)

#fit3 <- glmer(riqueza ~ piso*or_R+piso*pend_c+(1|grad),
 #             family=poisson(link="log"), 
  #            data=datos2)
#Anova(fit3)
#summary(fit3)

fit3.1 <- glmer(riqueza ~ piso*or_R+pend_c+(1|grad),
              family=poisson(link="log"), 
              data=datos2)
anova(fit3, fit3.1, test="Chisq")
summary(fit3.1)
Anova(fit3.1)

vif(fit3.1)

fit3.2 <- glmer(riqueza ~ piso+or_R+pend_c+(1|grad),
                family=poisson(link="log"), 
                data=datos2)
anova(fit3.2, fit3.1, test="Chisq")
summary(fit3.2)
Anova(fit3.2)

#Validation
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}

overdisp_fun(fit3.2) # no overdispersion

fitted_values <- fitted(fit3.2)
pearson_resid <- residuals(fit3.2, type = "pearson")
residual_data <- data.frame(Fitted = fitted_values, Residuals = pearson_resid)
Vm1 <- ggplot(residual_data, aes(x = Fitted, y = Residuals)) +
  geom_point(color = "black") +  # Scatter plot of residuals
  geom_smooth(method = "lm", color = "blue", se = FALSE) +  # Regression line
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Horizontal line at zero
  labs(title = "A",
       x = "Fitted Values",
       y = "Pearson Residuals") +
  theme_bw()  # no pattern

#post hoc
rph <- emmeans(fit3.2, ~ piso)
cld_results <- cld(rph, Letters = letters)

rph2 <- emmeans(fit3.2, ~ pend_c)
cld_results <- cld(rph2, Letters = letters)

rph3 <- emmeans(fit3.2, ~ or_R)
cld_results <- cld(rph3, Letters = letters)

### Model plots ###
pred.1=ggpredict(fit3.2, terms = c("piso", "pend_c"), ci.lvl = 0.95)
pred.2=ggpredict(fit3.2, terms = c("piso", "or_R"), ci.lvl = 0.95)

pred.3=ggpredict(fit3.2, terms = c("pend_c"), ci.lvl = 0.95)
pred.4=ggpredict(fit3.2, terms = c("or_R"), ci.lvl = 0.95)

windows() # optional
A <- plot(pred.1, facets = F, connect.lines = F, colors = "bw",ci.style = "ribbon", dodge = 150) + theme_bw() + xlab("Elevation (m.a.s.l.)") + ylab("Richness")

B <- plot(pred.2, facets = F, connect.lines = F, colors = "bw", ci.style = "ribbon", dodge = 150) + theme_bw() + xlab("Elevation (m.a.s.l.)") + ylab("Richness")

C <- plot(pred.3, facets = F, connect.lines = F, colors = "bw", ci.style = "ribbon") + theme_bw() + xlab("Rock inclination") + ylab("Richness")

d <- plot(pred.4, facets = F, connect.lines = F, colors = "bw", ci.style = "ribbon") + theme_bw() + xlab("Rock aspect") + ylab("Richness")

windows() # optional
grid.arrange(A, B, ncol=1)

# GLM for crustose richness
fitCr <- glmer(cr ~ piso*or_R+piso*pend_c+(1|grad), family=poisson(link="log"), data=DF4)
summary(fitCr)
Anova(fitCr)

fitCr2 <- glmer(cr ~ piso*or_R+piso+pend_c+(1|grad), family=poisson(link="log"), data=DF4)
anova(fitCr2, fitCr, test="Chisq")
summary(fitCr2)
Anova(fitCr2)

fitCr3 <- glmer(cr ~ piso*or_R+(1|grad), family=poisson(link="log"), data=DF4) # final model
anova(fitCr3, fitCr2, test="Chisq")
summary(fitCr3)
Anova(fitCr3)

#Validation
overdisp_fun(fitCr3) # overdispersion

fitted_values <- fitted(fitCr3)
pearson_resid <- residuals(fitCr3, type = "pearson")
residual_data <- data.frame(Fitted = fitted_values, Residuals = pearson_resid)
Vm2 <- ggplot(residual_data, aes(x = Fitted, y = Residuals)) +
  geom_point(color = "black") +  # Scatter plot of residuals
  geom_smooth(method = "lm", color = "blue", se = FALSE) +  # Regression line
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Horizontal line at zero
  labs(title = "B",
       x = "Fitted Values",
       y = "") +
  theme_bw()  # no pattern

# Posthoc
ph_cr <- emmeans(fitCr3, ~ piso * or_R)
cld_results <- cld(ph_cr, Letters = letters)

subset(cld_results, or_R == 'Flat')

flat <- subset(cld_results, or_R == 'Flat')
north <- subset(cld_results, or_R == 'North')
south <- subset(cld_results, or_R == 'South')

# Plot for crustose GLM
pred.3 = ggpredict(fitCr3, terms = c("piso", "or_R"), ci.lvl = 0.95)

c <- plot(pred.3, facets = F, connect.lines = F, colors = "bw", ci.style = "ribbon", dodge = 150) + theme_bw() + xlab("Elevation (m.a.s.l.)") + ylab("Predicted crustose richness")
c.1 <- plot(pred.3, facets = T, connect.lines = F, colors = "bw", ci.style = "ribbon", dodge = 150) + theme_bw() + xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

DF4$pend_c <- factor(DF4$pend_c, levels = c("Flat", "Low_inclination", "High_inclination"), ordered = F)

# GLM for foliose richness
fitfl <- glmer(fl ~ piso*or_R+piso*pend_c+(1|grad), family=poisson(link="log"), data=DF4)
summary(fitfl)
Anova(fitfl)

fitfl2 <- glmer(fl ~ piso+or_R+piso*pend_c+(1|grad), family=poisson(link="log"), data=DF4)
anova(fitfl, fitfl2, test="Chisq")
summary(fitfl2)
Anova(fitfl2)

fitfl3 <- glmer(fl ~ piso*pend_c+(1|grad), family=poisson(link="log"), data=DF4) # final model for foliose richness
anova(fitfl2, fitfl3, test="Chisq")
summary(fitfl3)
Anova(fitfl3)

#Validation
overdisp_fun(fitCr3) # overdispersion

fitted_values <- fitted(fitCr3)
pearson_resid <- residuals(fitCr3, type = "pearson")
residual_data <- data.frame(Fitted = fitted_values, Residuals = pearson_resid)
Vm3 <- ggplot(residual_data, aes(x = Fitted, y = Residuals)) +
  geom_point(color = "black") +  # Scatter plot of residuals
  geom_smooth(method = "lm", color = "blue", se = FALSE) +  # Regression line
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Horizontal line at zero
  labs(title = "C",
       x = "Fitted Values",
       y = "") +
  theme_bw()

# Post hoc
ph_cr2 <- emmeans(fitfl3, ~ piso)
cld_results2 <- cld(ph_cr2, Letters = letters)

ph_cr <- emmeans(fitfl3, ~ piso * pend_c)
cld_results <- cld(ph_cr, Letters = letters)

flat <- subset(cld_results, pend_c == 'Flat')
north <- subset(cld_results, pend_c == 'Low_inclination')
south <- subset(cld_results, pend_c == 'High_inclination')

# Plots for foliose GLM
pred.4=ggpredict(fitfl3, terms = c("piso", "pend_c"), ci.lvl = 0.95)
pred.4=ggpredict(fitfl3, terms = c("piso"), ci.lvl = 0.95)
windows() # GLM
d <- plot(pred.4, facets = F, connect.lines = F, colors = "bw", ci.style = "ribbon", dodge = 150) + theme_bw() + xlab("Elevation (m.a.s.l.)") + ylab("Predicted foliose richness")
d.1 <- plot(pred.4, facets = T, connect.lines = F, colors = "bw", ci.style = "ribbon", dodge = 150) + theme_bw() + xlab("Elevation (m.a.s.l.)") + ylab("Predicted foliose richness")+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

# GLM for fruticose richness
fitfr <- glmer(fr ~ piso*or_R+piso*pend_c+(1|grad), family=poisson(link="log"),  data=DF4)
summary(fitfr)
Anova(fitfr)

fitfr2 <- glmer(fr ~ piso+or_R+piso*pend_c+(1|grad), family=poisson(link="log"), data=DF4)
anova(fitfr2, fitfr, test="Chisq")
summary(fitfr2)
Anova(fitfr2)

fitfr3 <- glmer(fr ~ piso+or_R+pend_c+(1|grad), family=poisson(link="log"), data=DF4) # final model for fruticose richness
anova(fitfr3, fitfr2, test="Chisq")
summary(fitfr3)
Anova(fitfr3)

#Validation
overdisp_fun(fitfr3) # overdispersion

fitted_values <- fitted(fitfr3)
pearson_resid <- residuals(fitfr3, type = "pearson")
residual_data <- data.frame(Fitted = fitted_values, Residuals = pearson_resid)
Vm4 <- ggplot(residual_data, aes(x = Fitted, y = Residuals)) +
  geom_point(color = "black") +  # Scatter plot of residuals
  geom_smooth(method = "lm", color = "blue", se = FALSE) +  # Regression line
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Horizontal line at zero
  labs(title = "D",
       x = "Fitted Values",
       y = "") +
  theme_bw()  # no pattern

# Post hoc
ph_fr <- emmeans(fitfr3, ~ piso)
cld_results <- cld(ph_fr, Letters = letters)

ph_fr2 <- emmeans(fitfr3, ~ or_R)
cld_results2 <- cld(ph_fr2, Letters = letters)

ph_fr3 <- emmeans(fitfr3, ~ pend_c)
cld_results3 <- cld(ph_fr3, Letters = letters)

# GLM plots for fruticose richness
pred.5=ggpredict(fitfr3, terms = c("piso", "pend_c"), ci.lvl = 0.95)
pred.6=ggpredict(fitfr3, terms = c("piso", "or_R"), ci.lvl = 0.95)

pred.7=ggpredict(fitfr3, terms = c("pend_c"), ci.lvl = 0.95)
pred.8=ggpredict(fitfr3, terms = c("or_R"), ci.lvl = 0.95)

e <- plot(pred.5, facets = F, connect.lines = F, colors = "bw", ci.style = "ribbon", dodge = 150) + theme_bw() + xlab("Elevation (m.a.s.l.)") + ylab("Predicted fructicose richness")

f <- plot(pred.6, facets = F, connect.lines = F, colors = "bw", ci.style = "ribbon", dodge = 150) + theme_bw() + xlab("Elevation (m.a.s.l.)") + ylab("Predicted fructicose richness")

windows() # optional
grid.arrange(c,d)
grid.arrange(e,f)

g <- plot(pred.7, facets = F, connect.lines = F, colors = "bw", ci.style = "ribbon") + theme_bw() + xlab("") + ylab("Inclination")

h <- plot(pred.8, facets = F, connect.lines = F, colors = "bw", ci.style = "ribbon") + theme_bw() + xlab("") + ylab("Aspect")

#richness fitness validation plots
grid.arrange(Vm1, Vm2, Vm3, Vm4, ncol=2)

###############################

### Probando paquetes de rasgos ###

library(vegan)

rm(list=ls())
setwd("C:/R/tesis/comunidades_2023")


traits <- read.table("traits.csv",
                    header = T,
                    sep = ";",
                    dec = ".")

datos <- read.table("multivariado.csv",
                    header = T,
                    sep = ";",
                    dec = ".",
                    row.names = 1)

com <- datos[,c(23:158)]
env <- datos[,c(2,3,6,9)]

traits <- traits[,c(1:2)]

# CCA

# cca con variables ambientales específicas
cca1 <- cca(com ~ pend_c + piso + or_R + exp, data = env)
summary(cca1)
anova.cca(cca1)
anova.cca(cca1, by ="terms")

windows()
plot(cca1, display = c("species", "bp"))

### ggplot ###
library(ggplot2)

#CCA
cca_model<-cca(com ~ pend_c + piso + or_R + exp, data = env)
anova.cca(cca_model, by ="terms")

#Get CCA scores
df_species  <- data.frame(summary(cca_model)$species[,1:2])# get the species CC1 and CC2 scores
df_environ  <- scores(cca_model, display = 'bp') #get the environment vars CC1 and CC2 scores
df_environ <- as.data.frame(df_environ)

cca1_varex<-round(summary(cca_model)$cont$importance[2,1]*100,2) #Get percentage of variance explained by first axis
cca2_varex<-round(summary(cca_model)$cont$importance[2,2]*100,2) #Get percentage of variance explained by second axis

#Set a scaling variable to multiply the CCA values, in order to get a very similar plot to the the one generated by plot(cca_model). You can adjust it according to your data
scaling_factor <- 2

windows()
g1 <- ggplot(df_species, 
       aes(x=CCA1, y=CCA2)) + 
  #Draw lines on x = 0 and y = 0
  geom_hline(yintercept=0, 
             linetype="dashed") +
  geom_vline(xintercept=0, 
             linetype="dashed") +
  coord_fixed()+
  #Add species text
  geom_point(data=df_species, shape=17,
            aes(x=CCA1,#Score in CCA1 to add species text
                y=CCA2,#Score in CCA2 to add species text
                label=rownames(df_species),
                hjust=0.5*(1-sign(CCA1)),#Set the text horizontal alignment according to its position in the CCA plot
                vjust=0.5*(1-sign(CCA2))),#Set the text vertical alignment according to its position in the CCA plot
            color = "black")
#Add environmental vars arrows
rownames(df_environ) <- c("Flat rocks", "Vertical rocks","Elevation","Flat rocks (aspect)", "South-facing rocks",
                          "Exposition")

g2 <- g1 +  geom_segment(data=df_environ, 
               aes(x=0, #Starting coordinate in CCA1 = 0 
                   y=0, #Start in CCA2 = 0
                   xend=CCA1*scaling_factor,#Ending coordinate in CCA1  
                   yend=CCA2*scaling_factor), #Ending coordinate in CCA2 
               color="darkred", #set color
               arrow=arrow(length=unit(0.01,"npc"))#Set the size of the lines that form the tip of the arrow
  )
g2

#Add environmental vars text
g3 <- g2 +  geom_text(data=df_environ,
            aes(x=CCA1*scaling_factor, 
                y=CCA2*scaling_factor,
                label=rownames(df_environ),
                hjust=0.2*(1-sign(CCA1)),#Add the text of each environmental var at the end of the arrow
                vjust=0.5*(1-sign(CCA2))),#Add the text of each environmental var at the end of the arrow 
            color="darkred")+
  #Set bw theme
  theme_bw()+
  #Set x and y axis titles
  labs(x=paste0("CCA1 (",cca1_varex," %)"),
       y=paste0("CCA2 (",cca2_varex," %)"))
g3


############################################################

anova.cca(cca1)
anova(cca1, by = "terms") # todas significativas
RsquareAdj(cca1) # 4% de explicacion
summary(cca1)

# Scores de las especies
head(vegan::scores(cca1, display = "species")) 
head(env) # estos scores no se correlacionan directamente con el cca

# cca de cada amb importante + condicionales (covariables), 
#obligamos a cada variable a correlacionarse con el eje1
cca_p <- cca(com ~ pend + Condition(alt + exp + or_R), data = env)
cca_a <- cca(com ~ alt + Condition(pend + exp + or_R), data = env)
cca_e <- cca(com ~ exp + Condition(alt + pend + or_R), data = env)
cca_o <- cca(com ~ or_R + Condition(exp + alt + pend), data = env)

windows()
plot(cca_p, display = c("sp", "bp"))
text(cca_p, display = "species")

plot(cca_a, display = c("species", "bp"))
text(cca_a, display = "species")

plot(cca_e, display = c("species", "bp"))
text(cca_e, display = "species")

plot(cca_o, display = c("species", "bp"))
text(cca_a, display = "species")

####################################################
### CCA alt ###
cca_model <- cca_a
#plot(cca_p,choices=c(1,2), display=c('sp','bp'), scaling=2)

df_species  <- data.frame(summary(cca_model)$species[,1:2])# get the species CC1 and CC2 scores
df_environ  <- scores(cca_model, display = 'bp') #get the environment vars CC1 and CC2 scores
df_environ <- as.data.frame(df_environ)

cca1_varex<-round(summary(cca_model)$cont$importance[2,1]*100,2) #Get percentage of variance explained by first axis
cca2_varex<-round(summary(cca_model)$cont$importance[2,2]*100,2) #Get percentage of variance explained by second axis
scaling_factor <- 2

windows()
g1 <- ggplot(df_species, 
             aes(x=CCA1, y=CA1)) + 
  geom_hline(yintercept=0, 
             linetype="dashed") +
  geom_vline(xintercept=0, 
             linetype="dashed") +
  coord_fixed()+
  geom_point(data=df_species, shape=17,
             aes(x=CCA1,#Score in CCA1 to add species text
                 y=CA1,#Score in CCA2 to add species text
                 label=rownames(df_species),
                 hjust=0.5*(1-sign(CCA1)),#Set the text horizontal alignment according to its position in the CCA plot
                 vjust=0.5*(1-sign(CA1))),#Set the text vertical alignment according to its position in the CCA plot
             color = "black")
#rownames(df_environ) <- c("Rock slope", "Elevation","Flat rocks","South-facing rocks", "Exposition")

g2 <- g1 +  geom_segment(data=df_environ, 
                         aes(x=0, #Starting coordinate in CCA1 = 0 
                             y=0, #Start in CCA2 = 0
                             xend=CCA1*scaling_factor,#Ending coordinate in CCA1  
                             yend=CA1*scaling_factor), #Ending coordinate in CCA2 
                         color="darkred", #set color
                         arrow=arrow(length=unit(0.01,"npc"))#Set the size of the lines that form the tip of the arrow
)
g2

g3 <- g2 +  geom_text(data=df_environ,
                      aes(x=CCA1*scaling_factor, 
                          y=CA1*scaling_factor,
                          label=rownames(df_environ),
                          hjust=0.2*(1-sign(CCA1)),#Add the text of each environmental var at the end of the arrow
                          vjust=0.5*(1-sign(CA1))),#Add the text of each environmental var at the end of the arrow 
                      color="darkred")+
  theme_bw()+
  #Set x and y axis titles
  labs(x=paste0("CCA1 (",cca1_varex," %)"),
       y=paste0("CA1 (",cca2_varex," %)"))
g3


### Plot de cca pendiente ###

cca_model <- cca(com ~ pend + Condition(alt + exp + or_R), data = env)
plot(cca_p,choices=c(1,2), display=c('sp','bp'), scaling=2)

df_species  <- data.frame(summary(cca_model)$species[,1:2])# get the species CC1 and CC2 scores
df_environ  <- scores(cca_model, display = 'bp') #get the environment vars CC1 and CC2 scores
df_environ <- as.data.frame(df_environ)

cca1_varex<-round(summary(cca_model)$cont$importance[2,1]*100,2) #Get percentage of variance explained by first axis
cca2_varex<-round(summary(cca_model)$cont$importance[2,2]*100,2) #Get percentage of variance explained by second axis
scaling_factor <- 2

windows()
g1 <- ggplot(df_species, 
             aes(x=CCA1, y=CA1)) + 
  geom_hline(yintercept=0, 
             linetype="dashed") +
  geom_vline(xintercept=0, 
             linetype="dashed") +
  coord_fixed()+
  geom_point(data=df_species, shape=17,
             aes(x=CCA1,#Score in CCA1 to add species text
                 y=CA1,#Score in CCA2 to add species text
                 label=rownames(df_species),
                 hjust=0.5*(1-sign(CCA1)),#Set the text horizontal alignment according to its position in the CCA plot
                 vjust=0.5*(1-sign(CA1))),#Set the text vertical alignment according to its position in the CCA plot
             color = "black")
#rownames(df_environ) <- c("Rock slope", "Elevation","Flat rocks","South-facing rocks", "Exposition")

g2 <- g1 +  geom_segment(data=df_environ, 
                         aes(x=0, #Starting coordinate in CCA1 = 0 
                             y=0, #Start in CCA2 = 0
                             xend=CCA1*scaling_factor,#Ending coordinate in CCA1  
                             yend=CA1*scaling_factor), #Ending coordinate in CCA2 
                         color="darkred", #set color
                         arrow=arrow(length=unit(0.01,"npc"))#Set the size of the lines that form the tip of the arrow
)
g2

g3 <- g2 +  geom_text(data=df_environ,
                      aes(x=CCA1*scaling_factor, 
                          y=CA1*scaling_factor,
                          label=rownames(df_environ),
                          hjust=0.2*(1-sign(CCA1)),#Add the text of each environmental var at the end of the arrow
                          vjust=0.5*(1-sign(CA1))),#Add the text of each environmental var at the end of the arrow 
                      color="darkred")+
  theme_bw()+
  #Set x and y axis titles
  labs(x=paste0("CCA1 (",cca1_varex," %)"),
       y=paste0("CA1 (",cca2_varex," %)"))
g3

####################################################

#scores "corregidos"
head(vegan::scores(cca_p, choice = 1, display = "species"))
head(vegan::scores(cca_a, choice = 1, display = "species"))
head(vegan::scores(cca_e, choice = 1, display = "species"))
head(vegan::scores(cca_o, choice = 1, display = "species"))


# dendrogramas

dist_CCA <- dist(vegan::scores(cca1, display = "species"))
clust_CCA <- hclust(dist_CCA, method = "ward.D2")
windows()
plot(clust_CCA, cex = 0.6)

library("NbClust")
groups_CCA <- NbClust(diss = dist_CCA, distance = NULL, min.nc = 2, max.nc = 6,
                      method = "ward.D2", index = "silhouette")
groups_CCA$Best.partition

summary(lm(as.matrix(traits) ~ groups_CCA$Best.partition))

summary(lm(traits$cr ~ groups_CCA$Best.partition))
summary(lm(traits$fl ~ groups_CCA$Best.partition))
summary(lm(traits$fr ~ groups_CCA$Best.partition))
summary(lm(traits$di ~ groups_CCA$Best.partition))

# dCCA luego de obtener los scores de las especies agrupadas por sus traits hacemos un cluster para agrupar los rasgos funcionales

dCCA1 <- cca(com ~ pend + alt + or_R + exp, data = env) # no estoy seguro de esto

dist_dCCA <- dist(cca1$co)
clust_dCCA <- hclust(dist_dCCA, method = "ward.D2")
windows()
plot(clust_dCCA, cex = 0.6)

##############################################################################
######################   CWM       ###############
library("FD")
library("vegan")

#citation("FD")

DF <- read.table("multivariado.csv",
                   header = T,
                   sep = ";",
                   dec = ".",
                   row.names = 1) 

envxp <- DF[,c(1:22)]
spxp <- DF[,c(23:158)] 
spxp <- as.matrix(spxp)
spxt <- read.table("traits.csv",
                     header = T, row.names = 1,
                     sep = ";",
                     dec = ".")
spxt <- spxt[1]
row.names(spxt) <- colnames(spxp)

resCWM <- functcomp(spxt, spxp, CWM.type = "all")
head(resCWM)

windows()
g1 <- ggplot(envxp,aes(x = envxp$alt, y = resCWM$mft_Cr, shape=envxp$or_R, linetype=envxp$or_R))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x, color="black", se=F)+
  theme_bw() + labs(x="")

g2 <- ggplot(envxp,aes(x = envxp$alt, y = resCWM$mft_Di, shape=envxp$or_R, linetype=envxp$or_R))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x, color="black", se=F)+
  theme_bw() + labs(x="")

g3 <- ggplot(envxp,aes(x = envxp$alt, y = resCWM$mft_Fl, shape=envxp$or_R, linetype=envxp$or_R))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x, color="black", se=F)+
  theme_bw() + labs(x="")

g4 <- ggplot(envxp,aes(x = envxp$alt, y = resCWM$mft_Fr, shape=envxp$or_R, linetype=envxp$or_R))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x, color="black", se=F)+
  theme_bw()

library("gridExtra")
windows()
grid.arrange(g1, g2, g3, g4, ncol=2)

# GLMs con CWM_mft
DF <- cbind(resCWM,envxp)
summary(DF)

DF$mft_Crt <- with(DF, replace(mft_Cr, mft_Cr == 0.00000000, 0.00000001))
DF$mft_Flt <- with(DF, replace(mft_Fl, mft_Fl == 0.0000, 0.0001))
DF$mft_Dit <- with(DF, replace(mft_Di, mft_Di == 0.00000000, 0.00000001))
DF$mft_Frt <- with(DF, replace(mft_Fr, mft_Fr == 0.00000, 0.00001))
DF$mft_Frt <- with(DF, replace(mft_Frt, mft_Frt == 1.00000, 0.99999))
DF$piso <- as.factor(DF$piso)
DF$pend_c <- factor(DF$pend_c, levels = c("plana", "media", "vertical"), 
                        ordered = F)
DF$or_R <- factor(DF$or_R, levels = c("plana", "N", "S"), 
                    ordered = F)

p1 <- glmmTMB(mft_Crt ~ alt*or_R + alt*pend_c, data = DF, family = beta_family)
summary(p1)
Anova(p1)

p1.1 <- glmmTMB(mft_Crt ~ alt*or_R + alt+pend_c, data = DF, family = beta_family)
anova(p1, p1.1, test="Chisq")
summary(p1.1)
Anova(p1.1)

#####################################################################
p1.2 <- glmmTMB(mft_Crt ~ alt*or_R, data = DF, family = beta_family)
anova(p1.1, p1.2, test="Chisq")
summary(p1.2)
Anova(p1.2)

# Validation
library(DHARMa)
sim_residuals <- simulateResiduals(fittedModel = p1.2)
windows()
plot(sim_residuals)
sim_residuals$

testDispersion(sim_residuals)
testUniformity(sim_residuals)

fitted_values <- fitted(p1.2)
pearson_resid <- residuals(p1.2, type = "pearson")
residual_data <- data.frame(Fitted = fitted_values, Residuals = pearson_resid)
Vm5 <- ggplot(residual_data, aes(x = Fitted, y = Residuals)) +
  geom_point(color = "black") +  # Scatter plot of residuals
  geom_smooth(method = "lm", color = "blue", se = FALSE) +  # Regression line
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Horizontal line at zero
  labs(title = "A",
       x = "",
       y = "Pearson Residuals") +
  theme_bw()

#####################################################################

p2 <- glmmTMB(mft_Flt ~ alt*or_R + alt*pend_c, data = DF, family = beta_family)
summary(p2)
Anova(p2)

p2.1 <- glmmTMB(mft_Flt ~ alt+or_R + alt*pend_c, data = DF, family = beta_family)
anova(p2, p2.1, test="Chisq")
summary(p2.1)
Anova(p2.1)

###########################################################################
p2.2 <- glmmTMB(mft_Flt ~ alt*pend_c, data = DF, family = beta_family)
anova(p2.1, p2.2, test="Chisq")
summary(p2.1)
Anova(p2.2)

#Validation
fitted_values <- fitted(p2.2)
pearson_resid <- residuals(p2.2, type = "pearson")
residual_data <- data.frame(Fitted = fitted_values, Residuals = pearson_resid)
Vm6 <- ggplot(residual_data, aes(x = Fitted, y = Residuals)) +
  geom_point(color = "black") +  # Scatter plot of residuals
  geom_smooth(method = "lm", color = "blue", se = FALSE) +  # Regression line
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Horizontal line at zero
  labs(title = "B",
       x = "",
       y = "") +
  theme_bw()
###########################################################################

p3 <- glmmTMB(mft_Dit ~ alt*or_R +alt*pend_c , data = DF, family = beta_family)
summary(p3)
Anova(p3)

p3.1 <- glmmTMB(mft_Dit ~ alt+or_R +alt*pend_c , data = DF, family = beta_family)
anova(p3, p3.1, test="Chisq")
summary(p3.1)
Anova(p3.1)

#####################################################################
p3.2 <- glmmTMB(mft_Dit ~ alt*pend_c , data = DF, family = beta_family)
anova(p3.1, p3.2, test="Chisq")
summary(p3.2)
Anova(p3.2)

#Validation
fitted_values <- fitted(p3.2)
pearson_resid <- residuals(p3.2, type = "pearson")
residual_data <- data.frame(Fitted = fitted_values, Residuals = pearson_resid)
Vm7 <- ggplot(residual_data, aes(x = Fitted, y = Residuals)) +
  geom_point(color = "black") +  # Scatter plot of residuals
  geom_smooth(method = "lm", color = "blue", se = FALSE) +  # Regression line
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Horizontal line at zero
  labs(title = "C",
       x = "Fitted Values",
       y = "") +
  theme_bw()
#####################################################################


############################################################################
p4 <- glmmTMB(mft_Frt ~ alt*or_R + alt*pend_c, data = DF, family = beta_family)
summary(p4)
Anova(p4)

#Validation
fitted_values <- fitted(p4)
pearson_resid <- residuals(p4, type = "pearson")
residual_data <- data.frame(Fitted = fitted_values, Residuals = pearson_resid)
Vm8 <- ggplot(residual_data, aes(x = Fitted, y = Residuals)) +
  geom_point(color = "black") +  # Scatter plot of residuals
  geom_smooth(method = "lm", color = "blue", se = FALSE) +  # Regression line
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +  # Horizontal line at zero
  labs(title = "D",
       x = "",
       y = "") +
  theme_bw()
############################################################################
grid.arrange(Vm5, Vm6, Vm7, Vm8, ncol=2)

############################################################################

windows()
p1.2
summary(p1.2)
Anova(p1.2)

pred.Cr=ggemmeans(p1.2, terms = c("alt", "or_R"))
a <- plot(pred.Cr, facets = F, 
          connect.lines = F, 
          colors = "bw",
          ci.style = "ribbon", dodge = 150, ci = T) + 
  theme_bw() +
  xlab("Elevation (m.a.s.l.)") +
  ylab("Predicted crustose CWM")
a

pred.Fl=ggemmeans(p2.2, terms = c("alt", "pend_c"), ci.lvl = 0.95)
b <- plot(pred.Fl, facets = F, 
          connect.lines = F, 
          colors = "bw",
          ci.style = "ribbon", dodge = 150) + 
  theme_bw() +
  xlab("Elevation (m.a.s.l.)") +
  ylab("Predicted foliose CWM")

pred.Di=ggemmeans(p3.2, terms = c("alt", "pend_c"), ci.lvl = 0.95)
c <- plot(pred.Di, facets = F, 
          connect.lines = F, 
          colors = "bw",
          ci.style = "ribbon", dodge = 150) + 
  theme_bw() +
  xlab("Elevation (m.a.s.l.)") +
  ylab("Predicted Dimorphic CWM")
c

pred.Fr1=ggemmeans(p4, terms = c("alt", "or_R"), ci.lvl = 0.95)
d <- plot(pred.Fr1, facets = F, 
          connect.lines = F, 
          colors = "bw",
          ci.style = "ribbon", dodge = 150) + 
  theme_bw() +
  xlab("Elevation (m.a.s.l.)") +
  ylab("Fructicose")

pred.Fr2=ggemmeans(p4, terms = c("alt", "pend_c"), ci.lvl = 0.95)
e <- plot(pred.Fr2, facets = F, 
          connect.lines = F, 
          colors = "bw",
          ci.style = "ribbon", dodge = 150) + 
  theme_bw() +
  xlab("Elevation (m.a.s.l.)") +
  ylab("Fructicose")

windows()
e <- plot(pred.Fr, facets = F, 
          connect.lines = F, 
          colors = "bw",
          ci.style = "ribbon", dodge = 150) + 
  theme_bw() +
  xlab("Elevation (m.a.s.l.)") +
  ylab("Predicted fructicose richness")

f <- plot(pred.Fr2, facets = F, 
         connect.lines = F, 
         colors = "bw",
         ci.style = "ribbon", dodge = 150) + 
  theme_bw() +
  xlab("Elevation (m.a.s.l.)") +
  ylab("Predicted fructicose richness")
e
f

windows()
grid.arrange(a, b, c, d, e, ncol=2)
grid.arrange(a, b, ncol=2)
grid.arrange(e, f, ncol=2)
grid.arrange(c, d, e, f, ncol=2)

### RDA ###

rdaNEspain.all <- rda(resCWM ~alt + exp + pend + alt*or_R, data = envxp)
plot(rdaNEspain.all, type = "n", scaling = "sites")
text(rdaNEspain.all, dis = "cn", scaling = "sites")
text(rdaNEspain.all, dis = "sp", scaling = "sites", col = "red")

summary(rdaNEspain.all)

ordistep(rdaNEspain.all, direction = "both")

rdaNEspain.all <- rda(resCWM ~ alt + exp + pend + alt * or_R, data = envxp)

summary(rdaNEspain.all)

anova(rdaNEspain.all)
anova(rdaNEspain.all, by="terms")

RsquareAdj(rdaNEspain.all)

sumrda <- summary(rdaNEspain.all)

sumrda$species
sumrda$sites
sumrda$biplot
sumrda$cont
  
plot(resCWM$mft_Cr ~alt, data=envxp)
windows()

#######################################################################

# RDA CON GGPLOT

##########################################################
library("ggplot2")

df_species  <- data.frame(summary(rdaNEspain.all)$species[,1:2])# get the species CC1 and CC2 scores
df_environ  <- scores(rdaNEspain.all, display = 'bp') #get the environment vars CC1 and CC2 scores
df_environ <- as.data.frame(df_environ)

cca1_varex<-round(summary(rdaNEspain.all)$cont$importance[2,1]*100,2) #Get percentage of variance explained by first axis
cca2_varex<-round(summary(rdaNEspain.all)$cont$importance[2,2]*100,2) #Get percentage of variance explained by second axis
#scaling_factor <- 2
shapes <- c(15,18,16,17)
names(shapes) <- letters[1:4]

windows()
g1 <- ggplot(df_species, 
             aes(x=RDA1, y=RDA2)) + 
  geom_hline(yintercept=0, 
             linetype="dashed") +
  geom_vline(xintercept=0, 
             linetype="dashed") +
  coord_fixed()
  #geom_text(data=df_species, #shape=17,
             #aes(x=RDA1,#Score in CCA1 to add species text
                 #y=RDA2),#Score in CCA2 to add species text
                #label=rownames(df_species))#,
                 #hjust=0.5*(1-sign(rdaNEspain.all)),#Set the text horizontal alignment according to its position in the CCA plot
                 #vjust=0.5*(1-sign(rdaNEspain.all)),#Set the text vertical alignment according to its position in the CCA plot
             #color = "black")
g1
#rownames(df_environ) <- c("Rock slope", "Elevation","Flat rocks","South-facing rocks", "Exposition")

g2 <- g1 +  geom_segment(data=df_environ, 
                         aes(x=0, #Starting coordinate in CCA1 = 0 
                             y=0, #Start in CCA2 = 0
                             xend=RDA1,#Ending coordinate in CCA1  
                             yend=RDA2), #Ending coordinate in CCA2 
                         color="darkred", #set color
                         arrow=arrow(length=unit(0.01,"npc")))#Set the size of the lines that form the tip of the arrow

g3 <- g2 +  geom_text(data=df_environ,
                      aes(x=RDA1, 
                          y=RDA2,
                          label=rownames(df_environ),
                          hjust=0.2*(1-sign(RDA1)),#Add the text of each environmental var at the end of the arrow
                          vjust=0.5*(1-sign(RDA2))),#Add the text of each environmental var at the end of the arrow 
                      color="darkred")+
  theme_bw()+
  #Set x and y axis titles
  labs(x=paste0("RDA1 (",cca1_varex," %)"),
       y=paste0("RDA2 (",cca2_varex," %)"))

shapes = c(16, 17, 6, 15)
g4 <- g3 + geom_point(aes(x = RDA1, y = RDA2), data = df_species, shape = shapes, size= 3)
g4

################################################

############# Analisis con Generos #############

################################################

library("FD")
DF <- read.table("multivariado_generos.csv",
                 header = T,
                 sep = ";",
                 dec = ".",
                 row.names = 1) 

envxp <- DF[,c(1:22)]

spxp <- DF[,c(23:65)]

#spxp <- as.matrix(spxp)

spxp <- spxp * 5
presence_percentage <- (colSums(spxp > 0) / nrow(spxp)) * 100 # Calculate the percentage of presence for each species
species_to_remove <- names(presence_percentage[presence_percentage < 10]) # Identify species with less than 10% presence
spxp <- spxp[, !(colnames(spxp) %in% species_to_remove)] # Remove identified species
spxp <- spxp / 5 

spxt <- read.table("traits_generos.csv",
                   header = T, row.names = 1,
                   sep = ";",
                   dec = ".")

spxt <- spxt[1]
row.names(spxt) <- colnames(spxp)

resCWM <- functcomp(spxt, spxp, CWM.type = "all")

library("vegan")
rdaNEspain.all <- rda(spxp ~alt + exp + pend + alt*or_R, data = envxp)
windows()
plot(rdaNEspain.all, type = "n", scaling = "sites")
text(rdaNEspain.all, dis = "cn", scaling = "sites")
text(rdaNEspain.all, dis = "sp", scaling = "sites", col = "red")

ordistep(rdaNEspain.all, direction = "both")

rdaNEspain.all <- rda(spxp ~ alt + exp + pend + alt * or_R, data = envxp)

summary(rdaNEspain.all)

anova(rdaNEspain.all)
anova(rdaNEspain.all, by="terms")
RsquareAdj(rdaNEspain.all)

#######################################################################

# RDA CON GGPLOT

##########################################################
library("ggplot2")

df_species  <- data.frame(summary(rdaNEspain.all)$species[,1:2])# get the species CC1 and CC2 scores
df_environ  <- scores(rdaNEspain.all, display = 'bp') #get the environment vars CC1 and CC2 scores
df_environ <- as.data.frame(df_environ)

cca1_varex<-round(summary(rdaNEspain.all)$concont$importance[2,1]*100,2) #Get percentage of variance explained by first axis
cca2_varex<-round(summary(rdaNEspain.all)$concont$importance[2,2]*100,2) #Get percentage of variance explained by second axis
#scaling_factor <- 2
shapes <- c(15,18,16,17)
names(shapes) <- letters[1:4]

windows()
g1 <- ggplot(df_species, 
             aes(x=RDA1, y=RDA2)) + 
  geom_hline(yintercept=0, 
             linetype="dashed") +
  geom_vline(xintercept=0, 
             linetype="dashed") +
  coord_fixed()
#geom_text(data=df_species, #shape=17,
#aes(x=RDA1,#Score in CCA1 to add species text
#y=RDA2),#Score in CCA2 to add species text
#label=rownames(df_species))#,
#hjust=0.5*(1-sign(rdaNEspain.all)),#Set the text horizontal alignment according to its position in the CCA plot
#vjust=0.5*(1-sign(rdaNEspain.all)),#Set the text vertical alignment according to its position in the CCA plot
#color = "black")
g1
#rownames(df_environ) <- c("Rock slope", "Elevation","Flat rocks","South-facing rocks", "Exposition")

g2 <- g1 +  geom_segment(data=df_environ, 
                         aes(x=0, #Starting coordinate in CCA1 = 0 
                             y=0, #Start in CCA2 = 0
                             xend=RDA1,#Ending coordinate in CCA1  
                             yend=RDA2), #Ending coordinate in CCA2 
                         color="darkred", #set color
                         arrow=arrow(length=unit(0.01,"npc")))#Set the size of the lines that form the tip of the arrow

g3 <- g2 +  geom_text(data=df_environ,
                      aes(x=RDA1, 
                          y=RDA2,
                          label=rownames(df_environ),
                          hjust=0.2*(1-sign(RDA1)),#Add the text of each environmental var at the end of the arrow
                          vjust=0.5*(1-sign(RDA2))),#Add the text of each environmental var at the end of the arrow 
                      color="darkred")+
  theme_bw()+
  #Set x and y axis titles
  labs(x=paste0("RDA1 (",cca1_varex," %)"),
       y=paste0("RDA2 (",cca2_varex," %)"))

shapes = c(16, 17, 6, 15)
g4 <- g3 + geom_text(aes(x = RDA1+0.05, y = RDA2+0.05), data = df_species, 
                     label=rownames(df_species), size=2.5)
g5 <- g4 + geom_point(aes(x = RDA1, y = RDA2), data = df_species)
g5

df_species

################################################################################
######################################## NMS ###################################
################################################################################

setwd("C:/R/tesis/nms")
rm(list = ls())

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

##############################

g1 <- ggplot(combined_data, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(shape = as.factor(piso), color = as.factor(piso)),size=2) +  # Set shape within aes() and specify size
  #geom_text(data = names, aes(x = NMDS1, y = NMDS2, label = species)) +
  labs(x = "NMDS1", y = "NMDS2", title = "A") +
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

g2

g3 <- g2 +  geom_text(data=df_environ,
                      aes(x=NMDS1, 
                          y=NMDS2,
                          label=rownames(df_environ),
                          hjust=0.2*(1-sign(NMDS1)),#Add the text of each environmental var at the end of the arrow
                          vjust=0.5*(1-sign(NMDS2))),#Add the text of each environmental var at the end of the arrow 
                      color="darkred")  #Set x and y axis titles

windows()
A <- g3

############ B #############
combined_data$piso <- as.factor(combined_data$piso)

g1 <- ggplot(combined_data, aes(x = NMDS1, y = NMDS2)) + 
  #stat_ellipse() +
  labs(x = "NMDS1", y = "NMDS2", title = "B") +
  scale_color_manual(values = elevation_colors, name = "Elevation") +
  guides(linetype = guide_legend(title = "Elevation")) +
  theme_bw()+
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color=piso, linetype=piso))
g1

g2 <- g1 + geom_text(data = names, 
           aes(x = NMDS1, 
               y = NMDS2,
               label = species), size=2.5)
g2
B <- g2

#library("gridExtra")

################################3
datos <- read.table("C:/R/tesis/comunidades_2022/nmds.csv",
                    sep = ";",
                    dec = ".",
                    header = T)

C <- ggplot(datos,aes(x=alt, y=eje1, shape = pend_c, linetype=pend_c))+geom_point()+theme_bw()+
  geom_smooth(method=lm, se=F, color="black")+
  xlab("Altitud (m.s.n.m)")+ ylab("Eje 1 del NMDS")+
  labs(title="C")

grid.arrange(A, B, C, ncol=2)

########### REVISAR SI LOS MT SPP HABITAN OTRAS ALTURAS EN OTRAS TRANSECTAS ################

combined_data$piso <- as.factor(combined_data$piso)

mtop_spp_list <- c("R_geo", "Aca_1", "R_ino", "R_obs", "R_dis",
             "U_sax", "U_lut", "O_plla", "L_rup", "af_aet", 
             "D_dia", "R_cnsa", "L_gran", "A_lzi", "P_war")

mtop_spp <- bio[mtop_spp_list]
mtop_spp <- cbind(mtop_spp, env$piso)

############### GLM ##########
#plot(R_geo/25 ~ as.factor(env$piso), data = mtop_spp)

#library("glmmTMB")
#mtop_spp$R_geo.t <- mtop_spp$R_geo/25
#mtop_spp$R_geo.t <- (mtop_spp$R_geo.t * (length(mtop_spp$R_geo.t) - 1) + 
#                       0.5) / length(mtop_spp$R_geo.t)

#mtop_spp$piso <- as.factor(mtop_spp$`env$piso`)


#R_geo.fit <- glmmTMB(R_geo.t ~ piso, data = mtop_spp, 
#                 family = beta_family(link = "logit"))
#summary(R_geo.fit)
#library("car")
#Anova(R_geo.fit)
#############################
library("gridExtra")
mtop_spp$transect <- sub("_.*", "", row.names(mtop_spp))

mtop_spp$transect <- 
  ifelse(mtop_spp$transect == "T1", "Champaqui",
                            ifelse(mtop_spp$transect == "T2", "Gigantes",
                                   ifelse(mtop_spp$transect == "T3", "Condorito",
                                          ifelse(mtop_spp$transect == "T4", "Linderos", mtop_spp$transect)
                                          )))
windows()
names(mtop_spp)[names(mtop_spp) == "env$piso"] <- "piso"
mtop_spp$piso <- as.factor(mtop_spp$piso)
mtop_spp$transect <- factor(mtop_spp$transect, levels = c("Condorito", "Gigantes", "Champaqui", "Linderos"))


# Rhizocarpon #
geographicum_data1 <- mtop_spp %>%
  group_by(piso) %>%
  summarise(
    R_geo_mean = mean(R_geo, na.rm = TRUE),
    R_geo_sd = sd(R_geo, na.rm = TRUE),
    n = n(),  # number of observations
    R_geo_se = R_geo_sd / sqrt(n)  # standard error
  )

geographicum_data2 <- mtop_spp %>%
  group_by(transect) %>%
  summarise(
    R_geo_mean = mean(R_geo, na.rm = TRUE),
    R_geo_sd = sd(R_geo, na.rm = TRUE),
    n = n(),  # number of observations
    R_geo_se = R_geo_sd / sqrt(n)  # standard error
  )

A <- ggplot(data=geographicum_data1, aes(x=piso, y=R_geo_mean)) + 
  geom_col()+
  geom_errorbar(aes(ymin = R_geo_mean, ymax = R_geo_mean + R_geo_se), width = 0.2) +
  labs(x="Elevation (m.a.s.l.)", y="Rhizocarpon geographicum frequency", title="A")+
  theme_bw()

B <- ggplot(data=geographicum_data2, aes(x=transect, y=R_geo_mean)) + 
  geom_col()+
  geom_errorbar(aes(ymin = R_geo_mean, ymax = R_geo_mean + R_geo_se), width = 0.2) +
  labs(x="Transects", y="Rhizocarpon geographicum frequency", title="B")+
  theme_bw()

#Lorentzii #
A_lzi_data1 <- mtop_spp %>%
  group_by(piso) %>%
  summarise(
    A_lzi_mean = mean(A_lzi, na.rm = TRUE),
    A_lzi_sd = sd(A_lzi, na.rm = TRUE),
    n = n(),  # number of observations
    A_lzi_se = A_lzi_sd / sqrt(n)  # standard error
  )

A_lzi_data2 <- mtop_spp %>%
  group_by(transect) %>%
  summarise(
    A_lzi_mean = mean(A_lzi, na.rm = TRUE),
    A_lzi_sd = sd(A_lzi, na.rm = TRUE),
    n = n(),  # number of observations
    A_lzi_se = A_lzi_sd / sqrt(n)  # standard error
  )

C <- ggplot(data=A_lzi_data1, aes(x=piso, y=A_lzi_mean)) + 
  geom_col()+
  geom_errorbar(aes(ymin = A_lzi_mean, ymax = A_lzi_mean + A_lzi_se), width = 0.2) +
  labs(x="Elevation (m.a.s.l.)", y="Acarospora lorentzii frequency", title="C")+
  theme_bw()

D <- ggplot(data=A_lzi_data2, aes(x=transect, y=A_lzi_mean)) + 
  geom_col()+
  geom_errorbar(aes(ymin = A_lzi_mean, ymax = A_lzi_mean + A_lzi_se), width = 0.2) +
  labs(x="Transects", y="Rhizocarpon geographicum frequency", title="D")+
  theme_bw()
windows()
grid.arrange(A, B, C, D, ncol=2)

# not restricted species: R_cnsa L_gran A_lzi

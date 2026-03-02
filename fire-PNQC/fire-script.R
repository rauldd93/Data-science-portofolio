rm(list = ls())

# libraries
library("glmmTMB")
library("dplyr")
library("car")
library("ggplot2")
library("lattice")
library("FD")
library("vegan")
library("tidyverse")
library("tibble")

#df <- read.table("fire-PNQC/fire-DF.txt",
#                 head = T,
#                 sep=",")

df <- read.table("fire-PNQC/PNQC_df_mft.csv",
                 head = T,
                 sep=",")

head(df)

# Status is not applicable to covers values
df$status <- ifelse(df$spp_lichens %in% c("lichen_cov", "moss_cov", "rock_cov"), 
                    NA, 
                    df$status)

# I need to get the sun exposition, which is the angle in which the sun impact the microsite:
df$exp <- 180 - (df$E + df$W)

# Filter out the cover levels from spp_lichens, not important for richnness
df <- df %>%
  filter(!spp_lichens %in% c("lichen_cov", "moss_cov", "rock_cov"))

df <- df %>%
  filter(!status %in% c("zero"))

head(df)
df <- df %>%
  mutate(status = case_when(
    status %in% c("burned", "damaged") ~ "damaged",
    status %in% c("not_burned", "not_damaged") ~ "not_damaged",
    TRUE ~ status  # Keeps "zero" and any other values as they are
  ))
#levels(as.factor(df$status))

# Exploring full lichen data 
head(df)
boxplot(freq ~ status, data=df)
bwplot(freq ~ status|zone, data=df)

xyplot(freq ~ exp|status, data=df)
xyplot(freq ~ pend|status, data=df)

df <- df %>%
  mutate(pend_c = cut(pend,
                      breaks = c(-Inf, 20, 44, 90),
                      labels = c("flat", "mid", "vert"),
                      include.lowest = TRUE))

plot(freq ~ pend_c, data=df)
bwplot(freq ~ status|pend_c, data=df)

# Convert aspect degrees to ORNS index, -1 = NORTH, 1 = SOUTH

df$ORNS <- cos(df$or_R*pi/180)*-1

plot(freq ~ ORNS, data = df, xlab="ORNS -1=NORTH, 1=SOUTH")

# Create categories for north and south

df <- df %>%
  mutate(or = cut(ORNS,
                  breaks = c(-Inf, 0, Inf),  # Need 3 break points for 2 intervals
                  labels = c("north", "south"),
                  right = FALSE))  # [lower, upper) interval

plot(freq ~ or, data = df)
bwplot(freq ~ status|or, data = df) # aspect does not seems to be relevant

boxplot(freq ~ spp_lichens, data = df)
bwplot(freq ~ spp_lichens|status, data = df) # very hard to see patterns, I believe that using growth forms should be interesting

# Maybe using only Usnea??

df_usn <- df %>%
  filter(spp_lichens=="Usn_amb")

boxplot(freq ~ status, data=df_usn)
bwplot(freq ~ status|zone, data=df_usn, main="Usnea frequency")

bwplot(freq ~ status|pend_c, data=df_usn, main="Usnea frequency")

bwplot(freq ~ status|or, data=df_usn, main="Usnea frequency")

boxplot(freq ~ zone, data = df_usn)

# Model 1: Beta with interactions

# First, convert to proportion (0-1 scale)
df_usn$freq_prop <- df_usn$freq/25 # Assuming max is 25
  
df_usn$freq_prop[df_usn$freq_prop==1] <- 0.999
df_usn$freq_prop[df_usn$freq_prop==0] <- 0.001

summary(df_usn$freq_prop)

# Beta regression (values must be >0 and <1)
# If you have exact 0s or 1s, use zero/one-inflated beta

fit_usn_beta <- glmmTMB(
  freq_prop ~ status + or + pend_c + (1|zone), # Something's wrong with the random effect
  data = df_usn,
  family = beta_family()
)
summary(fit_usn_beta)
Anova(fit_usn_beta)

fit_usn <- lm(freq ~ status*or + status*pend_c, data = df_usn)
summary(fit_usn)
anova(fit_usn)
plot(fit_usn)

########################### RICHNESS #########################

# Calculate richness
# Calculate species richness per rock per grid

#richness_df <- df %>%
#  filter(freq > 0) %>%  # Only include species that were actually observed
#  group_by(nro, grid) %>%
#  summarise(richness = n_distinct(spp_lichens), .groups = "drop")

#head(richness_df)

richness_full <- df %>%
  filter(freq > 0) %>%
  group_by(nro, grid, or_R, pend, E, W, status, zone, obs, exp) %>%
  summarise(richness = n_distinct(spp_lichens), .groups = "drop")

windows()

boxplot(richness ~ status, data = richness_full) # problem with status

levels(as.factor(richness_full$status)) # lets solve this combining levels

richness_full <- richness_full %>%
  mutate(status = case_when(
    status %in% c("burned", "damaged") ~ "damaged",
    status %in% c("not_burned", "not_damaged") ~ "not_damaged",
    TRUE ~ status  # Keeps "zero" and any other values as they are
  ))

levels(as.factor(richness_full$status))

# Quick check
bwplot(richness ~ grid|nro, data = richness_full) # grids tend to be coherent along rocks

boxplot(richness ~ nro, data = richness_full) # in general richness of grids are "stable" within rocks

# are status and zone coherent with each other?

bwplot(richness ~ status|zone, data = richness_full)
# it seems that the damage was higher in not burned sites, 
#showing that lichens are highly sensitive to fire even if the fire did not
#burned those sites

xyplot(richness ~ pend|status, data = richness_full) # most of the damaged thalli is on vertical rocks

xyplot(richness ~ exp|status, data = richness_full)

bwplot(richness ~ status|zone, data = richness_full)

# are they differences between damaged and not damaged thalli in the zones?

fit1 <- lm(richness ~ status*zone, data = richness_full)
summary(fit1)
anova(fit1) # it seems that zone makes no difference in damage

fit1.1 <- lm(richness ~ status+zone, data = richness_full)
summary(fit1.1)
anova(fit1.1) # it seems that zone makes no difference in damage
plot(fit1.1) # let's try to improve this

fit2 <- glmmTMB(richness ~ status * zone + (1|nro/grid), 
                data = richness_full,
                family = poisson) # let's continue with this later

###############################################

head(richness_full)

richness_full <- richness_full %>%
  mutate(pend_c = cut(pend,
                      breaks = c(-Inf, 20, 44, 90),
                      labels = c("flat", "mid", "vert"),
                      include.lowest = TRUE))


boxplot(richness ~ pend, data = richness_full)
boxplot(richness ~ pend_c, data = richness_full)

xyplot(richness ~ pend|zone, data = richness_full)
bwplot(richness ~ pend_c|zone, data = richness_full)

# Convert aspect degrees to ORNS index, -1 = NORTH, 1 = SOUTH

richness_full$or_R
richness_full$or_R*pi/180

richness_full$ORNS <- cos(richness_full$or_R*pi/180)*-1

plot(richness ~ ORNS, data = richness_full, xlab="ORNS -1=NORTH, 1=SOUTH")

# Create categories for north and south

richness_full <- richness_full %>%
  mutate(or = cut(ORNS,
                  breaks = c(-Inf, 0, Inf),  # Need 3 break points for 2 intervals
                  labels = c("north", "south"),
                  right = FALSE))  # [lower, upper) interval

plot(richness ~ or, data = richness_full)
bwplot(richness ~ status|or, data = richness_full)

# Let's get back to the models, but first, I realized that I need to analyze the full data to see damage to individual level

########################## MULTIVARIATE ####################

# Reshape to wide format: rows = samples (nro + grid), columns = species
species_matrix <- df %>%
  select(nro, grid, spp_lichens, freq) %>%
  # Sum frequencies if same species appears multiple times per sample
  group_by(nro, grid, spp_lichens) %>%
  summarise(freq = sum(freq), .groups = "drop") %>%
  pivot_wider(
    names_from = spp_lichens,
    values_from = freq,
    values_fill = 0
  )

# Extract environmental variables
env_data <- df %>%
  distinct(nro, grid, .keep_all = TRUE) %>%
  select(nro, grid, status, zone, pend_c, exp, or) %>%
  mutate(across(c(status, zone, pend_c, or), as.factor))

# Merge to ensure same order
rownames(species_matrix) <- paste(species_matrix$nro, species_matrix$grid, sep = "_")
species_matrix <- as.data.frame(species_matrix[, -c(1:2)])  # Remove nro and grid columns

rownames(env_data) <- paste(env_data$nro, env_data$grid, sep = "_")
env_data <- env_data[, -c(1:2)]  # Remove nro and grid columns

# Permanova

# Calculate Bray-Curtis dissimilarity
dist_matrix <- vegdist(species_matrix, method = "bray")

# PERMANOVA - test which factors affect community composition
adonis_result <- adonis2(dist_matrix ~ status + zone + pend_c + exp + or, 
                         data = env_data,
                         permutations = 999)
print(adonis_result)

# Check for homogeneity of variances (important assumption)
dispersion <- betadisper(dist_matrix, group = env_data$status)
anova(dispersion)  # Should not be significant, 0.02 sig

# NMDS

# NMDS visualization
nmds <- metaMDS(species_matrix, distance = "bray", k = 2, trymax = 20)

# Plot NMDS
windows()
plot(nmds, type = "n", main = "NMDS of Lichen Communities")
points(nmds, display = "sites", 
       col = ifelse(env_data$zone == "burned", "red", "blue"),
       pch = ifelse(env_data$status == "damaged", 16, 17),
       cex = 1.2)

# Add environmental vectors
env_fit <- envfit(nmds ~ status + zone + pend_c + exp + or, 
                  data = env_data, permutations = 999)
plot(env_fit, col = "darkgreen", p.max = 0.05)  # Only significant vectors

# Add legend
legend("topright", 
       legend = c("Burned-Damaged", "Burned-Not Damaged", 
                  "Not Burned-Damaged", "Not Burned-Not Damaged"),
       col = c("red", "red", "blue", "blue"),
       pch = c(16, 17, 16, 17),
       cex = 0.8)

# RDA #

# Hellinger transformation for species data (recommended for RDA)
species_hell <- decostand(species_matrix, method = "hellinger")

# RDA with all environmental variables
rda_full <- rda(species_hell ~ status + zone + pend_c + exp + or, 
                data = env_data)

# Summary
summary(rda_full, scaling = 2)

# Variance inflation factors (check for multicollinearity)
vif.cca(rda_full)  # Should be <10

# Plot RDA
plot(rda_full, type = "n", main = "RDA of Lichen Communities")
points(rda_full, display = "sites", 
       col = ifelse(env_data$zone == "burned", "red", "blue"),
       pch = ifelse(env_data$status == "damaged", 16, 17))
text(rda_full, display = "species", col = "darkgreen", cex = 0.7)
text(rda_full, display = "bp", col = "black", cex = 0.8)

# Permutation test
anova(rda_full, permutations = 999)
anova(rda_full, by = "term", permutations = 999)  # Test each term
# status, zone pend or are significant, or is only 0.05

# Improve the RDA plot
plot(rda_full, type = "n", main = "RDA of Lichen Communities")

# Add points with colors and shapes for zone × status combinations
points(rda_full, display = "sites", 
       col = ifelse(env_data$zone == "burned", "black", "darkgrey"),
       pch = ifelse(env_data$status == "damaged", 16, 17),
       cex = 1.2)

# Add species labels
text(rda_full, display = "species", col = "black", cex = 0.7)

# Add environmental vectors
text(rda_full, display = "bp", col = "darkred", cex = 0.8)

# Add comprehensive legend
legend("topright", 
       legend = c("Burned - Damaged", "Burned - Not Damaged",
                  "Not Burned - Damaged", "Not Burned - Not Damaged",
                  "Environmental Vectors", "Species Labels"),
       col = c("black", "black", "darkgrey", "darkgrey", "darkred", "black"),
       pch = c(16, 17, 16, 17, NA, NA),
       lty = c(NA, NA, NA, NA, 1, NA),
       lwd = c(NA, NA, NA, NA, 2, NA),
       pt.cex = c(1.2, 1.2, 1.2, 1.2, NA, NA),
       cex = 0.8,
       bg = "white", box.col = "gray")

###### Growth forms and Functional Traits ######

head(df)

df$site_id <- paste(df$nro, df$grid, sep = "_")

community_matrix <- df %>%
  select(site_id, spp_lichens, freq) %>%
  pivot_wider(
    names_from = spp_lichens,
    values_from = freq,
    values_fill = 0
  )

df %>%
  count(site_id, spp_lichens) %>%
  filter(n > 1)


#errores:
#> df %>%
#  +   count(site_id, spp_lichens) %>%
#  +   filter(n > 1)
#site_id spp_lichens n
# Fixed?




g1 <- ggplot(envxp, aes(x = zone, 
                        y = resCWM$mft_crustose, 
                        fill = status)) +
  geom_boxplot() +
  theme_bw() + 
  labs(x = "", y = "Crustose CWM") +
  scale_fill_manual(values = c("damaged" = "red", "not_damaged" = "gray"))

g2 <- ggplot(envxp, aes(x = zone, 
                        y = resCWM$mft_foliose, 
                        fill = status)) +
  geom_boxplot() +
  theme_bw() + 
  labs(x = "", y = "Foliose CWM") +
  scale_fill_manual(values = c("damaged" = "red", "not_damaged" = "gray"))

g3 <- ggplot(envxp, aes(x = zone, 
                        y = resCWM$mft_fruticose, 
                        fill = status)) +
  geom_boxplot() +
  theme_bw() + 
  labs(x = "", y = "Fruticose CWM") +
  scale_fill_manual(values = c("damaged" = "red", "not_damaged" = "gray"))


#library("gridExtra")
grid.arrange(g1, g2, g3, ncol=2)

#####################################

# GLMs con CWM_mft
DF <- cbind(resCWM,envxp)
summary(DF)

DF$crt    <- with(DF, replace(mft_crustose, mft_crustose == 0, 0.001))
DF$crt <- with(DF, replace(mft_crustose, mft_crustose == 1, 0.99))
DF$flt <- with(DF, replace(mft_foliose, mft_foliose == 0, 0.001))
DF$frtt <- with(DF, replace(mft_fruticose, mft_fruticose  == 0, 0.001))

#DF$pend_c <- factor(DF$pend_c, levels = c("plana", "media", "vertical"), 
     #               ordered = F)
#DF$or_R <- factor(DF$or_R, levels = c("plana", "N", "S"), 
                #  ordered = F)
summary(DF)

p1 <- glmmTMB(crt ~ status*zone + status*pend_c + status*or, 
              data = DF, family = beta_family)
summary(p1)
Anova(p1)

p1.1 <- glmmTMB(crt ~ status*zone + status*pend_c, 
                data = DF, family = beta_family)
anova(p1, p1.1, test="Chisq")
summary(p1.1)
Anova(p1.1)

p1.2 <- glmmTMB(crt ~ status+zone + status*pend_c, 
                data = DF, family = beta_family)
anova(p1.1, p1.2, test="Chisq")
summary(p1.2)
Anova(p1.2)

p1.3 <- glmmTMB(crt ~ status+zone + pend_c, 
                data = DF, family = beta_family)
anova(p1.3, p1.2, test="Chisq")
summary(p1.3)
Anova(p1.3)

library("ggeffects")
library("emmeans")

pred.Cr=ggemmeans(p1.3, terms = c("zone", "status"))
pred.Cr2=ggemmeans(p1.3, terms = c("zone", "pend_c"))

a <- plot(pred.Cr) + theme_bw() +
  xlab("Zone") +
  ylab("Predicted crustose CWM")

b <- plot(pred.Cr2) + theme_bw() +
  xlab("Zone") +
  ylab("Predicted crustose CWM")


p2 <- glmmTMB(flt ~ status*zone + status*pend_c, 
                data = DF, family = beta_family)
summary(p2)
Anova(p2)

p2.1 <- glmmTMB(flt ~ status*zone, 
              data = DF, family = beta_family)
anova(p2.1, p2, test="Chisq")
summary(p2.1)
Anova(p2.1)

p2.2 <- glmmTMB(flt ~ status+zone, 
                data = DF, family = beta_family)
anova(p2.2, p2.1, test="Chisq")
summary(p2.2)
Anova(p2.2)


########################

# Extract functional diversity indices
functional_indices <- data.frame(
  site = rownames(fd_results$CWM),
  FRic = fd_results$FRic,     # Functional Richness
  FEve = fd_results$FEve,     # Functional Evenness
  FDiv = fd_results$FDiv,     # Functional Divergence
  RaoQ = fd_results$RaoQ      # Rao's quadratic entropy
)

functional_indices

# Add environmental data
functional_indices <- functional_indices %>%
  left_join(env_matrix %>% rownames_to_column("site"), by = "site")

# Community Weighted Means (CWM) for growth forms
cwm_mft <- fd_results$CWM %>%
  as.data.frame() %>%
  rownames_to_column("site")

# Functional diversity by status and zone

# 1. Test differences in functional diversity indices
functional_indices

windows()
boxplot(FRic ~ status, data = functional_indices)


library(lme4)
library(lmerTest)

# Functional Richness
mod_fric <- lmer(FRic ~ status * zone + (1|pend_c) + (1|or), 
                 data = functional_indices)
summary(mod_fric)
anova(mod_fric)

# Functional Evenness
mod_feve <- lmer(FEve ~ status * zone + (1|pend_c) + (1|or), 
                 data = functional_indices)
summary(mod_feve)

# Rao's Q (functional diversity)
mod_raoq <- lmer(RaoQ ~ status * zone + (1|pend_c) + (1|or), 
                 data = functional_indices)
summary(mod_raoq)

# 2. Visualize
library(ggplot2)

# Functional Richness plot
ggplot(functional_indices, aes(x = zone, y = FRic, fill = status)) +
  geom_boxplot() +
  labs(title = "Functional Richness by Zone and Damage Status",
       y = "Functional Richness (FRic)") +
  theme_minimal() +
  scale_fill_manual(values = c("damaged" = "orange", "not_damaged" = "steelblue"))

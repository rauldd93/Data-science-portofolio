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
library("ggeffects")
library("emmeans")
library("viridis")


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

#plot(freq ~ pend_c, data=df)
#bwplot(freq ~ status|pend_c, data=df)

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

## Usnea as model ##

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

windows()
bwplot(freq ~ status|mft, data = df)
bwplot(freq ~ status|reprod, data = df)

df$site_id <- paste(df$nro, df$grid, sep = "_")


## Create species matrix (sites × species)
community_matrix <- df %>%
  # Select relevant columns
  select(site_id, spp_lichens, freq) %>%
  # Pivot to wide format
  pivot_wider(
    id_cols = site_id,
    names_from = spp_lichens,
    values_from = freq,
    values_fill = 0
  ) %>%
  # Convert to matrix (remove site_id column for analysis)
  column_to_rownames("site_id") %>%
  as.matrix()

# Create environmental matrix (same order as community matrix)
env_matrix <- df %>%
  select(site_id, status, pend, exp, or, mft, zone) %>%
  distinct(site_id, .keep_all = TRUE) %>%
  column_to_rownames("site_id")

# Check alignment
all(rownames(community_matrix) == rownames(env_matrix))

# Calculate species richness and diversity
env_matrix$richness <- specnumber(community_matrix)
env_matrix$shannon <- diversity(community_matrix, index = "shannon")
env_matrix$simpson <- diversity(community_matrix, index = "simpson")

# Test differences in diversity
t.test(richness ~ status, data = env_matrix)
windows()
boxplot(richness ~ status, data = env_matrix, main = "Species Richness by Fire Status")
boxplot(shannon ~ status, data = env_matrix, main = "Species Richness by Fire Status")
boxplot(simpson ~ status, data = env_matrix, main = "Species Richness by Fire Status")

# NMS

# NMDS with Bray-Curtis dissimilarity
set.seed(123)  # for reproducibility
nmds <- metaMDS(community_matrix, 
                distance = "bray",
                k = 2,  # number of dimensions
                trymax = 100)

# Check stress (should be < 0.2)
nmds$stress

# Plot NMDS
plot(nmds, type = "n")
points(nmds, display = "sites", pch = 19, 
       col = ifelse(env_matrix$status == "burned", "red", "blue"))
text(nmds, display = "species", cex = 0.7, col = "darkgreen")
legend("topright", legend = c("Burned", "Damaged"), 
       col = c("red", "blue"), pch = 19)

# Add environmental vectors
ef <- envfit(nmds, env_matrix[, c("pend", "mft", "richness")], permutations = 999)
plot(ef, add = TRUE, col = "purple")

# CCA with environmental variables
cca_model <- cca(community_matrix ~ status + pend + exp + or, 
                 data = env_matrix)

# Check results
summary(cca_model)
anova(cca_model)  # Overall test
anova(cca_model, by = "term")  # Test each term
anova(cca_model, by = "axis")  # Test each axis

# Plot CCA
plot(cca_model, type = "n")
points(cca_model, display = "sites", pch = 19,
       col = ifelse(env_matrix$status == "burned", "red", "blue"))
text(cca_model, display = "species", cex = 0.7, col = "darkgreen")
text(cca_model, display = "bp", col = "purple", cex = 0.8)  # environmental variables

# Hellinger-transform species data (recommended for RDA)
com_hel <- decostand(community_matrix, method = "hellinger")

# RDA with environmental variables
rda_model <- rda(com_hel ~ status + pend + exp + or + mft, 
                 data = env_matrix)

# Variance partitioning
anova(rda_model)
RsquareAdj(rda_model)

# Plot
plot(rda_model, type = "n")
points(rda_model, display = "sites", pch = 19,
       col = ifelse(env_matrix$status == "burned", "red", "blue"))
text(rda_model, display = "species", cex = 0.7, col = "darkgreen")
text(rda_model, display = "bp", col = "purple", cex = 0.8)

# Create trait matrix (if you have multiple traits)
# For now, using mft as main trait
trait_matrix <- df %>%
  select(spp_lichens, mft) %>%
  distinct(spp_lichens, .keep_all = TRUE) %>%
  column_to_rownames("spp_lichens")

# Create a comprehensive visualization
library(ggplot2)
library(patchwork)

# NMDS scores for plotting
nmds_scores <- as.data.frame(scores(nmds)$sites)
nmds_scores$status <- env_matrix$status
nmds_scores$pend <- env_matrix$pend

p1 <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = status)) +
  geom_point(size = 3) +
  stat_ellipse() +
  theme_minimal() +
  ggtitle("Community Composition by Fire Status")

p2 <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = pend)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  ggtitle("Community Composition by Slope")

p1 + p2

# Fourth-corner analysis
# First, ensure all matrices align
com_order <- community_matrix[, rownames(trait_matrix)]

env <- as.data.frame(env_matrix[, c("status", "pend", "exp", "or")])

# First, let's check and prepare our data properly



# 2. Convert environmental data to data frame with correct types
env_df <- env_matrix[, c("status", "pend", "exp", "or")] %>%
  as.data.frame() %>%
  mutate(
    # Convert categorical variables to factors
    status = as.factor(status),
    exp = as.numeric(exp),
    or = as.factor(or),
    # Keep pend as numeric
    pend = as.numeric(pend)
  )

# 3. Ensure trait matrix is a data frame with species as rows
trait_df <- trait_matrix %>%
  as.data.frame() %>%
  mutate(mft = as.factor(mft))  # ensure mft is numeric

# 4. Ensure community matrix is numeric and has species as columns
# (it should already have this structure from our earlier step)
com_df <- as.data.frame(com_order)

# 5. Verify all species names match between community and trait matrices
all(colnames(com_df) %in% rownames(trait_df))
all(rownames(trait_df) %in% colnames(com_df))

# If they don't match perfectly, align them:
common_spp <- intersect(colnames(com_df), rownames(trait_df))
com_df <- com_df[, common_spp]
trait_df <- trait_df[common_spp, , drop = FALSE]

# 6. Now try fourth-corner again
library(ade4)  # sometimes fourthcorner comes from ade4

fourth_ade4 <- fourthcorner(
  env_df,           # environmental variables
  com_df,           # community data (sites × species)
  trait_df,         # trait data (species × traits)
  modeltype = "3",
  nrepet = 999
)

summary(fourth_ade4)
plot(fourth_ade4)

################ FD  #########

# Calculate CWM of mft for each site

head(df)

df$site_id <- paste(df$nro, df$grid, sep = "_")

# Create community matrix with sites as rows and species as columns
community_matrix <- df %>%
  # If you have multiple samples per site, you might want to aggregate
  group_by(site_id, spp_lichens) %>%
  summarise(abundance = sum(freq), .groups = "drop") %>%
  pivot_wider(names_from = spp_lichens, 
              values_from = abundance, 
              values_fill = 0) %>%
  column_to_rownames("site_id")

# View the community matrix
head(community_matrix[, 1:10]) # first 10 species

# Create trait matrix (species × traits)
# First, get unique species-trait combinations
trait_data <- df %>%
  select(spp_lichens, mft, reprod) %>%
  distinct() %>%
  # If mft is categorical, you might want to keep it as factor
  mutate(mft = as.factor(mft),
         reprod = as.factor(reprod)) %>%
  column_to_rownames("spp_lichens")

env_matrix <- df %>%
  select(site_id, status, pend, exp, or, zone) %>%
  distinct(site_id, .keep_all = TRUE) %>%
  column_to_rownames("site_id")
head(env_matrix)

trait_data

all(colnames(community_matrix) == rownames(trait_data))
community_matrix <- community_matrix[,match(rownames(trait_data), colnames(community_matrix))]
all(colnames(community_matrix) == rownames(trait_data))

trait_data <- as.matrix(trait_data)
community_matrix <- as.matrix(community_matrix)

resCWM <- functcomp(trait_data, community_matrix, CWM.type = "all")
head(resCWM)

colnames(env_matrix)

windows()
# Crustose CWM
g1 <- ggplot(env_matrix,aes(x = zone, 
                            y = resCWM$mft_crustose, 
                            fill=status))+
  geom_boxplot()+
  theme_bw() + labs(x="") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=1, alpha=0.5, width = 0.1)
g1

g1 <- ggplot(env_matrix,aes(x = pend, 
                            y = resCWM$mft_crustose, 
                            shape=status))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x, se=T)+
  theme_bw() + labs(x="") +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=1, alpha=0.5, width = 0.1)
g1

# Add to env_matrix
env_matrix$cwm_mft <- cwm$mft

env_matrix$cwm_mft
env_matrix$status

# Model CWM response
cwm_model <- lm(status ~ cwm_mft + pend + exp + zone, data = env_matrix)
summary(cwm_model)




#################################################################################

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

git config --global user.email "raulenriquedd@hotmail.com"
git config --global user.name "rauldd"

rm(list = ls())


# libraries
library("glmmTMB")
library("dplyr")
library("car")

df <- read.table("fire-PNQC/fire-DF.txt",
                 head = T,
                 sep=",")
head(df)

# Status is not applicable to covers values
df$status <- ifelse(df$spp_lichens %in% c("lichen_cov", "moss_cov", "rock_cov"), 
                    NA, 
                    df$status)
summary(df)

# I need to get the sun exposition, which is the angle in which the sun impact the microsite:
df$exp <- 180 - (df$E + df$W)
summary(df$exp)

head(df)

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
levels(as.factor(df$status))

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

# let's try a model with Usnea data

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
library(lattice)

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

library("glmmTMB")

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



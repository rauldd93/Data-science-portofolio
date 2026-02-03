rm(list = ls())

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

###### RICHNESS ########
library(dplyr)

head(df)

# Filter out the cover levels from spp_lichens, not important for richnness
df <- df %>%
  filter(!spp_lichens %in% c("lichen_cov", "moss_cov", "rock_cov"))

head(df)

# Calculate richness
# Calculate species richness per rock per grid

richness_df <- df %>%
  filter(freq > 0) %>%  # Only include species that were actually observed
  group_by(nro, grid) %>%
  summarise(richness = n_distinct(spp_lichens), .groups = "drop")

head(richness_df)

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

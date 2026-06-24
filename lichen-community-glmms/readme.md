# Lichen Community Analysis Along an Elevational Gradient

This repository contains the data processing, statistical modeling, and visualization code for the paper:

**[10.5252/cryptogamie-mycologie2026v47a3](http://cryptogamie.com/mycologie/47/3)**

## Analysis Overview

The script `lichen-community-analysis.R` performs all statistical analyses described in the paper, including:

- Data cleaning and processing
- Generalized Linear Mixed Models (GLMMs) for species richness
- Model selection and validation (VIF, diagnostic plots)
- Community Weighted Mean (CWM) trait calculations
- Beta regression GLMs for CWM analyses
- Type II Wald chi-square tests for significant terms
- Post-hoc comparisons using estimated marginal means
- Multivariate analyses (Canonical Correspondence Analysis, Non-metric Multidimensional Scaling)
- Permutational tests for CCA significance
- Canonical correlation analyses
- Result visualization and plotting

## Required Data Files

The following files must be present in the working directory:

| File | Description |
|------|-------------|
| `dat_comunidades.csv` | Main dataset with lichen species frequencies and environmental variables |
| `traits.csv` | Species trait data for Community Weighted Mean (CWM) calculations |
| `multivariado.csv` | Species matrix for multivariate analyses |

---

## Statistical Methods

All analyses were conducted in R (R Core Team 2020), with additional ordination analyses performed using PC-ORD (McCune & Mefford 2011).

### Species Richness Analysis

**Model:** Generalized Linear Mixed Model (GLMM) with Poisson distribution.

- **Response variable:** Lichen species richness
- **Fixed effects:** Elevation, rock aspect, rock slope
- **Random effects:** Sites, with rocks nested within sites (Inchausti 2023)
- **Interaction effects** between elevation, aspect, and slope were evaluated

**Model selection and validation:**
- Model selection was performed comparing candidate models
- Multicollinearity was assessed using Variance Inflation Factor (VIF)
- Model fitness was evaluated through diagnostic indicators (see Appendices 3 and 4 in the paper)
- For flat rocks, aspect was categorized as "flat" with no directional aspect considered, as supported by model selection

**Post-hoc tests:** Conducted using the `emmeans` package to identify differences between explanatory variable levels.

### Community Weighted Mean (CWM) Analysis

**CWM calculation:**
- A frequency matrix was constructed for each growth form (crustose, foliose, fruticose)
- CWMs weighted by relative abundance were calculated using the `FD` package (Laliberté & Legendre 2010)

**Model:** Beta regression GLMs.

- **Response variables:** CWM of each growth form (crustose, foliose, fruticose)
- **Fixed effects:** Elevation, rock aspect, rock inclination

**Significance testing:**
- Type II Wald chi-square tests were performed using the `car` package (Fox & Weisberg 2019)
- Results were plotted using the `ggeffects` package (Lüdecke 2018) and `ggplot2` (Wickham 2016)

**Post-hoc tests:** Applied to GLMs and GLMMs using the `emmeans` package.

### Multivariate Analyses

**Data preparation:**
- Species with less than 10% occurrence were excluded to reduce noise and improve interpretability (McCune & Grace 2003)

#### Canonical Correspondence Analysis (CCA)

- **Objective:** Examine relationships between community structure and environmental variables
- **Environmental variables:** Elevation, rock inclination, solar exposure, rock aspect
- **Significance testing:** Permutational test under a reduced model
- **Species–environment associations:** For each environmental vector, the five species with the smallest angular difference relative to the vector's direction were identified as most strongly associated (as visualized in Fig. 8 of the paper)
- Performed using the `vegan` package (Oksanen et al. 2022)

#### Non-metric Multidimensional Scaling (NMS)

- **Distance matrix:** Bray-Curtis
- **Iterations:** 100 runs, with stress levels stabilizing at iteration 40
- **Final configuration:** Three-axis solution
- The third axis was excluded from plots as no clear pattern emerged between the third and first or second axes
- **Environmental correlations:** Canonical correlation analyses were conducted between NMS axes and environmental variables to assess relationships with the first two axes
- Performed using PC-ORD (McCune & Mefford 2011) and the `vegan` package in R (Oksanen et al. 2022)
- Plotted using `ggplot2` (Wickham 2016)

**Note:** Rock inclination was treated as a quantitative variable in the NMS analysis.

---

## File Descriptions

### `dat_comunidades.csv`

Contains lichen species frequency data and associated environmental variables measured along an elevational gradient in a mountain system of central Argentina.

**Sampling design:** Frequency was measured on 424 rocks using five 20 × 20 cm grids placed equidistantly on each rock. Each grid was subdivided into 25 sub-squares (16 cm² each). Frequency was determined by quantifying the number of sub-squares occupied by each lichen species.

#### Column Descriptions

| Column | Description |
|--------|-------------|
| `Nro_roca` | Rock identification number |
| `ud_muestral` | Unique code indicating transect, elevational level, and rock number |
| `piso` | Elevational level (factor; ±100 m a.s.l. error) |
| `grad` | Transect name |
| `spp_liquenes` | Lichen species name |
| `genus` | Genus of the species |
| `morfotipo` | Growth form/morphotype: `cr` = crustose, `fr` = fruticose, `fl` = foliose |
| `fr1` – `fr5` | Frequency measurements for each of the 5 replicate grids |
| `fr_prom` | Mean frequency across the 5 replicates |
| `or_roca` | Rock face orientation relative to north (0°–360°) |
| `ORNS_roca` | Northness index of the rock face (−1 = fully south-facing, 1 = fully north-facing) |
| `or_R` | Rock face aspect category (`North`, `South`) |
| `exp` | Sun exposure estimate at the microsite, calculated as 180° minus the clinometer-measured angle to the horizon (east and west) |
| `sup` | Rock surface area (m²) |
| `pend` | Rock slope (°) |
| `pend_c` | Slope category: `flat` (0°–20°), `low_inclination` (22°–45°), `high_inclination` (46°–90°) |
| `tipo_roca` | Rock type (`granite` or `gneiss`) |
| `or_ladera` | Hillside aspect |
| `ORNS_lad` | Hillside northness index (−1 = south, 1 = north) |
| `or_L` | Hillside aspect category |
| `paj` | Grass cover (%) within 5 m radius |
| `rock` | Rock cover (%) within 5 m radius |
| `arb` | Shrub cover (%) within 5 m radius |
| `so` | Bare soil cover (%) within 5 m radius |
| `mg` | Moss cover (%) within 5 m radius |
| `rast` | Creeping vegetation cover (%) within 5 m radius |
| `ces` | Grass/sedge cover (%) within 5 m radius |
| `alt` | Elevation above sea level (m) |
| `obs` | Additional observations |

---

### `traits.csv`

Species-level trait data used to calculate Community Weighted Means (CWM) for each growth form via the `FD` package.

| Column | Description |
|--------|-------------|
| `spp` | Lichen species name |
| `mft` | Growth form: `cr` = crustose, `fl` = foliose, `fr` = fruticose |
| `genus` | Genus of the species |

---

### `multivariado.csv`

Community matrix for multivariate analyses (CCA and NMS), with species as columns and environmental variables included. Species with less than 10% occurrence were excluded prior to analysis.

| Column | Description |
|--------|-------------|
| `muestra` | Unique rock identifier (transect, elevational level, rock) |
| `grad` | Transect name |
| `or_R` | Rock face aspect category (`North`, `South`) |
| `pend_c` | Slope category: `flat` (0°–20°), `low_inclination` (22°–45°), `high_inclination` (46°–90°) |
| `tipo_roca` | Rock type (`granite` or `gneiss`) |
| `or_L` | Hillside aspect category |
| `exp` | Sun exposure estimate (see above) |
| `sup` | Rock surface area (m²) |
| `pend` | Rock slope (°) — treated as quantitative in NMS |
| `piso` | Elevational level (factor; ±100 m a.s.l. error) |
| `ORNS_L` | Hillside northness index (−1 = south, 1 = north) |
| `ORNS_R` | Rock face northness index (−1 = south-facing, 1 = north-facing) |
| `paj` | Grass cover (%) within 5 m radius |
| `rock` | Rock cover (%) within 5 m radius |
| `arb` | Shrub cover (%) within 5 m radius |
| `so` | Bare soil cover (%) within 5 m radius |
| `mg` | Moss cover (%) within 5 m radius |
| `rast` | Creeping vegetation cover (%) within 5 m radius |
| `ces` | Grass/sedge cover (%) within 5 m radius |
| `Lichen_cover` | Estimated total lichen cover on each grid (%) |
| `Moss_cover` | Estimated total moss cover on each grid (%) |
| `Rock_cover` | Estimated total rock cover on each grid (%) |
| `A_altoandina` – `Xanthoria_sp` | Individual lichen species columns (frequency values) |

---

## Software and Packages

| Software/Package | Version/Purpose | Reference |
|------------------|-----------------|-----------|
| R | 4.0.3 | R Core Team (2020) |
| `FD` | CWM calculations | Laliberté & Legendre (2010) |
| `car` | Type II Wald chi-square tests | Fox & Weisberg (2019) |
| `ggeffects` | Marginal effects plots | Lüdecke (2018) |
| `ggplot2` | Data visualization | Wickham (2016) |
| `emmeans` | Post-hoc comparisons | — |
| `vegan` | CCA and NMS ordinations | Oksanen et al. (2022) |
| PC-ORD | NMS ordination | McCune & Mefford (2011) |

---

## References

- Fox, J. & Weisberg, S. (2019). *An R Companion to Applied Regression.*
- Laliberté, E. & Legendre, P. (2010). A distance-based framework for measuring functional diversity from multiple traits. *Ecology*, 91(1), 299–305.
- Lüdecke, D. (2018). ggeffects: Tidy data frames of marginal effects from regression models. *Journal of Open Source Software*, 3(26), 772.
- Oksanen, J. et al. (2022). *vegan: Community Ecology Package.*
- R Core Team (2020). *R: A Language and Environment for Statistical Computing.*
- Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis.*

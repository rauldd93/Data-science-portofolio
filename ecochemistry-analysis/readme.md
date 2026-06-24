# Usnic Acid Concentration Analysis in *Usnea amblyoclada*

This repository contains the data and statistical analyses for the results published in:

**10.1007/s11557-023-01904-6**

## Analysis Overview

The script `ecochemistry-script.R` performs all statistical tests, analyses, and visualizations described in the paper, including:

- Descriptive statistics (median, standard deviation) of extract yield and usnic acid concentration
- Generalized Linear Models (GLMs) with Gaussian and Beta distributions
- Post-hoc comparisons (Tukey test and pairwise estimated marginal means)
- Model diagnostics (residuals vs. fitted, normal Q-Q, scale-location plots)
- Principal Component Analysis (PCA) incorporating water retention traits and environmental conditions

## Required Data

The file `cc_usnico2.csv` must be present in the working directory to run the analysis.

---

## File Description

### `cc_usnico2.csv`

Contains usnic acid concentration data from thalli of the lichen *Usnea amblyoclada* collected along an elevational gradient in a mountain system. Samples were taken from thalli growing on rocks under different microsite conditions. Total extract of secondary compounds was also measured as complementary data.

**Sampling design:** Each sample represents a composite of thalli collected at a specific elevation and microsite condition.

#### Column Descriptions

| Column | Description |
|--------|-------------|
| `muestra` | Sample ID number (a composite sample of thalli at a specific elevation and microsite condition) |
| `nro_uplc` | Sample number assigned by the UPLC operator |
| `area_prom` | Area under the curve detected by UPLC |
| `usnico` | Usnic acid concentration (% usnic acid per dry mass of thallus) |
| `peso_liquen` | Dry weight of the thallus/sample (mg) |
| `extracto` | Weight of the total extract obtained from each sample, including compounds other than usnic acid (mg) |
| `rendimiento` | Extract yield (% of total extract relative to sample dry mass) |
| `alt` | Elevational category (factor) |
| `orientación` | Rock face aspect category (factor) |
| `slope` | Rock inclination/slope (°) |
| `exp` | Estimated sun exposure at the microsite (°) |

---

## Statistical Methods

The following analyses were performed in R version 4.0.3:

### Descriptive Statistics
- Median extract yield (%) and usnic acid concentration (%), with corresponding standard deviations, were calculated per microsite condition and elevation (see Table 1 in the paper).

### Generalized Linear Models
1. **Usnic acid concentration ~ Elevation + Rock aspect**
   - Distribution: Gaussian
   - Link function: log
   - Post-hoc: Tukey test to identify which factor levels drove significant differences

2. **Extract yield ~ Elevation + Rock aspect**
   - Distribution: Beta
   - Post-hoc: Pairwise comparisons of estimated marginal means using the `multcomp` and `emmeans` packages

### Model Validation
Normality and homoscedasticity were assessed using diagnostic plots:
- Residuals vs. fitted
- Normal Q-Q
- Scale-location
(see Supplementary Materials S1 and S2; Inchausti 2022)

### Multivariate Analysis
A Principal Component Analysis (PCA) was conducted using the `vegan` package (Oksanen et al. 2022) to visualize relationships among:
- Water retention traits
- Usnic acid concentration
- Environmental conditions

This analysis incorporated data from a previous experiment (Díaz et al. 2022).

---

## References

- Díaz et al. (2022)
- Inchausti (2022)
- Oksanen et al. (2022) – `vegan` package

## Software Requirements

- R version 4.0.3 or compatible
- Required R packages: `vegan`, `multcomp`, `emmeans`

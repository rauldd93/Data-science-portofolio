# Lichen Community Analysis

This repository contains the data processing, statistical modeling, and visualization code for the paper:

**10.5252/cryptogamie-mycologie2026v47a3**

## Analysis Overview

The script `lichen-community-analysis.R` performs the following analyses:

- Data cleaning and processing
- Generalized Linear Models (GLMs)
- Generalized Linear Mixed Models (GLMMs)
- Multivariate analyses (Canonical Correspondence Analysis, Non-metric Multidimensional Scaling)
- Result visualization and plotting

## Required Data Files

The following files must be present in the working directory:

| File | Description |
|------|-------------|
| `dat_comunidades.csv` | Main dataset with lichen species frequencies and environmental variables |
| `traits.csv` | Species trait data for Community Weighted Mean (CWM) calculations |
| `multivariado.csv` | Species matrix for multivariate analyses |

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

Species-level trait data used to calculate Community Weighted Means (CWM) for each growth form.

| Column | Description |
|--------|-------------|
| `spp` | Lichen species name |
| `mft` | Growth form: `cr` = crustose, `fl` = foliose, `fr` = fruticose |
| `genus` | Genus of the species |

---

### `multivariado.csv`

Community matrix for multivariate analyses, with species as columns and environmental variables included.

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
| `pend` | Rock slope (°) |
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

## Notes

- Species names in the `multivariado.csv` file appear as individual columns from `A_altoandina` to `Xanthoria_sp`.
- Cover estimates for vegetation types (`paj`, `rock`, `arb`, `so`, `mg`, `rast`, `ces`) represent percentages within a 5-meter radius of each sampled rock.

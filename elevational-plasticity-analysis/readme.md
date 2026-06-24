# Transplant Experiment Reveals Phenotypic Plasticity in Saxicolous Lichens Along an Elevational Gradient

This repository contains the data and analysis code required to replicate the results published in:

**https://doi.org/10.1016/j.funeco.2026.101497**

## Study Summary

A transplant experiment was conducted along an elevational gradient in Central Argentina to evaluate phenotypic plasticity in physiological and biochemical traits of two saxicolous lichen species, *Parmotrema reticulatum* and *Usnea amblyoclada*. Thalli were transplanted to different elevations and microsites (differing rock aspect), and the following traits were measured:

- Photosynthetic performance (Fo, Fm, Fv/Fm)
- Hydrophobicity
- Water holding capacity (WHCtot, WHCint, WHCext)
- Usnic acid concentration (UA) — *U. amblyoclada* only
- Extract yield (EY) — *U. amblyoclada* only

---

## Analysis Overview

The script `CRH_transplante.R` performs all statistical analyses described in the paper, including:

- Generalized Linear Mixed Models (GLMMs) for physiological and biochemical traits
- Model validation (homoscedasticity, normality of residuals)
- Incidence Rate Ratio (IRR) calculations with 95% confidence intervals
- Post-hoc tests using general linear hypothesis tests and multiple comparisons
- Hydrophobicity curve modeling using survival analysis
- Result visualization and plotting

---

## Required Data Files

The following files must be present in the working directory:

| File | Description |
|------|-------------|
| `CRH_transplante.csv` | Water holding capacity and sample metadata |
| `hidrofobicidad.csv` | Hydrophobicity (water droplet absorption) time measurements |

---

## Statistical Methods

All analyses were performed using R software (R Core Team 2020).

### Generalized Linear Mixed Models (GLMMs)

Each physiological and biochemical parameter was analyzed individually as a response variable, with experimental treatments as fixed effects. The distribution family and link function for each model were selected based on the nature of the response variable:

| Parameter | Description | Distribution | Link Function | Transformation |
|-----------|-------------|--------------|---------------|----------------|
| **Tini** | Onset time of droplet absorption | Negative binomial | log | None (accounts for over-dispersed count data) |
| **WHCtot** | Total water holding capacity | Gaussian | identity | Log-transformed to meet homoscedasticity |
| **WHCint** | Internal water holding capacity | Gaussian | identity | None |
| **WHCext** | External water holding capacity | Gaussian | identity | None |
| **UA** | Usnic acid concentration | Gaussian | identity | None |
| **EY** | Extract yield (proportion) | Beta | logit | None (appropriate for continuous data bounded 0–1) |

**Note:** Photosynthetic parameters (Fo, Fm, Fv/Fm) were also evaluated but showed no significant differences.

### Random Effects Structure

To account for the nested sampling design and avoid pseudo-replication:

- **Individual rocks** were included as a random effect for all models to account for multiple fragments within the same transplant matrix.

- **For WHC metrics (WHCtot, WHCint, WHCext) and hydrophobicity parameters (Fo, Fm, Fv/Fm, Tini):** Sub-samples (sets of three thalli) were modeled as a nested random effect within rocks.

- **For UA and EY:** No nested random effect was required, as all fragments per rock were merged into a composite sample during processing (Inchausti 2022).

The subsample was treated as the primary unit of replication for trait measurement, nested within the transplantation unit, which was the independent experimental unit for treatments.

### Model Validation

Model assumptions (homoscedasticity, normality of residuals) were checked graphically and were satisfied for all final models presented (Inchausti 2022).

### Incidence Rate Ratios (IRRs)

To facilitate interpretation of treatment effects, Incidence Rate Ratios were calculated by exponentiating the fixed-effect coefficients of statistically significant models:

- For log-transformed responses (WHCtot): IRR yields multiplicative effects on the geometric mean.
- For untransformed responses: IRR provides the ratio of geometric means between treatment groups.

IRRs are reported with 95% confidence intervals derived from the Wald approximation (Halvorson et al. 2022). All treatment comparisons were made relative to transplanted control groups that experienced identical handling, transplantation, and monitoring procedures, controlling for shared transplantation stress.

### Post-hoc Tests

*A posteriori* tests were conducted using general linear hypothesis tests and multiple comparisons for the general mixed models from the `emmeans` package (Russell 2021).

### Hydrophobicity Curves

Hydrophobicity curves were calculated using droplet state as a covariate to visualize the behavior of droplet absorption over time, implemented through the `survival` package (Therneau & Grambsch 2000).

---

## Transplant Treatments

| Treatment Code | Description |
|----------------|-------------|
| **Low control** | Transplants from 900 m a.s.l. to another rock at the same elevation/microsite |
| **Low-to-mid** | Transplants from 900 to 1800 m a.s.l. |
| **Mid-to-high** | Transplants from 1800 to 2700 m a.s.l. |
| **North control** | Transplants between north-facing rocks (1800 m a.s.l.) |
| **South control** | Transplants between south-facing rocks (1800 m a.s.l.) |
| **North-to-south** | Transplants from north- to south-facing rocks (1800 m a.s.l.) |
| **South-to-north** | Transplants from south- to north-facing rocks (1800 m a.s.l.) |

---

## File Descriptions

### `CRH_transplante.csv`

Contains water holding capacity measurements and sample metadata for the transplant experiment.

| Column | Description |
|--------|-------------|
| `muestra` | Unique code indicating rock number and randomized position of the sample |
| `roca` | Rock identification number |
| `tmto` | Transplant treatment (see treatment table above) |
| `MS` | Dry mass of the sample (g) |
| `Pmesc` | Mass of the hydrated sample after removing excess water (g) |
| `area` | Surface area of the sample (cm²) |
| `CRHtot` | Total water holding capacity (g H₂O/mm²) |
| `CRHint` | Internal water holding capacity (g H₂O/mm²) |
| `CRHext` | External water holding capacity (g H₂O/mm²) |
| `CHint` | Internal water content (g) |
| `spp` | Lichen species |

---

### `hidrofobicidad.csv`

Contains hydrophobicity measurements based on water droplet absorption timing.

| Column | Description |
|--------|-------------|
| `roca` | Rock identification number |
| `posicion` | Randomized position of the sample |
| `tcontacto` | Time when the droplet was placed on the thallus (s) |
| `tinicio` | Time when the droplet began to be absorbed (s) |
| `tabs` | Time when the droplet was fully absorbed (s) |
| `spp` | Lichen species |
| `tmt` | Transplant treatment |
| `tabsparc` | Partial absorption time: `tabs` − `tinicio` (duration of absorption) |
| `tabstot` | Total absorption time: `tabs` − `tcontacto` (time from placement to full absorption) |
| `tini` | Onset time: time from droplet placement to change in droplet form |
| `status_final` | Droplet status at the time cutoff (90 s after the last droplet was placed) |
| `status_inicial` | Droplet status at the time the droplet was placed |

---

## Software and Packages

| Software/Package | Purpose | Reference |
|------------------|---------|-----------|
| R | Statistical computing | R Core Team (2020) |
| `emmeans` | Post-hoc comparisons and marginal means | Russell (2021) |
| `survival` | Hydrophobicity curve modeling | Therneau & Grambsch (2000) |

---

## References

- R Core Team (2020). *R: A Language and Environment for Statistical Computing.* R Foundation for Statistical Computing, Vienna, Austria.
- Russell, L. (2021). *emmeans: Estimated Marginal Means, aka Least-Squares Means.* R package.
- Therneau, T.M. & Grambsch, P.M. (2000). *Modeling Survival Data: Extending the Cox Model.* Springer, New York.

# Un Experimento de Trasplante Revela Plasticidad Fenotípica en Líquenes Saxícolas a lo Largo de un Gradiente Altitudinal

Este repositorio contiene los datos y el código de análisis necesarios para replicar los resultados publicados en:

**https://doi.org/10.1016/j.funeco.2026.101497**

## Resumen del Estudio

Se realizó un experimento de trasplante a lo largo de un gradiente altitudinal en el centro de Argentina para evaluar la plasticidad fenotípica en rasgos fisiológicos y bioquímicos de dos especies de líquenes saxícolas, *Parmotrema reticulatum* y *Usnea amblyoclada*. Se trasplantaron talos a diferentes altitudes y micrositios (con distinta orientación de la roca), y se midieron los siguientes rasgos:

- Rendimiento fotosintético (Fo, Fm, Fv/Fm)
- Hidrofobicidad
- Capacidad de retención hídrica (CRHtot, CRHint, CRHext)
- Concentración de ácido úsnico (AU) — solo *U. amblyoclada*
- Rendimiento del extracto (RE) — solo *U. amblyoclada*

---

## Resumen del Análisis

El script `CRH_transplante.R` realiza todos los análisis estadísticos descritos en el artículo, incluyendo:

- Modelos Lineales Generalizados Mixtos (GLMM) para rasgos fisiológicos y bioquímicos
- Validación de modelos (homocedasticidad, normalidad de residuos)
- Cálculo de Razones de Tasa de Incidencia (IRR) con intervalos de confianza del 95%
- Pruebas post-hoc mediante pruebas de hipótesis lineales generales y comparaciones múltiples
- Modelado de curvas de hidrofobicidad mediante análisis de supervivencia
- Visualización de resultados y gráficos

---

## Archivos de Datos Requeridos

Los siguientes archivos deben estar presentes en el directorio de trabajo:

| Archivo | Descripción |
|------|-------------|
| `CRH_transplante.csv` | Capacidad de retención hídrica y metadatos de las muestras |
| `hidrofobicidad.csv` | Mediciones de tiempo de hidrofobicidad (absorción de gotas de agua) |

---

## Métodos Estadísticos

Todos los análisis se realizaron con el software R (R Core Team 2020).

### Modelos Lineales Generalizados Mixtos (GLMM)

Cada parámetro fisiológico y bioquímico se analizó individualmente como variable respuesta, con los tratamientos experimentales como efectos fijos. La familia de distribución y la función de enlace para cada modelo se seleccionaron según la naturaleza de la variable respuesta:

| Parámetro | Descripción | Distribución | Función de Enlace | Transformación |
|-----------|-------------|--------------|-------------------|----------------|
| **Tini** | Tiempo de inicio de absorción de la gota | Binomial negativa | log | Ninguna (considera datos de conteo con sobredispersión) |
| **CRHtot** | Capacidad de retención hídrica total | Gaussiana | identidad | Transformación logarítmica para cumplir homocedasticidad |
| **CRHint** | Capacidad de retención hídrica interna | Gaussiana | identidad | Ninguna |
| **CRHext** | Capacidad de retención hídrica externa | Gaussiana | identidad | Ninguna |
| **AU** | Concentración de ácido úsnico | Gaussiana | identidad | Ninguna |
| **RE** | Rendimiento del extracto (proporción) | Beta | logit | Ninguna (apropiada para datos continuos acotados entre 0–1) |

**Nota:** Los parámetros fotosintéticos (Fo, Fm, Fv/Fm) también fueron evaluados, pero no mostraron diferencias significativas.

### Estructura de Efectos Aleatorios

Para considerar el diseño de muestreo anidado y evitar la pseudorreplicación:

- **Las rocas individuales** se incluyeron como efecto aleatorio en todos los modelos para considerar múltiples fragmentos dentro de la misma matriz de trasplante.

- **Para las métricas de CRH (CRHtot, CRHint, CRHext) y los parámetros de hidrofobicidad (Fo, Fm, Fv/Fm, Tini):** Las submuestras (conjuntos de tres talos) se modelaron como un efecto aleatorio anidado dentro de las rocas.

- **Para AU y RE:** No se requirió un efecto aleatorio anidado, ya que todos los fragmentos por roca se combinaron en una muestra compuesta durante el procesamiento (Inchausti 2022).

La submuestra se trató como la unidad primaria de replicación para la medición de rasgos, anidada dentro de la unidad de trasplante, que fue la unidad experimental independiente para los tratamientos.

### Validación del Modelo

Los supuestos del modelo (homocedasticidad, normalidad de los residuos) se verificaron gráficamente y se cumplieron para todos los modelos finales presentados (Inchausti 2022).

### Razones de Tasa de Incidencia (IRR)

Para facilitar la interpretación de los efectos del tratamiento, se calcularon las Razones de Tasa de Incidencia (Incidence Rate Ratios) exponenciando los coeficientes de efectos fijos de los modelos estadísticamente significativos:

- Para respuestas transformadas logarítmicamente (CRHtot): El IRR proporciona efectos multiplicativos sobre la media geométrica.
- Para respuestas no transformadas: El IRR proporciona la razón de medias geométricas entre los grupos de tratamiento.

Los IRR se reportan con intervalos de confianza del 95% derivados de la aproximación de Wald (Halvorson et al. 2022). Todas las comparaciones de tratamientos se realizaron en relación con los grupos de control trasplantados que experimentaron procedimientos idénticos de manipulación, trasplante y monitoreo, controlando así el estrés compartido del trasplante.

### Pruebas Post-hoc

Se realizaron pruebas a posteriori mediante pruebas de hipótesis lineales generales y comparaciones múltiples para los modelos mixtos generales del paquete `emmeans` (Russell 2021).

### Curvas de Hidrofobicidad

Se calcularon curvas de hidrofobicidad utilizando el estado de la gota como covariante para visualizar el comportamiento de la absorción de la gota a lo largo del tiempo, implementado a través del paquete `survival` (Therneau & Grambsch 2000).

---

## Tratamientos de Trasplante

| Código de Tratamiento | Descripción |
|----------------|-------------|
| **Control bajo** | Trasplantes desde 900 m s.n.m. a otra roca en la misma altitud/micrositio |
| **Bajo-a-medio** | Trasplantes desde 900 a 1800 m s.n.m. |
| **Medio-a-alto** | Trasplantes desde 1800 a 2700 m s.n.m. |
| **Control norte** | Trasplantes entre rocas con orientación norte (1800 m s.n.m.) |
| **Control sur** | Trasplantes entre rocas con orientación sur (1800 m s.n.m.) |
| **Norte-a-sur** | Trasplantes desde rocas con orientación norte a sur (1800 m s.n.m.) |
| **Sur-a-norte** | Trasplantes desde rocas con orientación sur a norte (1800 m s.n.m.) |

---

## Descripciones de Archivos

### `CRH_transplante.csv`

Contiene las mediciones de capacidad de retención hídrica y los metadatos de las muestras para el experimento de trasplante.

| Columna | Descripción |
|--------|-------------|
| `muestra` | Código único que indica el número de roca y la posición aleatorizada de la muestra |
| `roca` | Número de identificación de la roca |
| `tmto` | Tratamiento de trasplante (ver tabla de tratamientos arriba) |
| `MS` | Masa seca de la muestra (g) |
| `Pmesc` | Masa de la muestra hidratada tras eliminar el exceso de agua (g) |
| `area` | Área superficial de la muestra (cm²) |
| `CRHtot` | Capacidad de retención hídrica total (g H₂O/mm²) |
| `CRHint` | Capacidad de retención hídrica interna (g H₂O/mm²) |
| `CRHext` | Capacidad de retención hídrica externa (g H₂O/mm²) |
| `CHint` | Contenido hídrico interno (g) |
| `spp` | Especie de liquen |

---

### `hidrofobicidad.csv`

Contiene las mediciones de hidrofobicidad basadas en el tiempo de absorción de gotas de agua.

| Columna | Descripción |
|--------|-------------|
| `roca` | Número de identificación de la roca |
| `posicion` | Posición aleatorizada de la muestra |
| `tcontacto` | Tiempo en que la gota se colocó sobre el talo (s) |
| `tinicio` | Tiempo en que la gota comenzó a ser absorbida (s) |
| `tabs` | Tiempo en que la gota fue completamente absorbida (s) |
| `spp` | Especie de liquen |
| `tmt` | Tratamiento de trasplante |
| `tabsparc` | Tiempo de absorción parcial: `tabs` − `tinicio` (duración de la absorción) |
| `tabstot` | Tiempo de absorción total: `tabs` − `tcontacto` (tiempo desde la colocación hasta la absorción completa) |
| `tini` | Tiempo de inicio: tiempo desde la colocación de la gota hasta el cambio en la forma de la gota |
| `status_final` | Estado de la gota en el tiempo de corte (90 s después de colocar la última gota) |
| `status_inicial` | Estado de la gota en el momento de su colocación |

---

## Software y Paquetes

| Software/Paquete | Propósito | Referencia |
|------------------|---------|-----------|
| R | Computación estadística | R Core Team (2020) |
| `emmeans` | Comparaciones post-hoc y medias marginales | Russell (2021) |
| `survival` | Modelado de curvas de hidrofobicidad | Therneau & Grambsch (2000) |

---

## Referencias

- R Core Team (2020). *R: A Language and Environment for Statistical Computing.* R Foundation for Statistical Computing, Viena, Austria.
- Russell, L. (2021). *emmeans: Estimated Marginal Means, aka Least-Squares Means.* Paquete de R.
- Therneau, T.M. & Grambsch, P.M. (2000). *Modeling Survival Data: Extending the Cox Model.* Springer, Nueva York.

# Transplant Experiment Reveals Phenotypic Plasticity in Saxicolous Lichens Along an Elevational Gradient

This repository contains the data and analysis code required to replicate the results published in:

**https://doi.org/10.1016/j.funeco.2026.101497**

---

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


This repository archives the analyses included in the following paper:

> April Hedd, Katharine R. Studholme, Laura A. McFarlane Tranquilla, Joshua T. Cunningham, Gregory J. Robertson, Sydney M. Collins, David A. Fifield, Carina Gjerdrum, Rielle Hoeg, Patricia L. Jones, William A. Montevecchi, Ingrid L. Pollet, Isabeau Pratte, Jennifer C. Rock, Robert A. Ronconi, Dave Shutler, Sabina I. Wilhelm, Sarah N.P. Wong and Allison Patterson (2026). A species distribution model integrating at-sea behaviour, phenology, and population size to inform marine risk exposure for a threatened seabird during breeding. *Front. Mar. Sci.* 13:1950203.

## Interactive results

Explore the model predictions in the [**[Leach's storm-petrel at-sea species distribution model predictions](https://allison-patterson-lhsp-sdm-viz.share.connect.posit.cloud/)**].
Select one or more colonies, breeding stages, and behavioural risk categories, plus the month and time of day, to view rasters of the number of birds at sea.

## Contents

| Folder | What's inside |
|---|---|
| [`R/`](R/) | All scripts required for the analysis |
| [`raw_data/`](raw_data/) |  Foundation files required for the analyses (`.rds`, `.csv`, `.xlsx`, `.shp`) |
| [`input_data/`](input_data/) | Select interim files required to run the Hidden Markov Model (HMM) and platform proximity scripts, and to initiate the Species Distribution Model (SDM) analysis (`.rds`) |

## Analysis workflow

```mermaid
flowchart LR
    A[Hidden Markov Model] --> B[Species Distribution Model]
    B --> C[Density of birds at sea]
    C --> D[Risk exposure]
```

---

## Folder `R/`

Scripts are listed in the order they are run.

### Setup

Before using the scripts contained in this repository, you must open `template.config.yml`, replace the temporary filepaths and JDK version number with your own, then save the revised file as `config.yml`.

| File | Description |
|---|---|
| `functions.R` | Restore project environment (`renv`), set options, load commonly used values, create directories, and load custom helper functions sourced by the other scripts. |

### Step 1: Movement data and behavioural states (HMM)

| File | Description |
|---|---|
| `00_interpolate_trips_2h.Rmd` | Interpolates GPS trips to a 2-hour interval. Adds covariates for time of day.<br>*Note: for the data cleaning and trip identification steps that precede this one, see `R/_archived_rawdata_clean_tripsplit.Rmd`.* |
| `01_hidden_markov_model.Rmd` | Fits a 3-state Hidden Markov Model to classify at-sea behaviour. Compares models with all combinations of time of day, region, and breeding stage as predictors of the transition probabilities between states (no interactions). Selects best model by AIC, then generates 1) classified positions to feed into the SDM and 2) stationary state probabilities to use with breeding phenology data when estimating the number of birds at sea in each behavioural state. |

### Step 2: Species Distribution Model (SDM)

| File | Description |
|---|---|
| `02_pseudo_abs.Rmd` | Formats colony and tracking data, defines our study area, and gets static predictor layers. Assigns training:validation in our dataset at 80:20 within each colony and breeding stage. Observed locations are resampled at a constant rate within colony and breeding stage. Pseudo absences are then sampled 1:1 within 150% max foraging distance of each source colony. |
| `03_get_env_layers.Rmd` | Generates competition layers, downloads environmental predictor data sets, and extracts all static and dynamic predictor values for our observed and pseudo-absence data. |
| `04_model_fitting.Rmd` | Selects which predictors to model, identifies the best hyper-parameters to use based on leave-one-out cross-validation (where each colony is left out in turn), and fits the final Bayesian Additive Regression Trees (BART) model to our training data. |
| `05_extract_to_predict.Rmd` | Obtains monthly environmental data from the Copernicus Marine Service oceanographic models, calculates derived values as needed and creates standardized rasters for each prediction month, where each raster contains a layer for each predictor on the same grid. |
| `06_prediction.Rmd` | Uses the BART model to make monthly predictions (probability of presence) for each colony in four conditions: 1) on-water behavioural state during the day, 2) on-water behavioural state at night, 3) in-air behavioural state during the day, and 4) in-air behavioural state at night. |

### Step 3: Density of birds at sea

| File | Description |
|---|---|
| `07_nBirds_table.Rmd` | Estimates the number of birds at sea each month for each colony, time of day, breeding stage, and behavioural risk category. Uses hatch date estimates, stationary state probabilities from the HMM, and breeding population size estimates. |
| `08_nBirds_rast.Rmd` | Convert probability of presence rasters into rasters with the number of birds per pixel. |

### Step 4: Risk exposure

| File | Description |
|---|---|
| `09_prep_risk_polys.Rmd` | Generates 1) a shapefile of current offshore petroleum lease areas and 2) a shapefile including only production and exploration licenses for use in the risk exposure case studies. |
| `10_platform_proximity.Rmd` | Calculates for each foraging trip the single closest passage to an offshore petroleum production installation and all passages within any light or oil spill footprint of these facilities. Determines whether these passages occurred under conditions that could result in exposure to light or oil. |

### Step 5: Outputs

| File | Description |
|---|---|
| `11_results_summaries.Rmd` | Generates the summaries used in the paper. |
| `12_manuscript_figures.Rmd` | Produces the figures included in the paper. |

### Archived

| File | Description |
|---|---|
| `_archived_rawdata_clean_tripsplit.Rmd` | Loads and cleans raw GPS files (`.pos`), defines trips, and summarizes trip metrics. The raw GPS files required to run this script are not provided - the script is retained for methodological transparency only. The outputs of this script that are necessary to run all subsequent steps of our analyses are included. |

---

## Folder `raw_data/`

| File / folder | Description |
|---|---|
| `data_tripSplit_df.rds` | The cleaned and prepared GPS spatial dataset, with foraging trips identified |
| `tripSummaryLHSP.rds` | Summary of trip characteristics, including duration, maximum distance, trip type, and breeding stage |
| `MasterMetadata_LHSP_Tracking_NL_NB_NS_13Aug2024_GPSonly.csv` | Metadata (deployment, retrieval) for all GPS-tracked birds |
| `colony.csv` | Colony locations from which birds were tracked |
| `LHSP Atlantic population survey data 10Sep2026.xlsx` | Atlantic Leach's storm-petrel population sizes |
| `LHSP_allHatchEstimates.csv` | Hatch date estimates by colony, based on direct observation or chick measurements |
| `EEZ_Canada/200_mile.shp` | Canadian 200-nautical-mile EEZ boundary |
| `threats/nl_petroleum` | Newfoundland and Labrador petroleum licence layers (folder of `.shp`) |
| `threats/ns_petroleum` | Nova Scotia petroleum licence layers (folder of `.shp`, `.xls`) |
| `threats/threats_offshorePetro_productionInstallationPts.csv` | Locations of offshore petroleum production installations |

---

## Folder `input_data/`

This folder contains select interim files required to run the Hidden Markov Model (HMM) and platform proximity scripts, and to initiate the Species Distribution Model (SDM) analysis. These are output files from the GPS trip interpolation and Hidden Markov Modelling scripts.

| File | Used by | Description |
|---|---|---|
| `data_tripSplit_2h.rds` | `01_hidden_markov_model.Rmd` | Two-hour interpolated GPS trips |
| `HMM_classified_data.rds` | `02_pseudo_abs.Rmd`<br>`10_platform_proximity.Rmd` | Locations classified by the 3-state HMM |
| `HMM_stationary_states.rds` | `07_nBirds_table.Rmd` | Stationary-state probabilities from the 3-state HMM fit |

---

## Software requirements

- R version 4.3.1
- Key packages: Managed with `renv.lock`

## Citation

If you use this code or data, please cite the paper above, as well as the Zenodo repository and any relevant Movebank data packages. 

Movebank data packages include:

**Baccalieu Island:**
Hedd A, Cunningham JT, McFarlane Tranquilla LA, Robertson GJ, Studholme KR, Wilhelm SI. 2026. Data from: Study "Leach's storm-petrel GPS - Baccalieu Is. - Atlantic Canada - 2016-2023". Movebank Data Repository. <https://doi.org/10.5441/001/1.757>

**Bon Portage Island:**
Hoeg R, Pollet I, Shutler D, Hedd A. 2026. Data from: Study "Leach's storm-petrel GPS - Bon Portage Is. - Atlantic Canada - 2017-2022". Movebank Data Repository. <https://doi.org/10.5441/001/1.759>

**Country Island:**
Rock JC, Gjerdrum C, Ronconi RA, Wilhelm SI, Wong SNP, Hedd A. 2026. Data from: Study "Leach's storm-petrel GPS - Country Is. - Atlantic Canada - 2017-2022". Movebank Data Repository. <https://doi.org/10.5441/001/1.760>

**Gull Island**
Hedd A, Collins SM, Cunningham JT, Montevecchi WA, Robertson GJ, Wilhelm SI. 2026. Data from: Study "Leach's storm-petrel GPS - Gull Is. (Witless Bay) - Atlantic Canada - 2016-2023". Movebank Data Repository. <https://doi.org/10.5441/001/1.761>

**Kent Island:**
Jones PL, Ronconi RA, McFarlane Tranquilla LA, Wilhelm SI, Wong SNP, Hedd A. 2026. Data from: Study "Leach's storm-petrel GPS - Kent Is. - Atlantic Canada - 2018-2022". Movebank Data Repository. <https://doi.org/10.5441/001/1.762>

**Middle Lawn Island**
Studholme KR, Cunningham JT, Hedd A. 2026. Data from: Study "Leach's storm-petrel GPS - Middle Lawn Is. - Atlantic Canada - 2022-2023". Movebank Data Repository. <https://doi.org/10.5441/001/1.763>

## Contact

Allison Patterson <allison.patterson@ec.gc.ca>

Katharine Studholme <katharine.studholme@ec.gc.ca>

## License

[CCBY](https://creativecommons.org/licenses/by/4.0/)

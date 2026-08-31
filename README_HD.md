# Failure to invest below-ground may limit the northern expansion of invasive knotweed: lessons from a two-phase transplant experiment

This repository stores the data and scripts related to the manuscript 
*Failure to invest below-ground may limit the northern expansion of invasive knotweed: lessons from a two-phase transplant experiment*, 
including the raw biomass data and the R script used for data cleaning, 
modelling and figures.

## Contents

The following materials are available in this repository:

- [`CombinedBiomass_250717.csv`](CombinedBiomass_250717.csv) - Biomass data from the second-phase (T2) re-transplant experiment, 
including population of origin, latitude of origin, garden of origin (2022) and re-transplant garden (2023), block, 
initial rhizome weight, sprouting and survival records, and dry biomass (above-ground, below-ground, total) for each individual.
- [`Knotweed_260821.R`](Knotweed_260821.R) - R script used to prepare and analyse the biomass data.
- [`Reynoutria_SDM.R`](Reynoutria_SDM.R) - R script for the species distribution model (SDM), 
producing the 4-panel figure: (A) PCA biplot of the bioclimatic variables, (B) PCA climate space comparing the European background with European and Japanese (octoploid) occurrences, (C) a 6-class MaxEnt map of predicted presence probability across Europe 
(occurrences cleaned from GBIF/iNaturalist records, background points drawn with `megaSDM`, model tuned with `ENMeval` and fitted with `dismo::maxent`), and (D) distributions and boxplots of the main predictors, mean annual temperature and annual temperature range, in grid cells with predicted probability of occurrence > 0.5, with climatic limits indicated as 2.5th and 97.5th quantiles.
- [`PloidyLevel_R_japo_Japan.xlsx`](PloidyLevel_R_japo_Japan.xlsx) - Ploidy level records (4x/8x) 
for *Reynoutria japonica* populations in Japan, with collection locality, prefecture, coordinates 
and literature source, used to restrict the Japanese occurrence data in the SDM to the octoploid 
cytotype comparable to the invasive European populations. Additional details are provided in the manuscript.
- [`Reynoutria_japonica_occurrences_00_25.rar`](Reynoutria_japonica_occurrences_00_25.rar) - *Reynoutria japonica* occurrences used for the SDM.
- CHELSA bioclimatic rasters used in the analyses.


## Abstract

The ecological and evolutionary processes determining species limits remain poorly understood. Ultimately, range 
limits depend on the species' abilities to persist under heterogeneous conditions 
by adaptive differentiation and phenotypic plasticity, including transgenerational effects. 
To investigate ecological differentiation and transgenerational effects in the clonal invasive 
knotweed, *Reynoutria japonica*, in Europe, we carried out a two-phase transplant experiment: 
plants sampled along a latitudinal gradient were planted at three sites located at the northern range margin, 
mid-range and near the southern range margin, and then re-transplanted among all three sites after two years. 
Biomass production and allocation were generally not associated with latitude of origin, 
and previous growth at the same site did not promote performance. We therefore found no evidence that 
adaptive differentiation or adaptive transgenerational effects contribute to the wide distribution of 
*R. japonica* in Europe. However, at the northern site, with a 25% shorter season, 
knotweed plants invested much less biomass below-ground, 
even though overwintering below-ground rhizomes are essential for the species' survival and spread. 
This pattern was further strengthened in plants that had grown at the northern site in the previous generation. 
We further explored limiting climate conditions in a species distribution model for the European range and 
found that mean annual temperature and temperature annual range are the main predictors of the European distribution of 
*R. japonica*. Taken together, our study suggests that low temperatures and associated short growing seasons may 
pose a limit to the broad environmental tolerance of *R. japonica* and restrict its northward spread 
by reducing below-ground biomass accumulation.

**Keywords:** below-ground investment, clonal plants, invasive species, reciprocal transplant experiment, 
species distribution model, species range limits, transgenerational effects

## Acknowledgments
This study was supported by the German Federal Ministry of Education and Research (BMBF; MOPGA Project 306055 to 
Christina L. Richards), the German Research Foundation (DFG; grant 431595342 to Oliver Bossdorf and Christina L. Richards) 
the European Union's Horizon 2020 research and innovation program under the Marie Skłodowska-Curie grant agreement No 101033168 (to Ramona E. Irimia), 
and the MSCA4Ukraine program of the Alexander von Humboldt Foundation (doctoral fellowship No. 1233613 to Hanna Danko).


## Citation

Please cite the repository, dataset and article as (to complete once the manuscript has been accepted for publication).



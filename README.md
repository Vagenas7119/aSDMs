# aquatic Species Distribution Model (aSDMs)

The files included herein represent the complete set of outputs and the scripts required for the analysis.

![SIBECOL_Vagenas_scroped](https://github.com/user-attachments/assets/2a56da26-b665-4091-b802-3541cb3330a8)

```
# Project structure
├── data/                     # Local data files (bioclimatic variables, hydrological baseline layers etc.)
├── figures/                  # Output plots and figures
├── supplementary/            # Output supplementary plots and figures
├── outputs/                  # Output intermediate figures 
└── Scripts.R                 # Primary analysis script (R Markdown)
```

## Abstract:
Species Distribution Models (SDMs) have traditionally been developed for terrestrial and open-land systems, yet their application to aquatic ecosystems presents distinct challenges. These include predicting and projecting species distributions across spatial restricted structures, incorporating key environmental drives such as climate and water availability and the lack of a standardized modelling framework to address these specificities. We explore methodological solutions enabling accurate high-resolution SDMs for freshwater organisms using presence-only records. Specifically, we focused on the Ichthyofauna of the Iberian Peninsula and we applied models trained globally and projected regionally (i.e., global-to-regional models) for non-endemic species and models trained and projected regionally (i.e., strictly regional models) for endemic species. Our study systematically compares two modelling approaches for ensemble aquatic SDMs: unconstrained models trained across the entire freshwater species ranges, and constrained models, in which model training is limited to watersheds where species occur. Additionally, we tested different combinations of environmental variables, from individual predictors to  hierarchical combinations of climatic and hydrological predictors. The results demonstrate that the pre-constrained spatial modelling strategy is more effective. Moreover, models trained with climate predictors consistently outperformed models trained with hydrological covariates alone. We conclude that all proposed modelling stages should be followed to accurately predict aquatic species distributions. This ensures comprehensive representation of spatial variability, appropriate covariate selection, and optimal model configuration—accounting for the unique complexity of aquatic ecosystems.

### Keywords: 
SDMs, freshwater ecosystems, fish, climate, hydrology

#### Citation (APA):
 Vagenas, G., Matias, M., Araujo M.B. (2025). Beyond land: a framework for modelling aquatic species distributions.
 
#### DOI:  [Pending]

##### Figures - Main Manuscript

![Figure1_cropped_jpg](https://github.com/user-attachments/assets/b655c959-c405-4b5d-a443-8130f5ef321d)
Figure 1. Spatial distribution of species richness of the dataset used for the development of the aSDMs through a (A) global (GBIF | 50 arc-minute grid) to (B) regional (MITECO, SNIPAD, GBIF | 10 arc-minute grid) for the endemics and widespread freshwater fish species of the study area. Colors represent gradients of species richness (low = yellow; high = red). The finer resolution (B) highlights richness patterns in the Iberian Peninsula.

![Figure_1_v13](https://github.com/user-attachments/assets/1c169e20-9673-494c-a66c-65a2674d9639)
Figure 2. Flowchart illustrating the implementation and evaluation workflow for aquatic Species Distribution Models (aSDMs), comprising nine sequential stages (i.e., I-IX), from input data preparation and modelling through to performance evaluation and the generation of stacked suitability maps.

![Figure3_jpg](https://github.com/user-attachments/assets/fdedcc87-2717-442b-8759-81a66303ade9)
Figure 3. Classification flowchart based on recursive partitioning. Each tree uses aSDM performance metrics (AUC, CBI, TSS) as the classification target. Iterative branch splitting identifies high-performing model classes (highlighted in green). The label alien and clima+ indicate the factors of non-endemic/introduced species and all the variables incorporating climate excluding hydrological covariates, respectively.

![Figure4_jpg_cropped](https://github.com/user-attachments/assets/eee40ddf-cef8-4e8b-8110-b4582e1864c1)
Figure 4. Stacked aSDMs for endemic (N=39) and non-endemic/introduced species (N=53) across the study area. The maps represent stacked outputs derived through aSDMs using the dominant pre-constrained h5 spatial strategy, based on climate (left) and hydrological (right) covariate sets. Boxes indicate the high-resolution (~1x1 km) aSDMs projected across the hydrographic network. Highlighted zoomed-in areas are illustrative examples for visual comparison and do not represent specific ecological patterns.

#Author: Georgios Vagenas
Affiliation: Biogeography and Global Change Department, National Museum of Natural Sciences, CSIC, C/ Jose Gutierrez Abascal, 2, Madrid 28006, Spain





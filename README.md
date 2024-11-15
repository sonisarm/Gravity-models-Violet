### Code for the following publication: "*Connectivity in two endangered violets in the Canary Islands: how topography, temperature and water availability shape gene flow and diversity in Viola cheiranthifolia and Viola guaxarensis*"

This GitHub page provides a step-by-step guide to the analyses presented in the publication, "*Connectivity in two endangered violets in the Canary Islands: how topography, temperature, and water availability shape gene flow and diversity in Viola cheiranthifolia and Viola guaxarensis.*" This study is part of the GENCLIMA project, which aims to evaluate functional connectivity among high-mountain and laurel forest species in the Canary Islands, with a focus on preserving genetic diversity and mitigating the impacts of climate change.

The project's specific objectives are to (1) identify correlations between genetic diversity and environmental variables, (2) assess the functional connectivity of V. cheiranthifolia and V. guaxarensis at El Teide Summit, Tenerife, Spain, (3) determine which landscape variables influence gene flow in each species, and (4) evaluate changes in connectivity under projected climate change scenarios.

Remark: The input data for the scripts in the following sections can be found in folder *00-data*.


![framework](https://github.com/sonisarm/grav-models-violet/blob/main/framework-connectivity-analysis.jpg)
### 00-Sampling
The sampling scheme used here is individual-based, differing from the population-level analyses commonly applied in gravity models for metapopulations. This approach enhances our ability to detect complex interactions of environmental variables on gene flow patterns and allows for a more precise identification of landscape factors independently influencing gene flow. 


### 01-Genetic diversity
This script evaluates the relationship between genetic diversity and environmental variables for Teide Violets. It processes spatial and genetic data, extracts environmental variables from raster files, and calculates population-level means. Linear models are applied to test the influence of each variable on genetic diversity, identifying significant predictors. Collinearity between variables is assessed, and results are visualized with annotated plots to support ecological interpretations.



### 02-Gravity models
This section provides a comprehensive script for building and evaluating gravity models to study genetic flow between individuals. The code integrates spatial and genetic data, extracts relevant environmental variables, assesses collinearity, and constructs gravity models to assess functional connectivity and identify key landscape factors influencing genetic connectivity. Additionally, it includes tools for comparing and refining models.



### 03-Climate modelling
This script evaluates the impact of climate change on connectivity and gene flow across future scenarios. It integrates climate data, gene flow values, and site coordinates, processing them to derive environmental statistics for nodes and edges while incorporating topographic variables. Gravity models are applied to analyze connectivity, with comparisons between current and future scenarios to assess changes and inform potential conservation strategies.











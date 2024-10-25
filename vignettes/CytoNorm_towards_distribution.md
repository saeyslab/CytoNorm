Use case 2: Normalize towards specific distribution
================

<!-- github markdown built using
rmarkdown::render("vignettes/CytoNorm_towards_distribution.Rmd", output_format = "github_document")
-->

#### In this vignette, you learn how CytoNorm can be used to normalize a dataset towards a specific distribution.

Often, data is analyzed over a longer period of time, especially in
clinical trials. Instead of waiting for all data to start the
computational data analysis, one can already analyze a first part of the
data and integrate the later data afterwards. A pitfall here is that the
later data should not be normalized towards a median over the batches
because otherwise the earlier data analysis, clustering for example,
needs to be re-done. Rather, we would like to normalize towards the goal
distribution of a previous model. In another experiment, it might be the
case that manual gating is done on a single batch. In order to apply
this gating strategy to a second batch, the second batch needs to be
normalized towards the first batch. It is clear that some experimental
settings or designs require the normalization algorithm to be flexible,
specifically regarding the normalization goal. In most normalization
algorithms, this is not the case.

By default, CytoNorm normalizes the data towards the median over the
batches, but it is also convenient to specify the goal distribution. It
is possible to assign one of the batches as “goal batch”, or a goal
distribution from a previous CytoNorm model can be retrieved and passed
along in the call when training a new CytoNorm model.

## Prepare the CytoNorm normalization step

In the [Use case 1: Run CytoNorm on a dataset without dedicated controls
vignette](CytoNorm_without_controls.md),the two first cohort batches of
our dataset were normalized towards each other.Now, we want to normalize
the data of the second cohort (P1_C2 and P2_C2) towards the level of the
first.

To do this, we extract the CytoNorm goal quantiles of the model that
wasgenerated for the [first use case](CytoNorm_without_controls.md) and
pass these along while training the CytoNorm model over here. This
CytoNorm model will use the same FlowSOM model as the first model and
will normalize the same (backbone) markers.

``` r
library(CytoNorm)
library(flowCore)
library(FlowSOM)
library(ggpubr)
```

``` r
dir_prepr <- "Data/Preprocessed"

# Organize preprocessed data in batches
files_P1_C1 <- list.files(path = dir_prepr,
                          pattern = "ID[1-4]_Panel1_TP[1-3].fcs",
                          full.names = TRUE)
files_P1_C2 <- list.files(path = dir_prepr,
                          pattern = "ID[5-8]_Panel1_TP[1-3].fcs",
                          full.names = TRUE)
files_P2_C1 <- list.files(path = dir_prepr,
                          pattern = "ID[1-4]_Panel2_TP[1-3].fcs",
                          full.names = TRUE)
files_P2_C2 <- list.files(path = dir_prepr,
                          pattern = "ID[5-8]_Panel2_TP[1-3].fcs",
                          full.names = TRUE)

# Read in 1 file per panel to get channel information
ff_P1 <- read.FCS(files_P1_C1[1])
ff_P2 <- read.FCS(files_P2_C1[1])

# Retrieve channel information
cellTypeMarkers <- c("CD27", "CD38", "IgD", "IgM") 
cellTypeChannels <- GetChannels(object = ff_P1, markers = cellTypeMarkers)
P1_cellStateMarkers <- c("CD5", "CD268", "CD24", "PDL1", "PD1", "CD86", "KI67")
P1_cellStateChannels <- GetChannels(object = ff_P1, markers = P1_cellStateMarkers)
P2_cellStateMarkers <- c("CD40", "CD21", "HLA-DR", "ICOSL", "PD1", "CD86", "KI67")
P2_cellStateChannels <- GetChannels(object = ff_P2, markers = P2_cellStateMarkers)

# Read in manualLabels object that was generated in the CytoNorm_without_controls vignette
manualLabels <- readRDS("Data/Preprocessed/attachments/manualLabels.rds")

# Normalized data
files_P1_C1_norm <- list.files(path = "Data/Normalized_cellType", 
                               pattern = "ID[1-4]_Panel1_TP..fcs", full.names = TRUE)
files_P2_C1_norm <- list.files(path = "Data/Normalized_cellType", 
                               pattern = "ID[1-4]_Panel2_TP..fcs", full.names = TRUE)
```

## Train the CytoNorm model

### Extract the goal distribution

``` r
previousModel <- readRDS("RDS/model1_withoutControl.rds")
goal_q <- getCytoNormQuantiles(previousModel)
```

### Train the CytoNorm model

``` r
model2_towardsDistribution <- CytoNorm.train(files = c(files_P1_C2, files_P2_C2),
                                             labels = c(rep(x = "P1",times = length(files_P1_C2)),
                                                        rep(x = "P2",times = length(files_P2_C2))),
                                             channels = cellTypeChannels,
                                             transformList = NULL,
                                             seed = 1,
                                             plot = TRUE,
                                             verbose = TRUE,
                                             normParams = list("goal" = goal_q))
#> Warning in CytoNorm.train(files = c(files_P1_C2, files_P2_C2), labels = c(rep(x = "P1", : Reusing FlowSOM result previously saved at ./tmp/CytoNorm_FlowSOM.RDS. 
#>          If this was not intended, one can either specify another outputDir, 
#>          make use of the recompute parameter or move the FlowSOM object in the 
#>          file manager.
#> Splitting Data/Preprocessed/ID5_Panel1_TP1.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 967 cells (7.8%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID5_Panel1_TP2.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1867 cells (8.19%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID5_Panel1_TP3.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 2332 cells (7.53%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID6_Panel1_TP1.fcs
#> Splitting Data/Preprocessed/ID6_Panel1_TP2.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1153 cells (5.36%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID6_Panel1_TP3.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1167 cells (5.71%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID7_Panel1_TP1.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 2216 cells (11.97%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID7_Panel1_TP2.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1808 cells (9.71%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID7_Panel1_TP3.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1519 cells (10.23%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID8_Panel1_TP1.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 4863 cells (6.92%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID8_Panel1_TP2.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 4936 cells (8.58%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID8_Panel1_TP3.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 5572 cells (8.21%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID5_Panel2_TP1.fcs
#> Splitting Data/Preprocessed/ID5_Panel2_TP2.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1784 cells (8.56%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID5_Panel2_TP3.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 2060 cells (6.87%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID6_Panel2_TP1.fcs
#> Splitting Data/Preprocessed/ID6_Panel2_TP2.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1448 cells (6.18%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID6_Panel2_TP3.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1284 cells (6.42%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID7_Panel2_TP1.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 2637 cells (13.36%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID7_Panel2_TP2.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1987 cells (10.4%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID7_Panel2_TP3.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1759 cells (10.37%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID8_Panel2_TP1.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 4814 cells (7.08%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID8_Panel2_TP2.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 5119 cells (7.71%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID8_Panel2_TP3.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 6582 cells (8.87%) seem far from their cluster centers.
#> Processing cluster 1
#> Computing Quantiles
#>   P1 (FileID 1,2,3,4,5,6,7,8,9,10,11,12)
#> Reading ./tmp/Data_Preprocessed_ID5_Panel1_TP1.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID5_Panel1_TP2.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID5_Panel1_TP3.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel1_TP1.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel1_TP2.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel1_TP3.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel1_TP1.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel1_TP2.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel1_TP3.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel1_TP1.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel1_TP2.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel1_TP3.fcs_fsom1.fcs
#>   P2 (FileID 13,14,15,16,17,18,19,20,21,22,23,24)
#> Reading ./tmp/Data_Preprocessed_ID5_Panel2_TP1.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID5_Panel2_TP2.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID5_Panel2_TP3.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel2_TP1.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel2_TP2.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel2_TP3.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel2_TP1.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel2_TP2.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel2_TP3.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel2_TP1.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel2_TP2.fcs_fsom1.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel2_TP3.fcs_fsom1.fcs
#> Computing Splines
#>   P1
#>   P2
#> Processing cluster 2
#> Computing Quantiles
#>   P1 (FileID 1,2,3,4,5,6,7,8,9,10,11,12)
#> Reading ./tmp/Data_Preprocessed_ID5_Panel1_TP1.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID5_Panel1_TP2.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID5_Panel1_TP3.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel1_TP1.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel1_TP2.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel1_TP3.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel1_TP1.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel1_TP2.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel1_TP3.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel1_TP1.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel1_TP2.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel1_TP3.fcs_fsom2.fcs
#>   P2 (FileID 13,14,15,16,17,18,19,20,21,22,23,24)
#> Reading ./tmp/Data_Preprocessed_ID5_Panel2_TP1.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID5_Panel2_TP2.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID5_Panel2_TP3.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel2_TP1.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel2_TP2.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel2_TP3.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel2_TP1.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel2_TP2.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel2_TP3.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel2_TP1.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel2_TP2.fcs_fsom2.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel2_TP3.fcs_fsom2.fcs
#> Computing Splines
#>   P1
#>   P2
#> Processing cluster 3
#> Computing Quantiles
#>   P1 (FileID 1,2,3,4,5,6,7,8,9,10,11,12)
#> Reading ./tmp/Data_Preprocessed_ID5_Panel1_TP1.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID5_Panel1_TP2.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID5_Panel1_TP3.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel1_TP1.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel1_TP2.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel1_TP3.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel1_TP1.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel1_TP2.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel1_TP3.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel1_TP1.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel1_TP2.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel1_TP3.fcs_fsom3.fcs
#>   P2 (FileID 13,14,15,16,17,18,19,20,21,22,23,24)
#> Reading ./tmp/Data_Preprocessed_ID5_Panel2_TP1.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID5_Panel2_TP2.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID5_Panel2_TP3.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel2_TP1.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel2_TP2.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID6_Panel2_TP3.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel2_TP1.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel2_TP2.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID7_Panel2_TP3.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel2_TP1.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel2_TP2.fcs_fsom3.fcs
#> Reading ./tmp/Data_Preprocessed_ID8_Panel2_TP3.fcs_fsom3.fcs
#> Computing Splines
#>   P1
#>   P2
```

``` r
                                             #normParams = list("goal" = "P1")) #E.g. to normalize towards one of the batches
saveRDS(object = model2_towardsDistribution, file = "RDS/model2_towardsDistribution.rds")
```

The `normParams` argument can be used to pass parameters to the
normalization method call, `QuantileNorm.train` by default. Here, you
can specify for example `nQ`, the number of quantiles that will be
computed and aligned, or goal, which are the goal quantiles. Note that
no `FlowSOM.params` are specified, since CytoNorm will look in the
`outputDir` (default = “./tmp”) for a previous FlowSOM object. If this
is the case, it will throw a warning saying CytoNorm will reuse the
saved FlowSOM result. Here, it will thus reuse the FlowSOM clustering
from the previous CytoNorm model ([use case
1](CytoNorm_without_controls.md)). If this is not the desired behavior,
one can either specify another `outputDir` or move the FlowSOM object in
the file manager.

## Apply the CytoNorm model

``` r
CytoNorm.normalize(model = model2_towardsDistribution,
                   files = c(files_P1_C2, files_P2_C2),
                   labels = c(rep(x = "P1",times = length(files_P1_C2)),
                              rep(x = "P2",times = length(files_P2_C2))),
                   transformList = NULL,
                   verbose = TRUE,
                   prefix = "",
                   transformList.reverse = NULL, 
                   outputDir = "Data/Normalized_cellType")
#> Splitting Data/Preprocessed/ID5_Panel1_TP1.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 967 cells (7.8%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID5_Panel1_TP2.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1867 cells (8.19%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID5_Panel1_TP3.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 2332 cells (7.53%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID6_Panel1_TP1.fcs
#> Splitting Data/Preprocessed/ID6_Panel1_TP2.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1153 cells (5.36%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID6_Panel1_TP3.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1167 cells (5.71%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID7_Panel1_TP1.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 2216 cells (11.97%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID7_Panel1_TP2.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1808 cells (9.71%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID7_Panel1_TP3.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1519 cells (10.23%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID8_Panel1_TP1.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 4863 cells (6.92%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID8_Panel1_TP2.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 4936 cells (8.58%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID8_Panel1_TP3.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 5572 cells (8.21%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID5_Panel2_TP1.fcs
#> Splitting Data/Preprocessed/ID5_Panel2_TP2.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1784 cells (8.56%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID5_Panel2_TP3.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 2060 cells (6.87%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID6_Panel2_TP1.fcs
#> Splitting Data/Preprocessed/ID6_Panel2_TP2.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1448 cells (6.18%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID6_Panel2_TP3.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1284 cells (6.42%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID7_Panel2_TP1.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 2637 cells (13.36%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID7_Panel2_TP2.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1987 cells (10.4%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID7_Panel2_TP3.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 1759 cells (10.37%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID8_Panel2_TP1.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 4814 cells (7.08%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID8_Panel2_TP2.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 5119 cells (7.71%) seem far from their cluster centers.
#> Splitting Data/Preprocessed/ID8_Panel2_TP3.fcs
#> Warning in FlowSOM::NewData(fsom, ff): 6582 cells (8.87%) seem far from their cluster centers.
#> Processing cluster 1
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel1_TP1.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel1_TP2.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel1_TP3.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel1_TP1.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel1_TP2.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel1_TP3.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel1_TP1.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel1_TP2.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel1_TP3.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel1_TP1.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel1_TP2.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel1_TP3.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel2_TP1.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel2_TP2.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel2_TP3.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel2_TP1.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel2_TP2.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel2_TP3.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel2_TP1.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel2_TP2.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel2_TP3.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel2_TP1.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel2_TP2.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel2_TP3.fcs_fsom1.fcs (P2)
#> Normalizing P2
#> Processing cluster 2
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel1_TP1.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel1_TP2.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel1_TP3.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel1_TP1.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel1_TP2.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel1_TP3.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel1_TP1.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel1_TP2.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel1_TP3.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel1_TP1.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel1_TP2.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel1_TP3.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel2_TP1.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel2_TP2.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel2_TP3.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel2_TP1.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel2_TP2.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel2_TP3.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel2_TP1.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel2_TP2.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel2_TP3.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel2_TP1.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel2_TP2.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel2_TP3.fcs_fsom2.fcs (P2)
#> Normalizing P2
#> Processing cluster 3
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel1_TP1.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel1_TP2.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel1_TP3.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel1_TP1.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel1_TP2.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel1_TP3.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel1_TP1.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel1_TP2.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel1_TP3.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel1_TP1.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel1_TP2.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel1_TP3.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel2_TP1.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel2_TP2.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID5_Panel2_TP3.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel2_TP1.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel2_TP2.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID6_Panel2_TP3.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel2_TP1.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel2_TP2.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID7_Panel2_TP3.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel2_TP1.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel2_TP2.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID8_Panel2_TP3.fcs_fsom3.fcs (P2)
#> Normalizing P2
#> Rebuilding Data/Preprocessed/ID5_Panel1_TP1.fcs
#> Rebuilding Data/Preprocessed/ID5_Panel1_TP2.fcs
#> Rebuilding Data/Preprocessed/ID5_Panel1_TP3.fcs
#> Rebuilding Data/Preprocessed/ID6_Panel1_TP1.fcs
#> Rebuilding Data/Preprocessed/ID6_Panel1_TP2.fcs
#> Rebuilding Data/Preprocessed/ID6_Panel1_TP3.fcs
#> Rebuilding Data/Preprocessed/ID7_Panel1_TP1.fcs
#> Rebuilding Data/Preprocessed/ID7_Panel1_TP2.fcs
#> Rebuilding Data/Preprocessed/ID7_Panel1_TP3.fcs
#> Rebuilding Data/Preprocessed/ID8_Panel1_TP1.fcs
#> Rebuilding Data/Preprocessed/ID8_Panel1_TP2.fcs
#> Rebuilding Data/Preprocessed/ID8_Panel1_TP3.fcs
#> Rebuilding Data/Preprocessed/ID5_Panel2_TP1.fcs
#> Rebuilding Data/Preprocessed/ID5_Panel2_TP2.fcs
#> Rebuilding Data/Preprocessed/ID5_Panel2_TP3.fcs
#> Rebuilding Data/Preprocessed/ID6_Panel2_TP1.fcs
#> Rebuilding Data/Preprocessed/ID6_Panel2_TP2.fcs
#> Rebuilding Data/Preprocessed/ID6_Panel2_TP3.fcs
#> Rebuilding Data/Preprocessed/ID7_Panel2_TP1.fcs
#> Rebuilding Data/Preprocessed/ID7_Panel2_TP2.fcs
#> Rebuilding Data/Preprocessed/ID7_Panel2_TP3.fcs
#> Rebuilding Data/Preprocessed/ID8_Panel2_TP1.fcs
#> Rebuilding Data/Preprocessed/ID8_Panel2_TP2.fcs
#> Rebuilding Data/Preprocessed/ID8_Panel2_TP3.fcs
```

## Evauate the quality

### Density plots

``` r
# List files
files_P1_C2_norm <- list.files(path = "Data/Normalized_cellType", 
                               pattern = "ID[5-8]_Panel1_TP[1-3].fcs", full.names = TRUE)
files_P2_C2_norm <- list.files(path = "Data/Normalized_cellType", 
                               pattern = "ID[5-8]_Panel2_TP[1-3].fcs", full.names = TRUE)

# Make plots
p <- plotDensities(input = list("P1_C2" = files_P1_C2,
                                "P2_C2" = files_P2_C2,
                                "P1_C2_norm" = files_P1_C2_norm,
                                "P2_C2_norm" = files_P2_C2_norm),
                   channels = cellTypeChannels,
                   colors = c("P1_C1" = "#900002","P1_C2" = "#FF3437",
                              "P2_C1" = "#0046A1","P2_C2" = "#6FCFFF",
                              "C1_goal" = "#7209b7", "Goal distribution" = "#7209b7"),
                   model = model2_towardsDistribution,
                   show_goal = TRUE)
#> Reading Data/Preprocessed/ID5_Panel1_TP1.fcs
#> Reading Data/Preprocessed/ID5_Panel1_TP2.fcs
#> Reading Data/Preprocessed/ID5_Panel1_TP3.fcs
#> Reading Data/Preprocessed/ID6_Panel1_TP1.fcs
#> Reading Data/Preprocessed/ID6_Panel1_TP2.fcs
#> Reading Data/Preprocessed/ID6_Panel1_TP3.fcs
#> Reading Data/Preprocessed/ID7_Panel1_TP1.fcs
#> Reading Data/Preprocessed/ID7_Panel1_TP2.fcs
#> Reading Data/Preprocessed/ID7_Panel1_TP3.fcs
#> Reading Data/Preprocessed/ID8_Panel1_TP1.fcs
#> Reading Data/Preprocessed/ID8_Panel1_TP2.fcs
#> Reading Data/Preprocessed/ID8_Panel1_TP3.fcs
#> Reading Data/Normalized_cellType/ID5_Panel1_TP1.fcs
#> Reading Data/Normalized_cellType/ID5_Panel1_TP2.fcs
#> Reading Data/Normalized_cellType/ID5_Panel1_TP3.fcs
#> Reading Data/Normalized_cellType/ID6_Panel1_TP1.fcs
#> Reading Data/Normalized_cellType/ID6_Panel1_TP2.fcs
#> Reading Data/Normalized_cellType/ID6_Panel1_TP3.fcs
#> Reading Data/Normalized_cellType/ID7_Panel1_TP1.fcs
#> Reading Data/Normalized_cellType/ID7_Panel1_TP2.fcs
#> Reading Data/Normalized_cellType/ID7_Panel1_TP3.fcs
#> Reading Data/Normalized_cellType/ID8_Panel1_TP1.fcs
#> Reading Data/Normalized_cellType/ID8_Panel1_TP2.fcs
#> Reading Data/Normalized_cellType/ID8_Panel1_TP3.fcs
#> Reading Data/Preprocessed/ID5_Panel2_TP1.fcs
#> Reading Data/Preprocessed/ID5_Panel2_TP2.fcs
#> Reading Data/Preprocessed/ID5_Panel2_TP3.fcs
#> Reading Data/Preprocessed/ID6_Panel2_TP1.fcs
#> Reading Data/Preprocessed/ID6_Panel2_TP2.fcs
#> Reading Data/Preprocessed/ID6_Panel2_TP3.fcs
#> Reading Data/Preprocessed/ID7_Panel2_TP1.fcs
#> Reading Data/Preprocessed/ID7_Panel2_TP2.fcs
#> Reading Data/Preprocessed/ID7_Panel2_TP3.fcs
#> Reading Data/Preprocessed/ID8_Panel2_TP1.fcs
#> Reading Data/Preprocessed/ID8_Panel2_TP2.fcs
#> Reading Data/Preprocessed/ID8_Panel2_TP3.fcs
#> Reading Data/Normalized_cellType/ID5_Panel2_TP1.fcs
#> Reading Data/Normalized_cellType/ID5_Panel2_TP2.fcs
#> Reading Data/Normalized_cellType/ID5_Panel2_TP3.fcs
#> Reading Data/Normalized_cellType/ID6_Panel2_TP1.fcs
#> Reading Data/Normalized_cellType/ID6_Panel2_TP2.fcs
#> Reading Data/Normalized_cellType/ID6_Panel2_TP3.fcs
#> Reading Data/Normalized_cellType/ID7_Panel2_TP1.fcs
#> Reading Data/Normalized_cellType/ID7_Panel2_TP2.fcs
#> Reading Data/Normalized_cellType/ID7_Panel2_TP3.fcs
#> Reading Data/Normalized_cellType/ID8_Panel2_TP1.fcs
#> Reading Data/Normalized_cellType/ID8_Panel2_TP2.fcs
#> Reading Data/Normalized_cellType/ID8_Panel2_TP3.fcs
#> [1] "original"
#> Warning in FlowSOM::NewData(model$fsom, as.matrix(dfs[[type]][model$fsom$map$colsUsed])): 19321 cells (8.05%) seem far from their cluster
#> centers.
#> [1] "normalized"
```

``` r

# Arrange figure
p_ <- ggarrange(ggarrange(plotlist = p[1:length(p)-1], 
                          ncol = (length(p)-1)/(2*length(cellTypeChannels)), 
                          nrow = 2*length(cellTypeChannels)),
                widths = c(1,3),
                p$legend, nrow = 2, heights = c(10,1))
print(p_)
```

![](CytoNorm_towards_distribution_files/figure-gfm/Evaluate%20densities-1.png)<!-- -->

### Spline plots

``` r
plotSplines(model = model2_towardsDistribution, channels = cellTypeChannels, groupClusters = TRUE)
#> $P1
#> Warning: Removed 8 rows containing missing values or values outside the scale range (`geom_line()`).
```

![](CytoNorm_towards_distribution_files/figure-gfm/Evaluate%20splines-1.png)<!-- -->

    #> 
    #> $P2
    #> Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_line()`).

![](CytoNorm_towards_distribution_files/figure-gfm/Evaluate%20splines-2.png)<!-- -->

### EMD and MAD

#### Make aggregates and get aggregate manual labels

``` r
dir.create("tmp/before")
#> Warning in dir.create("tmp/before"): 'tmp/before' already exists
```

``` r
dir.create("tmp/after")
#> Warning in dir.create("tmp/after"): 'tmp/after' already exists
```

``` r

aggregates <- list("agg_P1_C2" = files_P1_C2,
                   "agg_P2_C2" = files_P2_C2)

for (agg in names(aggregates)){
  writeLines(agg)
  files <- aggregates[[agg]]
  set.seed(2023)
  before <- AggregateFlowFrames(fileNames = files, silent = TRUE,
                                cTotal = length(files) * 10000)

  write.FCS(before, paste0("tmp/before/", agg, ".fcs"))  
  after <- before
  
  manual <- c()
  for (i_file in unique(before@exprs[,"File"])){
    ff <- read.FCS(paste0("Data/Normalized_cellType/", basename(aggregates[[agg]][i_file])))
    i_IDs <- before@exprs[before@exprs[,"File"] == i_file,"Original_ID2"]
    after@exprs[before@exprs[,"File"] == i_file,1:(ncol(after@exprs)-3)] <- ff@exprs[i_IDs, ]
    labels <- as.character(manualLabels[[basename(aggregates[[agg]][i_file])]])
    manual <- c(manual, labels[i_IDs])
  }    
  
  manualLabels[[paste0(agg, ".fcs")]] <- factor(manual, levels = levels(manualLabels[[1]]))
  write.FCS(after, paste0("tmp/after/", agg, ".fcs"))
}
#> agg_P1_C2
#> agg_P2_C2
```

#### Define cell types and files for model evaluation

``` r
cellTypes <- c("unlabeled",
               "All B cells",
               "Transitional B cells",
               "Naive B cells",
               "IgD+ memory B cells",
               "IgM+ memory B cells", 
               "Class switched memory B cells", 
               "Double negative B cells")

models <- list("Model2_before" = list("files" = c(files_P1_C1_norm, files_P2_C1_norm, 
                                                   "tmp/after/agg_P1_C1.fcs", "tmp/after/agg_P2_C1.fcs",
                                                   "tmp/before/agg_P1_C2.fcs", "tmp/before/agg_P2_C2.fcs", 
                                                   files_P1_C2, files_P2_C2),
                                       "channels" = cellTypeChannels),
                "Model2_after" = list("files" = c(files_P1_C1_norm, files_P2_C1_norm, 
                                                  "tmp/after/agg_P1_C1.fcs", "tmp/after/agg_P2_C1.fcs",
                                                  "tmp/after/agg_P1_C2.fcs", "tmp/after/agg_P2_C2.fcs", 
                                                  files_P1_C2_norm, files_P2_C2_norm),
                                       "channels" = cellTypeChannels))
```

#### Calculate EMD and MAD

``` r
EMDs <- list()
MADs <- list()

for (model in names(models)){
  print(model)
  files <- models[[model]][["files"]]
  channels <- models[[model]][["channels"]]
  
  EMDs[[model]] <- emdEvaluation(files = files,
                                 channels = channels,
                                 manual = manualLabels, return_all = TRUE, 
                                 manualThreshold = 20)
  MADs[[model]] <- CytoNorm:::madEvaluation(files = files[!startsWith(files, "tmp")], 
                                 channels = channels, return_all = TRUE,
                                 manual = manualLabels, transformList = NULL)
}
#> [1] "Model2_before"
#> [1] "Model2_after"
```

``` r

EMD_scores <- list()
for (subset in names(models)){
  dist <- EMDs[[subset]]$distances
  EMD_scores[[subset]] <- matrix(NA,
                                 nrow = length(dist), ncol = length(dist[[1]]),
                                 dimnames = list(names(dist),
                                                 names(dist[[1]])))
  for (cellType in names(dist)){
      for (marker in names(dist[[cellType]])){
        EMD_scores[[subset]][cellType, marker] <- max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T)
      }
  }
}
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
#> Warning in max(dist[[cellType]][[marker]][25:26, 27:28], na.rm = T): no non-missing arguments to max; returning -Inf
```

#### Make figures

``` r
markerlevels <- NULL
plotlist <- list()

# EMD plot
EMD_before <- EMD_scores[["Model2_before"]][-2,]
EMD_after <- EMD_scores[["Model2_after"]][-2,]

EMD_df <- data.frame("before" = c(EMD_before),
                     "after" = c(EMD_after),
                     "feature" = paste(rownames(EMD_before), 
                                       rep(colnames(EMD_before), 
                                           each = nrow(EMD_before)), sep = "_"))
if (is.null(markerlevels)){
  markerlevels <- unique(unique(sub(".*_", "", EMD_df$feature)))
}
EMD_df$marker <- factor(sub(".*_", "", EMD_df$feature), levels = markerlevels)
EMD_df$celltype <- factor(sub("_.*", "", EMD_df$feature), levels = c("AllCells", cellTypes))
max <- max(EMD_df[,1:2], na.rm = T)
p_emd <- ggplot(EMD_df, aes(x = after, y = before)) +
  xlim(0,max) + ylim(0,max) +
  geom_point(aes(color = marker, shape = celltype)) +
  geom_abline() +
  scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 6)) +
  ggtitle("EMD") +
  theme_minimal()

leg_emd <- get_legend(p_emd)

p_emd <- p_emd  +
    theme(legend.position = "none")
plotlist[[length(plotlist)+1]] <- p_emd    

names(markerlevels) <- sub(".*<(.*)>.*", "\\1", markerlevels)

# MAD plot
MAD_before <- MADs[["Model2_before"]][["comparison"]][-2,]
MAD_after <- MADs[["Model2_after"]][["comparison"]][-2,]
MAD_df <- data.frame("before" = c(MAD_before),
                     "after" = c(MAD_after),
                     "feature" = paste(rownames(MAD_before), 
                                       rep(colnames(MAD_before), 
                                           each = nrow(MAD_before)), sep = "_"))
MAD_df$marker <- factor(markerlevels[sub(".*_", "", MAD_df$feature)], levels = markerlevels)
MAD_df$celltype <- factor(sub("_.*", "", MAD_df$feature), c("AllCells", cellTypes))
max <- max(MAD_df[,1:2], na.rm = T)
p_mad <- ggplot(MAD_df, aes(x = after, y = before)) +
  xlim(0,max) + ylim(0,max) +
  geom_point(aes(color = marker, shape = celltype)) +
  geom_abline() +
  scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 6)) +
  ggtitle("MAD") +
  theme_minimal()

leg_mad <- get_legend(p_mad)

p_mad <- p_mad  +
    theme(legend.position = "none")

plotlist[[length(plotlist)+1]] <- p_mad

ggarrange(ggarrange(plotlist = plotlist),
          as_ggplot(leg_emd), ncol=2, widths = c(4,1.5))
```

![](CytoNorm_towards_distribution_files/figure-gfm/EMD%20and%20MAD%20figures-1.png)<!-- -->

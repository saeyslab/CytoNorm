Use case 1: Run CytoNorm on a dataset without dedicated controls
================

<!-- github markdown built using
rmarkdown::render("vignettes/CytoNorm_without_controls.Rmd", output_format = "github_document")
-->

#### In this vignette, you learn how CytoNorm can be used to normalize a dataset where there are no controls available.

Ideally, controls should have been taken along over the course of the
study to have a better understanding of the technical variation and to
train and/or evaluate normalization models eventually. These controls
should be the same sample measured multiple times so that the variation
present is only the result of technical issues. However, sometimes it is
not feasible to include dedicated controls. For example, if a study runs
over a long period of time and concerns fresh samples, it is often not
possible to measure an additional control each time. Other times, it
might be the case that there is simply no appropriate control available;
think of rare samples or studies with a specific stimulation or
treatment. In the original CytoNorm publication, it was explained that
CytoNorm requires acontrol sample per batch on which the normalization
model is trained. Although we believe that using dedicated controls to
build and validate the normalization model is preferable, we understand
that control samples are not always available and that this limits the
applicability of CytoNorm.

Here, we explain how CytoNorm can still be used even without controls.
The key is to treat an aggregate of the batch as a proxy for a control
sample of that batch. Either the aggregates are built internally when
multiple files per label are given when training the model, or they can
be created prior to the CytoNorm call. This has the advantage of more
control on the number of cells per file.

## Prepare the CytoNorm normalization step

### Load CytoNorm and other necessary packages and create directories

``` r
library(CytoNorm)
library(flowCore)
library(FlowSOM)
#> Loading required package: igraph
#> 
#> Attaching package: 'igraph'
#> The following object is masked from 'package:flowCore':
#> 
#>     normalize
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> The following object is masked from 'package:base':
#> 
#>     union
#> Thanks for using FlowSOM. From version 2.1.4 on, the scale 
#> parameter in the FlowSOM function defaults to FALSE
```

``` r
library(ggpubr)
#> Loading required package: ggplot2
#> Want to understand how all the pieces fit together? Read R for Data Science: https://r4ds.hadley.nz/
```

``` r

dir.create("RDS")
```

### Download and organize the data

As example dataset, we will use flow cytometry PBMC data that was
collected in the context of a clinical study to evaluate adjuvant
dendritic cell based vaccination for malignant peritoneal mesothelioma.
The raw fcs files were compensated and transformed, and margin events,
low quality events and non-B cells were removed using a preprocessing
workflow that was described in [Analyzing high-dimensional cytometry
data using FlowSOM](https://rdcu.be/cndgZ).The preprocessed and cleaned
dataset can be downloaded from flowRepository dataset
[FR-FCM-Z6UT](http://flowrepository.org/experiments/7133). Once the data
is downloaded, put it in a folder in your working directory,
e.g. “Data/Preprocessed”.

The dataset consists of 48 fcs files (8 patients, 2 panels and 3
consecutive time points per patient-panel combination) and there are two
levels of batch effects:

- Panel, both panels have common set of backbone markers
  - B cell panel 1 (P1)
  - B cell panel 2 (P2)
- Patient cohort, corresponds with analysis day
  - Inclusion cohort 1 (C1)
  - Inclusion cohort 2 (C2)

This results in four batches of data with 12 fcs files each (P1_C1,
P1_C2, P2_C1 and P2_C2).

``` r
dir_prepr <- "Data/Preprocessed"

# Organize preprocessed data in batches
files_P1_C1 <- list.files(path = dir_prepr,
                          pattern = "ID[1-4]_Panel1_TP[1-3].fcs",
                          full.names = TRUE)
files_P2_C1 <- list.files(path = dir_prepr,
                          pattern = "ID[1-4]_Panel2_TP[1-3].fcs",
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
```

In this vignette, we will normalize the cell type markers of the two
panels of the first patient cohort towards each other, i.e. P1_C1 and
P2_C1.

## Train the CytoNorm model

### Create aggregate files that will be used as a proxi for a control file

``` r
# Create aggregate per batch with 10000 cells from each file
set.seed(1)
agg_P1_C1 <- AggregateFlowFrames(fileNames = files_P1_C1,
                                 cTotal = length(files_P1_C1)*10000) 
#> Reading Data/Preprocessed/ID1_Panel1_TP1.fcs
#> Reading Data/Preprocessed/ID1_Panel1_TP2.fcs
#> Reading Data/Preprocessed/ID1_Panel1_TP3.fcs
#> Reading Data/Preprocessed/ID2_Panel1_TP1.fcs
#> Reading Data/Preprocessed/ID2_Panel1_TP2.fcs
#> Reading Data/Preprocessed/ID2_Panel1_TP3.fcs
#> Reading Data/Preprocessed/ID3_Panel1_TP1.fcs
#> Reading Data/Preprocessed/ID3_Panel1_TP2.fcs
#> Reading Data/Preprocessed/ID3_Panel1_TP3.fcs
#> Reading Data/Preprocessed/ID4_Panel1_TP1.fcs
#> Reading Data/Preprocessed/ID4_Panel1_TP2.fcs
#> Reading Data/Preprocessed/ID4_Panel1_TP3.fcs
```

``` r
agg_P2_C1 <- AggregateFlowFrames(fileNames = files_P2_C1,
                                 cTotal = length(files_P2_C1)*10000)
#> Reading Data/Preprocessed/ID1_Panel2_TP1.fcs
#> Reading Data/Preprocessed/ID1_Panel2_TP2.fcs
#> Reading Data/Preprocessed/ID1_Panel2_TP3.fcs
#> Reading Data/Preprocessed/ID2_Panel2_TP1.fcs
#> Reading Data/Preprocessed/ID2_Panel2_TP2.fcs
#> Reading Data/Preprocessed/ID2_Panel2_TP3.fcs
#> Reading Data/Preprocessed/ID3_Panel2_TP1.fcs
#> Reading Data/Preprocessed/ID3_Panel2_TP2.fcs
#> Reading Data/Preprocessed/ID3_Panel2_TP3.fcs
#> Reading Data/Preprocessed/ID4_Panel2_TP1.fcs
#> Reading Data/Preprocessed/ID4_Panel2_TP2.fcs
#> Reading Data/Preprocessed/ID4_Panel2_TP3.fcs
```

``` r

# Save aggregates
write.FCS(agg_P1_C1, "Data/Preprocessed/agg_P1_C1.fcs")
#> [1] "Data/Preprocessed/agg_P1_C1.fcs"
```

``` r
write.FCS(agg_P2_C1, "Data/Preprocessed/agg_P2_C1.fcs")
#> [1] "Data/Preprocessed/agg_P2_C1.fcs"
```

### Train the CytoNorm model

``` r
model1_withoutControl <- CytoNorm.train(files = flowSet(agg_P1_C1, agg_P2_C1),
                                        labels = c("P1", "P2"),
                                        channels = cellTypeChannels,
                                        transformList = NULL,
                                        seed = 1,
                                        verbose = TRUE,
                                        plot = TRUE,
                                        FlowSOM.params = list(nCells = 1e+06, 
                                                              xdim = 7, ydim = 7, 
                                                              nClus = 3, 
                                                              scale = FALSE))
#> Plot FlowSOM trees
#> Plot file distribution
#> Calculate t-SNE
#> Plot cluster per metacluster distribution
#> Plot heatmap
#> Make tables
#> Printing
#> Splitting V1
#> Splitting V2
#> Processing cluster 1
#> Computing Quantiles
#>   P1 (FileID 1)
#> 
#>   P2 (FileID 2)
#> 
#> Computing Splines
#>   P1
#>   P2
#> Processing cluster 2
#> Computing Quantiles
#>   P1 (FileID 1)
#> 
#>   P2 (FileID 2)
#> 
#> Computing Splines
#>   P1
#>   P2
#> Processing cluster 3
#> Computing Quantiles
#>   P1 (FileID 1)
#> 
#>   P2 (FileID 2)
#> 
#> Computing Splines
#>   P1
#>   P2
```

``` r
saveRDS(object = model1_withoutControl, file = "RDS/model1_withoutControl.rds")
```

The `files` parameter in the `CytoNorm.train` call requires the paths to
the control fcs files or a flowSet with the control samples in it and
the `labels` argument is used to define the batch labels for the control
files. Here however, an aggregate or concatenated fcs file is given per
batch instead of a control file. The `channels` that should be
normalized should be specified in the channels argument and -if needed-
a transformation object to transform the data can be given at
`transformList`. Optionally, a seed can be determined to make the
results reproducible, progress updates can be printed via `verbose` and
`plot = TRUE` will generate and save visualizations of the internal
FlowSOMmodel and CytoNorm splines. The dimensions of the self-organizing
map (SOM) of the FlowSOM model, number of cells to use and number of
meta-clusters, should be specified in the `FlowSOM.params` argument.
Here, we had a moderate number of markers for the clustering and
normalization and based on our panel we did not expect many distinct
cell populations. Therefore, and after inspection of the visualizations
listed below, we opted for a seven-by-seven SOM grid and three
meta-clusters. In contrast, the `files` argument in
`CytoNorm.normalize()` takes all the files that need to be normalized
with the trained model that is provided to the model parameter. Batch
identities per file should be specified in labels and similar to the
function to train the model, progress updates can be printed via
verbose. If a transformation object is provided at `transformList`, an
object to reverse the transformation should also be given at
`transformList.reverse` so that the fcs files can be saved in the
original space. A prefix can be specified to mark the normalized fcs
files when they are saved in the `outputDir`.

### Evaluate the CytoNorm assumption: check the CVs

``` r
# Evaluate the CVs
CVs1 <- testCV(fsom = model1_withoutControl$fsom, 
               cluster_values = 3:20, plot = FALSE)

PlotOverviewCV(fsom = model1_withoutControl$fsom, cv_res = CVs1,
               show_cv = 0.8, max_cv = 1.5)
```

![](CytoNorm_without_controls_files/figure-gfm/Evaluate%20CVs-1.png)<!-- -->![](CytoNorm_without_controls_files/figure-gfm/Evaluate%20CVs-2.png)<!-- -->![](CytoNorm_without_controls_files/figure-gfm/Evaluate%20CVs-3.png)<!-- -->

## Apply CytoNorm model

``` r
# Normalize files
CytoNorm.normalize(model = model1_withoutControl,
                   files = c(files_P1_C1, files_P2_C1),
                   labels = c(rep("P1", length(files_P1_C1)),
                              rep("P2", length(files_P2_C1))),
                   transformList = NULL,
                   verbose = TRUE,
                   prefix = "",
                   transformList.reverse = NULL, 
                   outputDir = "Data/Normalized_cellType")
#> Splitting Data/Preprocessed/ID1_Panel1_TP1.fcs
#> Splitting Data/Preprocessed/ID1_Panel1_TP2.fcs
#> Splitting Data/Preprocessed/ID1_Panel1_TP3.fcs
#> Splitting Data/Preprocessed/ID2_Panel1_TP1.fcs
#> Splitting Data/Preprocessed/ID2_Panel1_TP2.fcs
#> Splitting Data/Preprocessed/ID2_Panel1_TP3.fcs
#> Splitting Data/Preprocessed/ID3_Panel1_TP1.fcs
#> Splitting Data/Preprocessed/ID3_Panel1_TP2.fcs
#> Splitting Data/Preprocessed/ID3_Panel1_TP3.fcs
#> Splitting Data/Preprocessed/ID4_Panel1_TP1.fcs
#> Splitting Data/Preprocessed/ID4_Panel1_TP2.fcs
#> Splitting Data/Preprocessed/ID4_Panel1_TP3.fcs
#> Splitting Data/Preprocessed/ID1_Panel2_TP1.fcs
#> Splitting Data/Preprocessed/ID1_Panel2_TP2.fcs
#> Splitting Data/Preprocessed/ID1_Panel2_TP3.fcs
#> Splitting Data/Preprocessed/ID2_Panel2_TP1.fcs
#> Splitting Data/Preprocessed/ID2_Panel2_TP2.fcs
#> Splitting Data/Preprocessed/ID2_Panel2_TP3.fcs
#> Splitting Data/Preprocessed/ID3_Panel2_TP1.fcs
#> Splitting Data/Preprocessed/ID3_Panel2_TP2.fcs
#> Splitting Data/Preprocessed/ID3_Panel2_TP3.fcs
#> Splitting Data/Preprocessed/ID4_Panel2_TP1.fcs
#> Splitting Data/Preprocessed/ID4_Panel2_TP2.fcs
#> Splitting Data/Preprocessed/ID4_Panel2_TP3.fcs
#> Processing cluster 1
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel1_TP1.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel1_TP2.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel1_TP3.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel1_TP1.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel1_TP2.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel1_TP3.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel1_TP1.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel1_TP2.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel1_TP3.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel1_TP1.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel1_TP2.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel1_TP3.fcs_fsom1.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel2_TP1.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel2_TP2.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel2_TP3.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel2_TP1.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel2_TP2.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel2_TP3.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel2_TP1.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel2_TP2.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel2_TP3.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel2_TP1.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel2_TP2.fcs_fsom1.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel2_TP3.fcs_fsom1.fcs (P2)
#> Normalizing P2
#> Processing cluster 2
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel1_TP1.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel1_TP2.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel1_TP3.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel1_TP1.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel1_TP2.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel1_TP3.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel1_TP1.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel1_TP2.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel1_TP3.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel1_TP1.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel1_TP2.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel1_TP3.fcs_fsom2.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel2_TP1.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel2_TP2.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel2_TP3.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel2_TP1.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel2_TP2.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel2_TP3.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel2_TP1.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel2_TP2.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel2_TP3.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel2_TP1.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel2_TP2.fcs_fsom2.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel2_TP3.fcs_fsom2.fcs (P2)
#> Normalizing P2
#> Processing cluster 3
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel1_TP1.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel1_TP2.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel1_TP3.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel1_TP1.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel1_TP2.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel1_TP3.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel1_TP1.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel1_TP2.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel1_TP3.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel1_TP1.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel1_TP2.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel1_TP3.fcs_fsom3.fcs (P1)
#> Normalizing P1
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel2_TP1.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel2_TP2.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID1_Panel2_TP3.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel2_TP1.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel2_TP2.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID2_Panel2_TP3.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel2_TP1.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel2_TP2.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID3_Panel2_TP3.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel2_TP1.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel2_TP2.fcs_fsom3.fcs (P2)
#> Normalizing P2
#>   Data/Normalized_cellType/Data_Preprocessed_ID4_Panel2_TP3.fcs_fsom3.fcs (P2)
#> Normalizing P2
#> Rebuilding Data/Preprocessed/ID1_Panel1_TP1.fcs
#> Rebuilding Data/Preprocessed/ID1_Panel1_TP2.fcs
#> Rebuilding Data/Preprocessed/ID1_Panel1_TP3.fcs
#> Rebuilding Data/Preprocessed/ID2_Panel1_TP1.fcs
#> Rebuilding Data/Preprocessed/ID2_Panel1_TP2.fcs
#> Rebuilding Data/Preprocessed/ID2_Panel1_TP3.fcs
#> Rebuilding Data/Preprocessed/ID3_Panel1_TP1.fcs
#> Rebuilding Data/Preprocessed/ID3_Panel1_TP2.fcs
#> Rebuilding Data/Preprocessed/ID3_Panel1_TP3.fcs
#> Rebuilding Data/Preprocessed/ID4_Panel1_TP1.fcs
#> Rebuilding Data/Preprocessed/ID4_Panel1_TP2.fcs
#> Rebuilding Data/Preprocessed/ID4_Panel1_TP3.fcs
#> Rebuilding Data/Preprocessed/ID1_Panel2_TP1.fcs
#> Rebuilding Data/Preprocessed/ID1_Panel2_TP2.fcs
#> Rebuilding Data/Preprocessed/ID1_Panel2_TP3.fcs
#> Rebuilding Data/Preprocessed/ID2_Panel2_TP1.fcs
#> Rebuilding Data/Preprocessed/ID2_Panel2_TP2.fcs
#> Rebuilding Data/Preprocessed/ID2_Panel2_TP3.fcs
#> Rebuilding Data/Preprocessed/ID3_Panel2_TP1.fcs
#> Rebuilding Data/Preprocessed/ID3_Panel2_TP2.fcs
#> Rebuilding Data/Preprocessed/ID3_Panel2_TP3.fcs
#> Rebuilding Data/Preprocessed/ID4_Panel2_TP1.fcs
#> Rebuilding Data/Preprocessed/ID4_Panel2_TP2.fcs
#> Rebuilding Data/Preprocessed/ID4_Panel2_TP3.fcs
```

## Evaluate the quality

### Density plots

``` r
# List files
files_P1_C1_norm <- list.files(path = "Data/Normalized_cellType", 
                               pattern = "ID[1-4]_Panel1_TP[1-3].fcs", full.names = TRUE)
files_P2_C1_norm <- list.files(path = "Data/Normalized_cellType", 
                               pattern = "ID[1-4]_Panel2_TP[1-3].fcs", full.names = TRUE)

# Make plots
p <- plotDensities(input = list("P1_C1" = agg_P1_C1,
                                "P2_C1" = agg_P2_C1,
                                "P1_C1_norm" = files_P1_C1_norm,
                                "P2_C1_norm" = files_P2_C1_norm),
                   channels = cellTypeChannels,
                   colors = c("P1_C1" = "#900002", 
                              "P2_C1" = "#0046A1"), 
                   model = model1_withoutControl)
#> Reading Data/Normalized_cellType/ID1_Panel1_TP1.fcs
#> Reading Data/Normalized_cellType/ID1_Panel1_TP2.fcs
#> Reading Data/Normalized_cellType/ID1_Panel1_TP3.fcs
#> Reading Data/Normalized_cellType/ID2_Panel1_TP1.fcs
#> Reading Data/Normalized_cellType/ID2_Panel1_TP2.fcs
#> Reading Data/Normalized_cellType/ID2_Panel1_TP3.fcs
#> Reading Data/Normalized_cellType/ID3_Panel1_TP1.fcs
#> Reading Data/Normalized_cellType/ID3_Panel1_TP2.fcs
#> Reading Data/Normalized_cellType/ID3_Panel1_TP3.fcs
#> Reading Data/Normalized_cellType/ID4_Panel1_TP1.fcs
#> Reading Data/Normalized_cellType/ID4_Panel1_TP2.fcs
#> Reading Data/Normalized_cellType/ID4_Panel1_TP3.fcs
#> Reading Data/Normalized_cellType/ID1_Panel2_TP1.fcs
#> Reading Data/Normalized_cellType/ID1_Panel2_TP2.fcs
#> Reading Data/Normalized_cellType/ID1_Panel2_TP3.fcs
#> Reading Data/Normalized_cellType/ID2_Panel2_TP1.fcs
#> Reading Data/Normalized_cellType/ID2_Panel2_TP2.fcs
#> Reading Data/Normalized_cellType/ID2_Panel2_TP3.fcs
#> Reading Data/Normalized_cellType/ID3_Panel2_TP1.fcs
#> Reading Data/Normalized_cellType/ID3_Panel2_TP2.fcs
#> Reading Data/Normalized_cellType/ID3_Panel2_TP3.fcs
#> Reading Data/Normalized_cellType/ID4_Panel2_TP1.fcs
#> Reading Data/Normalized_cellType/ID4_Panel2_TP2.fcs
#> Reading Data/Normalized_cellType/ID4_Panel2_TP3.fcs
#> [1] "original"
#> Warning in FlowSOM::NewData(model$fsom, as.matrix(dfs[[type]][model$fsom$map$colsUsed])): 1255 cells (1.14%) seem far from their cluster centers.
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

![](CytoNorm_without_controls_files/figure-gfm/Evaluate%20densities-1.png)<!-- -->

### Spline plots

``` r
plotSplines(model = model1_withoutControl, channels = cellTypeChannels, 
            groupClusters = TRUE)
#> $P1
#> Warning: Removed 42 rows containing missing values or values outside the scale range (`geom_line()`).
```

![](CytoNorm_without_controls_files/figure-gfm/Evaluate%20splines-1.png)<!-- -->

    #> 
    #> $P2
    #> Warning: Removed 13 rows containing missing values or values outside the scale range (`geom_line()`).

![](CytoNorm_without_controls_files/figure-gfm/Evaluate%20splines-2.png)<!-- -->

``` r
plotSplines(model = model1_withoutControl, channels = cellTypeChannels)
#> $P1
#> Warning: Removed 2 rows containing missing values or values outside the scale range (`geom_line()`).
```

![](CytoNorm_without_controls_files/figure-gfm/Evaluate%20splines-3.png)<!-- -->

    #> 
    #> $P2
    #> Warning: Removed 11 rows containing missing values or values outside the scale range (`geom_line()`).

![](CytoNorm_without_controls_files/figure-gfm/Evaluate%20splines-4.png)<!-- -->

### EMD and MAD

#### Make aggregates and get aggregate manual labels

``` r
# Get manual labels per preprocessed file
gating <- readRDS("Data/Preprocessed/attachments/gating.rds")
manualLabels <- list()

dictionary <- c("Unlabeled" = "unlabeled",
                "Trans" = "Transitional B cells",
                "Q1: CD27- , IgD+ (naive B)" = "Naive B cells",
                "Q2: CD27+ , IgD+ (mem)" = "IgD+ memory B cells",
                "CD27, IgM+ (mem)" = "IgM+ memory B cells", 
                "CD27, IgM- (class switched mem)" = "Class switched memory B cells", 
                "Q4: CD27- , IgD-" = "Double negative B cells")

for (file in names(gating)){
  manual <- ManualVector(manualMatrix = gating[[file]], 
                         cellTypes = names(dictionary)[-1])
  manualLabels[[paste0(file, ".fcs")]] <- factor(unname(dictionary[manual]),                                                 
                                                 levels = unname(dictionary))
}
saveRDS(manualLabels, "Data/Preprocessed/attachments/manualLabels.rds")

# Aggregates
dir.create("tmp/before")
dir.create("tmp/after")
manual_agg <- list()
manualLabels <- readRDS("Data/Preprocessed/attachments/manualLabels.rds")

aggregates <- list("agg_P1_C1" = files_P1_C1,
                   "agg_P2_C1" = files_P2_C1)

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
#> agg_P1_C1
#> agg_P2_C1
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

models <- list("Model1_before" = list("files" = c(files_P1_C1, "tmp/before/agg_P1_C1.fcs", 
                                                   "tmp/before/agg_P2_C1.fcs", files_P2_C1),
                                       "channels" = cellTypeChannels),
                "Model1_after" = list("files" = c(files_P1_C1_norm, "tmp/after/agg_P1_C1.fcs", 
                                                  "tmp/after/agg_P2_C1.fcs", files_P2_C1_norm),
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
#> [1] "Model1_before"
#> [1] "Model1_after"
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
      EMD_scores[[subset]][cellType, marker] <- dist[[cellType]][[marker]][13, 14]
      }
  }
}
```

#### Make figures

``` r
markerlevels <- NULL
plotlist <- list()

# EMD plot
EMD_before <- EMD_scores[["Model1_before"]][-2,]
EMD_after <- EMD_scores[["Model1_after"]][-2,]

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
MAD_before <- MADs[["Model1_before"]][["comparison"]][-2,]
MAD_after <- MADs[["Model1_after"]][["comparison"]][-2,]
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

![](CytoNorm_without_controls_files/figure-gfm/EMD%20and%20MAD%20figures-1.png)<!-- -->

#' prepareFlowSOM
#'
#' Aggregate files, transform them and run FlowSOM.
#' This is used as a first step in the normalization procedure to detect
#' groups of similar cells. Typically you will not call this function, but use
#' the wrapper function \code{\link{CytoNorm.train}} instead.
#'
#' @param files     Files to use (with path)
#' @param colsToUse IDs or column names of the columns to use for clustering
#' @param nCells    The total number of cells to use for the FlowSOM clustering.
#'                  This number is divided by the number of files to determine
#'                  the amount to select from each individual file.
#'                  Default = 1 000 000.
#' @param FlowSOM.params List with parameters to pass to the FlowSOM algorithm.
#'                       Default = \code{list(xdim = 15, ydim = 15, nclus = 30,
#'                       scale = FALSE)}.
#' @param transformList Transformation list to pass to the flowCore
#'                      \code{transform} function. Default = NULL.
#' @param plot      If TRUE, the FlowSOM tree is plotted. Default = FALSE.
#' @param verbose   If TRUE, extra output is printed while running.
#' @param seed      If not NULL, set.seed is called with this argument for
#'                  reproducable results. Default = NULL.
#'
#' @return FlowSOM object describing the FlowSOM clustering
#'
#' @examples
#'
#'
#' @importFrom dplyr '%>%' filter
#' @importFrom flowCore read.FCS transformList colnames
#' @importFrom stringr str_match
#'
#' @export
prepareFlowSOM <- function(files,
                           colsToUse,
                           nCells = 1000000,
                           FlowSOM.params = list(xdim = 15,
                                                 ydim = 15,
                                                 nClus = 30,
                                                 scale = FALSE),
                           transformList = NULL,
                           plot = FALSE,
                           verbose = FALSE,
                           seed = NULL){

    if (verbose) message("Aggregating files ... ")

    if(!is.null(seed)) set.seed(seed)
    o <- capture.output( ff <- FlowSOM::AggregateFlowFrames(files, nCells))
    if(!is.null(transformList)) ff <- flowCore::transform(ff, transformList)

    FlowSOM.params <- c(FlowSOM.params,
                        list(input = ff,
                             colsToUse = colsToUse,
                             seed = seed))
    FlowSOM.params <- FlowSOM.params[unique(names(FlowSOM.params))]

    if (verbose) message("Running the FlowSOM algorithm ... ")
    fsom <- do.call(FlowSOM::FlowSOM, FlowSOM.params)

    if (plot) {
        # Plot result
        FlowSOM::PlotStars(fsom$FlowSOM,
                           backgroundValues = fsom$metaclustering)
    }

    fsom
}

#' CytoNorm.train
#'
#' CytoNorm learns batch effects, using control samples which should be the
#' the same per batch. The cells are first split into rough groups using the
#' FlowSOM algorithm, after which splines are trained to map the quantiles for
#' each control group to an average distribution
#' Typically, you will use the function \code{\link{CytoNorm.normalize}} after
#' this function to reverse the batch effects in other files.
#'
#' Temporary fcs files splitting the data into clusters are written to
#' \code{outputDir} and will be removed again by default (depending on
#' \code{clean}).
#'
#' @param files       Full paths of to the fcs files of the control samples.
#' @param labels      A label for every file, indicating to which batch it
#'                    belongs, e.g. the plate ID.
#' @param channels    Column names of the channels that need to be normalized
#' @param transformList   Transformation list to pass to the flowCore
#'                        \code{transform} function
#' @param outputDir   Directory to put the temporary files in. Default = "./tmp"
#' @param clean       Whether to remove the temporary files again at the end.
#'                    Default = TRUE.
#' @param plot        If TRUE, plots are saved to the current dir.
#'                    Default = FALSE.
#' @param verbose     If TRUE, progress updates are printed. Default = FALSE.
#' @param seed        Set a seed for reproducable results.
#' @param FlowSOM.params    Extra parameters to be passed to the FlowSOM
#'                          function, such as the number of cells
#'                          (\code{nCells}), the \code{xdim} and \code{y dim} or
#'                          the number of clusters (\code{nClus})
#' @param normMethod.train Normalization method to use for each cluster.
#'                         Default = \code{\link{QuantileNorm.train}}
#' @param normParams Parameters to pass to the normalization method. Default,
#'                   assuming \code{\link{QuantileNorm.train}}:
#'                   list(nQ = 21)). nQ is the number of quantiles
#'                   to use.
#'
#' @return A list containing two elements: the FlowSOM clustering and a list
#'         containing all the splines per cluster. This can be used as input
#'         for the \code{\link{CytoNorm.normalize}} function.
#' @seealso   \code{\link{QuantileNorm.train}}, \code{\link{CytoNorm.normalize}},
#'            \code{\link{prepareFlowSOM}}
#'
#' @examples
#'
#' @export
CytoNorm.train <- function(files,
                            labels,
                            channels,
                            transformList,
                            outputDir = "./tmp",
                            FlowSOM.params = list(nCells = 1000000,
                                                  xdim = 15,
                                                  ydim = 15,
                                                  nClus = 30,
                                                  scale = FALSE),
                            normMethod.train = QuantileNorm.train,
                            normParams = list(nQ = 21),
                            seed = NULL,
                            clean = TRUE,
                            plot = FALSE,
                            verbose = FALSE){

    if (length(labels) != length(files)) {
        stop("Input parameters 'labels' and 'files'",
             " should have the same length")
    }

    # Create output directory
    dirCreated = FALSE
    if (!dir.exists(outputDir)) {
        dirCreated = dir.create(outputDir)
    }

    if(!file.exists(file.path(outputDir, "CytoNorm_FlowSOM.RDS"))){
        # Compute clustering
        if (plot) {
            grDevices::pdf(file.path(outputDir, "CytoNorm_FlowSOM.pdf"),
                           height=15, width=20)
        }


        nCells <- FlowSOM.params[["nCells"]]
        if(is.null(FlowSOM.params[["channels"]])){
            FlowSOM.channels <- channels
        } else {
            FlowSOM.channels <- FlowSOM.params[["channels"]]
        }
        FlowSOM.params <- FlowSOM.params[grep("nCells|channels",
                                              names(FlowSOM.params),
                                              invert = TRUE)]
        fsom <- prepareFlowSOM(files = files,
                               nCells = nCells,
                               FlowSOM.params = FlowSOM.params,
                               transformList = transformList,
                               colsToUse = FlowSOM.channels,
                               plot = plot,
                               seed = seed)

        if (plot) {
            grDevices::dev.off()
            saveRDS(fsom, file.path(outputDir, "CytoNorm_FlowSOM.RDS"))
        }
    } else {
        fsom <- readRDS(file.path(outputDir, "CytoNorm_FlowSOM.RDS"))
        warning("Reusing previously saved FlowSOM result.")
    }

    # Split files by clusters
    for(file in files){
        message("Splitting ",file)
        ff <- flowCore::read.FCS(file)
        if (!is.null(transformList)) {
            ff <- flowCore::transform(ff,transformList)
        }

        # Map the file to the FlowSOM clustering
        fsom_file <- FlowSOM::NewData(fsom$FlowSOM, ff)

        # Get the metacluster label for every cell
        cellClusterIDs <- fsom$metaclustering[fsom_file$map$mapping[,1]]
        for (cluster in unique(fsom$metaclustering)) {
            if (sum(FlowSOM::GetMetaclusters(fsom_file,
                                             fsom$metaclustering) == cluster) > 0) {
                suppressWarnings(
                    flowCore::write.FCS(
                        ff[cellClusterIDs == cluster,],
                        file=file.path(outputDir,
                                       paste0(gsub("/","_",file),
                                              "_fsom",cluster,".fcs"))))
            }
        }
    }

    # Learn quantiles for each cluster
    clusterRes <- list()
    for (cluster in unique(fsom$metaclustering)) {
        message("Processing cluster ",cluster)
        if (plot) {
            grDevices::pdf(file.path(outputDir,
                                     paste0("CytoNorm_norm_Cluster",
                                            cluster, ".pdf")),
                           height = 3*(2*length(files)+2),
                           width = 3*(length(channels)+1))
        }

        normParams_tmp <- c(normParams,
                            list(files = file.path(outputDir,
                                                   paste0(gsub("/", "_", files),
                                                          "_fsom", cluster, ".fcs")),
                                 labels = as.character(labels),
                                 channels = channels,
                                 transformList = NULL,
                                 verbose = verbose,
                                 plot = plot))
        normParams_tmp <- normParams_tmp[unique(names(normParams_tmp))]
        clusterRes[[cluster]] <- do.call(normMethod.train,
                                         normParams_tmp)

        if (plot) { grDevices::dev.off() }
    }

    if(clean){
        for(cluster in unique(fsom$metaclustering)){
            tmp_files <- file.path(outputDir,
                                   paste0(gsub("/","_",files),"_fsom",cluster,".fcs"))

            file.remove(tmp_files[file.exists(tmp_files)])
        }
        if(dirCreated & !plot){
            unlink(outputDir, recursive=TRUE)
        }
    }
    named.list(fsom, clusterRes)
}


#' normalizeClustered
#'
#' Normalize data, given the batch effects learned from control samples per
#' cell type/cluster (output from \code{\link{CytoNorm.train}}). New fcs files
#' are written to the given output directory.
#'
#' @param model       Model of the batch effercts, as computed by
#'                    \code{\link{CytoNorm.train}}
#' @param files       Full paths of to the fcs files of the samples.
#' @param labels      A label for every file, indicating to which batch it
#'                    belongs, e.g. the plate ID.
#' @param transformList   Transformation list to pass to the flowCore
#'                    \code{transform} function
#' @param transformList.reverse   Transformation list with the reverse functions,
#'                    so the normalized files can be saved in the untransformed
#'                    space
#' @param outputDir   Directory to put the temporary files in. Default = "."
#' @param prefix      Prefix to put in front of the normalized file names.
#'                    Default = "Norm_"
#' @param verbose     If TRUE, progress updates are printed. Default = FALSE.
#' @param clean       If FALSE, temporary files describing the FlowSOM clusters
#'                    seperately are not removed at the end. Default = TRUE.
#' @param normMethod.normalize Normalization method to use.
#' @return Nothing is returned, but the new FCS files are written to the output
#'         directory
#' @seealso   \code{\link{CytoNorm.train}}
#'
#' @examples
#'
#' @export
CytoNorm.normalize <- function(model,
                                files,
                                labels,
                                transformList,
                                transformList.reverse,
                                outputDir = ".",
                                prefix = "Norm_",
                                clean = TRUE,
                                verbose = FALSE,
                                normMethod.normalize = QuantileNorm.normalize){
    if(is.null(model$fsom) |
       is.null(model$clusterRes)){
        stop("The 'model' paramter should be the result of using the
             trainQuantiles function.")
    }

    if(length(labels) != length(files)){
        stop("Input parameters 'labels' and 'files' should have the same length")
    }

    # Create output directory
    if(!dir.exists(outputDir)){
        dir.create(outputDir)
    }

    fsom <- model$fsom
    clusterRes <- model$clusterRes

    # Split files by clusters
    cellClusterIDs <- list()
    meta <- list()
    cluster_files <- list()
    for(file in files){
        if(verbose) message("Splitting ",file)
        ff <- flowCore::read.FCS(file)

        if(!is.null(transformList)){
            ff <- flowCore::transform(ff, transformList)
            # meta[[file]] <- list()
            # meta[[file]][["description_original"]] <- ff@description
            # meta[[file]][["parameters_original"]] <- ff@parameters
        }

        fsom_file <- FlowSOM::NewData(fsom$FlowSOM,ff)

        cellClusterIDs[[file]] <- fsom$metaclustering[fsom_file$map$mapping[,1]]

        for(cluster in unique(fsom$metaclustering)){
            if (sum(FlowSOM::GetMetaclusters(fsom_file,
                                             fsom$metaclustering) == cluster) > 0) {
                f <- file.path(outputDir,
                               paste0(gsub(".*/","",file),
                                      "_fsom", cluster, ".fcs"))
                suppressWarnings(
                    flowCore::write.FCS(ff[cellClusterIDs[[file]] == cluster],
                                        file = f)
                )
            }
        }
    }

    # Apply normalization on each cluster
    for(cluster in unique(fsom$metaclustering)){
        if(verbose) message("Processing cluster ",cluster)
        files_tmp <- file.path(outputDir,
                               paste0(gsub(".*/",
                                           "",
                                           files),
                                      "_fsom",
                                      cluster,
                                      ".fcs"))
        labels_tmp <- labels[file.exists(files_tmp)]
        files_tmp <- files_tmp[file.exists(files_tmp)]
        normMethod.normalize(model = clusterRes[[cluster]],
                               files = files_tmp,
                               labels = labels_tmp,
                               outputDir = file.path(outputDir),
                               prefix = "Norm_",
                               transformList = NULL,
                               transformList.reverse = NULL,
                               removeOriginal = TRUE,
                               verbose = verbose)
    }

    # Combine clusters into one final fcs file
    for(file in files){
        if(verbose) message("Rebuilding ",file)

        ff <- flowCore::read.FCS(file)
        for(cluster in unique(fsom$metaclustering)){
            file_name <- file.path(outputDir,
                                   paste0("Norm_",gsub(".*/","",file),
                                          "_fsom",cluster,".fcs"))
            if (file.exists(file_name)) {
                ff_subset <- flowCore::read.FCS(file_name)
                flowCore::exprs(ff)[cellClusterIDs[[file]] == cluster,] <- flowCore::exprs(ff_subset)
            }
        }

        if(!is.null(transformList.reverse)){
            ff <- flowCore::transform(ff, transformList.reverse)
            # ff@description <- meta[[file]][["description_original"]]
            # ff@parameters <- meta[[file]][["parameters_original"]]
        }

        if(clean){
            file.remove(file.path(outputDir,
                      paste0("Norm_",gsub(".*/","",file),
                             "_fsom",unique(fsom$metaclustering),".fcs")))
        }

        suppressWarnings(flowCore::write.FCS(ff,
                                             file=file.path(outputDir,
                                                    paste0(prefix,gsub(".*/","",file)))))
    }
}

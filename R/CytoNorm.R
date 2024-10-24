#' prepareFlowSOM
#'
#' Aggregate files, transform them and run FlowSOM.
#' This is used as a first step in the normalization procedure to detect
#' groups of similar cells. Typically you will not call this function, but use
#' the wrapper function \code{\link{CytoNorm.train}} instead.
#'
#' @param files     path to FCS file or a flowSet containing the samples
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
#' @param verbose   If TRUE, extra output is printed while running.
#' @param seed      If not NULL, set.seed is called with this argument for
#'                  reproducable results. Default = NULL.
#' @param ...         Additional arguments to pass to read.FCS
#'
#' @return FlowSOM object describing the FlowSOM clustering
#'
#' @examples
#'
#' dir <- system.file("extdata", package = "CytoNorm")
#' files <- list.files(dir, pattern = "fcs$")
#' data <- data.frame(File = files,
#'                    Path = file.path(dir, files),
#'                    Type = stringr::str_match(files, "_([12]).fcs")[,2],
#'                    Batch = stringr::str_match(files, "PTLG[0-9]*")[,1],
#'                    stringsAsFactors = FALSE)
#' data$Type <- c("1" = "Train", "2" = "Validation")[data$Type]
#' train_data <- dplyr::filter(data, Type == "Train")
#'
#' ff <- flowCore::read.FCS(data$Path[1])
#' channels <- grep("Di$", flowCore::colnames(ff), value = TRUE)
#' transformList <- flowCore::transformList(channels,
#'                                          cytofTransform)
#'
#' fsom <- prepareFlowSOM(train_data$Path,
#'                        channels,
#'                        nCells = 10000, #1000000
#'                        FlowSOM.params = list(xdim = 15,
#'                                              ydim = 15,
#'                                              nClus = 10,
#'                                              scale = FALSE),
#'                        transformList = transformList,
#'                        seed = 1)
#' FlowSOM::PlotStars(fsom)
#'
#' @importFrom dplyr '%>%' filter
#' @importFrom flowCore read.FCS transformList colnames
#' @importFrom stringr str_match
#' @importFrom pheatmap pheatmap
#' @importFrom stats density
#' @importFrom utils capture.output
#' @importFrom gridExtra grid.arrange
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
                           verbose = FALSE,
                           seed = NULL,
                           ...){

    if (verbose) message("Aggregating files ... ")

    if(!is.null(seed)) set.seed(seed)
    o <- capture.output( ff <- FlowSOM::AggregateFlowFrames(files, nCells,
                                                            channels = colsToUse,
                                                            ...))
    if(!is.null(transformList)) ff <- flowCore::transform(ff, transformList)

    FlowSOM.params <- c(FlowSOM.params,
                        list(input = ff,
                             colsToUse = colsToUse,
                             seed = seed))
    FlowSOM.params <- FlowSOM.params[unique(names(FlowSOM.params))]

    if (verbose) message("Running the FlowSOM algorithm ... ")
    fsom <- do.call(FlowSOM::FlowSOM, FlowSOM.params)

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
#' @param files       Full paths to the fcs files or a flowSet of the control
#'                    samples.
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
#'                   list(nQ = 99)). nQ is the number of quantiles
#'                   to use, 0.01 to 0.99 by default.
#' @param recompute If FALSE, will try to reuse previously saved FlowSOM model.
#'                  If so, a warning message will be printed. Default = FALSE.
#' @param ...         Additional arguments to pass to read.FCS
#'
#' @return A list containing two elements: the FlowSOM clustering and a list
#'         containing all the splines per cluster. This can be used as input
#'         for the \code{\link{CytoNorm.normalize}} function.
#' @seealso   \code{\link{QuantileNorm.train}}, \code{\link{CytoNorm.normalize}},
#'            \code{\link{prepareFlowSOM}}
#'
#' @examples
#'
#' dir <- system.file("extdata", package = "CytoNorm")
#' files <- list.files(dir, pattern = "fcs$")
#' data <- data.frame(File = files,
#'                    Path = file.path(dir, files),
#'                    Type = stringr::str_match(files, "_([12]).fcs")[,2],
#'                    Batch = stringr::str_match(files, "PTLG[0-9]*")[,1],
#'                    stringsAsFactors = FALSE)
#' data$Type <- c("1" = "Train", "2" = "Validation")[data$Type]
#' train_data <- dplyr::filter(data, Type == "Train")
#' validation_data <- dplyr::filter(data, Type == "Validation")
#'
#' ff <- flowCore::read.FCS(data$Path[1])
#' channels <- grep("Di$", flowCore::colnames(ff), value = TRUE)
#' transformList <- flowCore::transformList(channels,
#'                                          cytofTransform)
#' transformList.reverse <- flowCore::transformList(channels,
#'                                                  cytofTransform.reverse)
#'
#' model <- CytoNorm.train(files = train_data$Path,
#'                         labels = train_data$Batch,
#'                         channels = channels,
#'                         transformList = transformList,
#'                         FlowSOM.params = list(nCells = 10000, #1000000
#'                                               xdim = 10,
#'                                               ydim = 10,
#'                                               nClus = 10,
#'                                               scale = FALSE),
#'                         normParams = list(nQ = 99),
#'                         seed = 1)
#'
#' CytoNorm.normalize(model = model,
#'                    files = validation_data$Path,
#'                    labels = validation_data$Batch,
#'                    transformList = transformList,
#'                    transformList.reverse = transformList.reverse)
#'
#' @importFrom methods is
#' @importFrom flowCore sampleNames
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
                                                 nClus = 10,
                                                 scale = FALSE),
                           normMethod.train = QuantileNorm.train,
                           normParams = list(nQ = 99),
                           seed = NULL,
                           clean = TRUE,
                           plot = FALSE,
                           verbose = FALSE,
                           recompute = FALSE,
                           ...){

    if (length(labels) != length(files)) {
        stop("Input parameters 'labels' and 'files'",
             " should have the same length")
    }
  
    # Create output directory
    dirCreated = FALSE
    if (!dir.exists(outputDir)) {
        dirCreated = dir.create(outputDir)
    }

    if(recompute | !file.exists(file.path(outputDir, "CytoNorm_FlowSOM.RDS"))){

        if(!"nCells" %in% names(FlowSOM.params)){
            stop("FlowSOM.params should contain the parameter nCells.")
        }
        nCells <- FlowSOM.params[["nCells"]]
        if(is.null(c(FlowSOM.params[["channels"]], FlowSOM.params[["colsToUse"]]))){
          FlowSOM.channels <- channels
        } else if(!is.null(FlowSOM.params[["channels"]])){
          FlowSOM.channels <- FlowSOM.params[["channels"]]
        } else if(!is.null(FlowSOM.params[["colsToUse"]])){
          FlowSOM.channels <- FlowSOM.params[["colsToUse"]]
        }
        FlowSOM.params <- FlowSOM.params[grep("nCells|channels|colsToUse",
                                              names(FlowSOM.params),
                                              invert = TRUE)]
        fsom <- prepareFlowSOM(files = files,
                               nCells = nCells,
                               FlowSOM.params = FlowSOM.params,
                               transformList = transformList,
                               colsToUse = FlowSOM.channels,
                               seed = seed,
                               ...)

        saveRDS(fsom, file.path(outputDir, "CytoNorm_FlowSOM.RDS"))

        if (plot) {
            FlowSOM::FlowSOMmary(fsom,
                                 plotFile = file.path(outputDir, "CytoNorm_FlowSOM.pdf"))
        }
    } else {
        fsom <- readRDS(file.path(outputDir, "CytoNorm_FlowSOM.RDS"))
        warning("Reusing FlowSOM result previously saved at ",
    file.path(outputDir, "CytoNorm_FlowSOM.RDS. 
         If this was not intended, one can either specify another outputDir, 
         make use of the recompute parameter or move the FlowSOM object in the 
         file manager.\n"))
    }
    
    # Error when goal does not correspond to FlowSOM model
    if ("goal" %in% names(normParams) & is.list(normParams[["goal"]])){
      if (length(normParams[["goal"]]) != FlowSOM::NMetaclusters(fsom)) {
        stop(paste0(length(normParams[["goal"]]), " clusters in the goal",
                    " distribution and ", FlowSOM::NMetaclusters(fsom), 
                    " clusters in the FlowSOM model.
  This should be the same."))
      }
      if ("nQ" %in% names(normParams)){
        if (nrow(normParams[["goal"]][[1]]) != normParams[["nQ"]]) {
          stop(paste0(nrow(normParams[["goal"]][[1]]), " quantiles in the goal",
                      " distribution and ", normParams[["nQ"]], 
                      " quantiles in the CytoNorm call.
  This should be the same."))
        }
      }
    }

    # Split files by clusters
    for(i in seq_along(files)) {
      if(is(files, "flowSet")) {
        file <- sampleNames(files)[i]
        ff <- files[[i]]
      } else {
        file <- files[i]
        ff <- flowCore::read.FCS(file, ...)
      }
      if(verbose) message("Splitting ", file)
      if (!is.null(transformList)) {
        ff <- flowCore::transform(ff,transformList)
      }
      # Map the file to the FlowSOM clustering
      fsom_file <- FlowSOM::NewData(fsom, ff)
      # Get the metacluster label for every cell
      cellClusterIDs <- FlowSOM::GetMetaclusters(fsom_file) #fsom$metaclustering[GetClusters(fsom_file)]
      for (cluster in unique(fsom$metaclustering)) {
        if (sum(cellClusterIDs == cluster) > 0) {
          suppressWarnings(
            flowCore::write.FCS(
              ff[cellClusterIDs == cluster,],
              file=file.path(outputDir,
                             paste0(gsub("[:/]","_",file),
                                    "_fsom",cluster,".fcs"))))
        }
      }
    }

    # file names
    if(is(files, "flowSet")) {
      file_names <- sampleNames(files)
    } else {
      file_names <- files
    }

    # Learn quantiles for each cluster
    clusterRes <- list()
    for (cluster in unique(fsom$metaclustering)) {
        if(verbose) message("Processing cluster ",cluster)
        if (plot) {
            grDevices::pdf(file.path(outputDir,
                                     paste0("CytoNorm_norm_Cluster",
                                            cluster, ".pdf")),
                           height = 3*(2*length(files)+2),
                           width = 3*(length(channels)+1))
        }

        normParams_tmp <- c(normParams,
                            list(files = file.path(outputDir,
                                                   paste0(gsub("[:/]", "_", file_names),
                                                          "_fsom", cluster, ".fcs")),
                                 labels = as.character(labels),
                                 channels = channels,
                                 transformList = NULL,
                                 verbose = verbose,
                                 plot = plot))
        normParams_tmp <- normParams_tmp[unique(names(normParams_tmp))]
        if(is.list(normParams[["goal"]])){
            normParams_tmp[["goal"]] <- normParams[["goal"]][[cluster]]
        }
        clusterRes[[cluster]] <- do.call(normMethod.train,
                                         normParams_tmp)

        if (plot) { grDevices::dev.off() }
    }

    if(clean){
        for(cluster in unique(fsom$metaclustering)){
            tmp_files <- file.path(outputDir,
                                   paste0(gsub("[:/]", "_", file_names),
                                          "_fsom", cluster, ".fcs"))

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
#' @param files       Full paths of the fcs files or a flowSet of the samples.
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
#' @param write       logical indicating whether the normalised samples should
#'                    be written to new FCS files in a \code{Normalized}
#'                    directory within \code{outputDir}, set to TRUE by default.
#' @param ...         Additional arguments to pass to read.FCS
#' @return a flowSet containing the normalised samples and optionally write FCS
#'                    files to \code{Normalized} directory in \code{outputDir}.
#' @seealso   \code{\link{CytoNorm.train}}
#'
#' @examples
#'
#'
#' dir <- system.file("extdata", package = "CytoNorm")
#' files <- list.files(dir, pattern = "fcs$")
#' data <- data.frame(File = files,
#'                    Path = file.path(dir, files),
#'                    Type = stringr::str_match(files, "_([12]).fcs")[,2],
#'                    Batch = stringr::str_match(files, "PTLG[0-9]*")[,1],
#'                    stringsAsFactors = FALSE)
#' data$Type <- c("1" = "Train", "2" = "Validation")[data$Type]
#' train_data <- dplyr::filter(data, Type == "Train")
#' validation_data <- dplyr::filter(data, Type == "Validation")
#'
#' ff <- flowCore::read.FCS(data$Path[1])
#' channels <- grep("Di$", flowCore::colnames(ff), value = TRUE)
#' transformList <- flowCore::transformList(channels,
#'                                          cytofTransform)
#' transformList.reverse <- flowCore::transformList(channels,
#'                                                  cytofTransform.reverse)
#'
#' model <- CytoNorm.train(files = train_data$Path,
#'                         labels = train_data$Batch,
#'                         channels = channels,
#'                         transformList = transformList,
#'                         FlowSOM.params = list(nCells = 10000, #1000000
#'                                               xdim = 15,
#'                                               ydim = 15,
#'                                               nClus = 10,
#'                                               scale = FALSE),
#'                         normParams = list(nQ = 99),
#'                         seed = 1,
#'                         verbose = TRUE)
#'
#' CytoNorm.normalize(model = model,
#'                    files = validation_data$Path,
#'                    labels = validation_data$Batch,
#'                    transformList = transformList,
#'                    transformList.reverse = transformList.reverse,
#'                    verbose = TRUE)
#'
#' @importFrom methods is
#' @importFrom flowCore sampleNames flowSet
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
                                normMethod.normalize = QuantileNorm.normalize,
                               write = TRUE,
                               ...){
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

    # file names
    if(is(files, "flowSet")) {
      file_names <- sampleNames(files)
    } else {
      file_names <- files
    }

    # Split files by clusters
    cellClusterIDs <- list()
    meta <- list()
    cluster_files <- list()
    for(i in seq_along(files)) {
      if(is(files, "flowSet")) {
        file <- file_names[i]
        ff <- files[[i]]
      } else {
        file <- files[i]
        ff <- flowCore::read.FCS(file, ...)
      }
      if(verbose) message("Splitting ",file)
      if(!is.null(transformList)){
        ff <- flowCore::transform(ff, transformList)
        # meta[[file]] <- list()
        # meta[[file]][["description_original"]] <- ff@description
        # meta[[file]][["parameters_original"]] <- ff@parameters
      }

      fsom_file <- FlowSOM::NewData(fsom,ff)

      cellClusterIDs[[file]] <- FlowSOM::GetMetaclusters(fsom_file)

      for(cluster in unique(fsom$metaclustering)){
        if (sum(cellClusterIDs[[file]] == cluster) > 0) {
          f <- file.path(outputDir,
                         paste0(gsub("[:/]","_",file),
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
                               paste0(gsub("[:/]",
                                           "_",
                                           file_names),
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
    res <- lapply(
        seq_along(files),
        function(i) {

            if(is(files, "flowSet")) {
                file <- file_names[i]
                ff <- files[[i]]
            } else {
                file <- files[i]
                ff <- flowCore::read.FCS(file, ...)
            }
            if(!is.null(transformList)){
              ff <- flowCore::transform(ff, transformList)
            }
            if(verbose) message("Rebuilding ",file)
            for(cluster in unique(fsom$metaclustering)){
                file_name <- file.path(outputDir,
                                       paste0("Norm_",gsub("[:/]","_",file),
                                              "_fsom",cluster,".fcs"))
                if (file.exists(file_name)) {
                    ff_subset <- flowCore::read.FCS(file_name, ...)
                    flowCore::exprs(ff)[cellClusterIDs[[file]] == cluster,] <- flowCore::exprs(ff_subset)
                }
            }
            if(!is.null(transformList.reverse)){
                ff <- flowCore::transform(ff, transformList.reverse)
                # ff@description <- meta[[file]][["description_original"]]
                # ff@parameters <- meta[[file]][["parameters_original"]]
            }


            # Adapt to real min and max because this gets strange values otherwise

                # params <- flowCore::parameters(fcs)
                # pd <-NULL
                # cols <- as.vector(pd$name)
                # idxs <- match(cols, pd$name)
                # if (any(is.na(idxs))) {
                #     stop("Invalid column specifier")
                # }
                # keyval <- list()
                # for (channel_number in 1:ncol(exprs)){
                #     channel_name<-colnames(exprs)[channel_number]
                #     if (is.null(desc1))
                #         desc1<-colnames(exprs)[channel_number]
                #     channel_id <- paste("$P", channel_number, sep = "")
                #     channel_range <- max(exprs[,channel_number]) + 1
                #     channel_min<-min(0,min(exprs[,channel_number])-1)
                #     plist <- matrix(c(channel_name, desc1[channel_number], channel_range, channel_min, channel_range - 1))
                #     rownames(plist) <- c("name", "desc", "range", "minRange", "maxRange")
                #     colnames(plist) <- c(channel_id)
                #     pd <- rbind(pd, t(plist))
                #     keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
                #     keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
                #     keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
                #     keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
                #     keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
                # }




            ff@parameters@data[,"minRange"] <- floor(pmin(0, apply(ff@exprs, 2, min)))
            ff@parameters@data[,"maxRange"] <- ceiling(apply(ff@exprs, 2, max))
            ff@parameters@data[,"range"] <- pmax(ff@parameters@data[,"maxRange"] + 1 ,
                                                 ff@parameters@data[,"maxRange"] - ff@parameters@data[,"minRange"])
            for(i in seq_len(flowCore::ncol(ff))){
                ff@description[[paste0("$P",i,"R")]] <- ff@parameters@data[i,"range"]
            }

            if(clean){
                file.remove(file.path(outputDir,
                                      paste0("Norm_",gsub("[:/]","_",file),
                                             "_fsom",unique(fsom$metaclustering),".fcs")))
            }

            if(write) {
                suppressWarnings(
                    flowCore::write.FCS(
                        ff,
                        file= file.path(outputDir,paste0(prefix,gsub(".*/","",file)))
                    )
                )
            }
            return(ff)
        })

    # dremove empty output directory
    if(length(list.files(outputDir)) == 0){
        unlink(outputDir)
    }

    # normalized flowSet
    if(is(files, "flowSet")) {
        res <- flowSet(res)
        return(res)
    }
}


#' Extract quantiles used by a previous CytoNorm model
#'
#' Can be of interest to pass as goal quantiles to a new model.
#'
#' @param model Model trained with CytoNorm.train
#' @param type "ref" for goal quantiles or batch label
#'
#' @examples
#' dir <- system.file("extdata", package = "CytoNorm")
#' files <- list.files(dir, pattern = "fcs$")
#' data <- data.frame(File = files,
#'                    Path = file.path(dir, files),
#'                    Type = stringr::str_match(files, "_([12]).fcs")[,2],
#'                    Batch = stringr::str_match(files, "PTLG[0-9]*")[,1],
#'                    stringsAsFactors = FALSE)
#' data$Type <- c("1" = "Train", "2" = "Validation")[data$Type]
#' train_data <- dplyr::filter(data, Type == "Train")
#' validation_data <- dplyr::filter(data, Type == "Validation")
#'
#' ff <- flowCore::read.FCS(data$Path[1])
#' channels <- grep("Di$", flowCore::colnames(ff), value = TRUE)
#' transformList <- flowCore::transformList(channels,
#'                                          cytofTransform)
#' transformList.reverse <- flowCore::transformList(channels,
#'                                                  cytofTransform.reverse)
#'
#' model <- CytoNorm.train(files = train_data$Path,
#'                         labels = train_data$Batch,
#'                         channels = channels,
#'                         transformList = transformList,
#'                         FlowSOM.params = list(nCells = 10000, #1000000
#'                                               xdim = 10,
#'                                               ydim = 10,
#'                                               nClus = 10,
#'                                               scale = FALSE),
#'                         normParams = list(nQ = 99),
#'                         seed = 1)
#'
#' quantiles <- getCytoNormQuantiles(model)
#'
#' @export
getCytoNormQuantiles <- function(model, type = "ref"){
    if(type == "ref"){
        lapply(model$clusterRes, function(x) x$refQuantiles)
    } else {
        lapply(model$clusterRes, function(x) x$quantiles[[type]])
    }
}

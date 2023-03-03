#' getQuantiles
#'
#' @param files       Full paths of to the fcs files of the samples
#' @param channels    Names of the channels to compute the quantiles for
#' @param transformList   Transformation list to pass to the flowCore
#'                    \code{transform} function
#' @param nQ          Number of quantiles to compute Default = 101, which
#'                    results in quantiles for every percent of the data.
#'                    Ignored if quantileValues is given.
#' @param minCells    Minimum number of cells required to compute trust-worthy
#'                    quantiles. Otherwise NA is returned. Default = 50.
#' @param quantileValues Vector of length with values between 0 and 1, giving
#'                        the percentages at which the quantiles should be
#'                        computed. If NULL (default), the quantiles will be
#'                        evenly distributed, including 0 and 1.
#' @param labels      A label for every file, indicating to which group it
#'                    belongs. If multiple files have the same label, they
#'                    get aggregated. If NULL, all files are handled separately.
#' @param selection   List with indexation vector for every file.
#' @param verbose     If TRUE, extra output is printed. Default = FALSE
#' @param plot        If TRUE, plots are generated showing all quantiles.
#'                    Default = FALSE.
#' @param ...         Additional arguments to pass to read.FCS
#'
#' @examples
#' dir <- system.file("extdata", package = "CytoNorm")
#' files <- list.files(dir, pattern = "fcs$")
#'
#' ff <- flowCore::read.FCS(file.path(dir, files[1]))
#' channels <- grep("Di$", flowCore::colnames(ff), value = TRUE)
#' transformList <- flowCore::transformList(channels,
#'                                          cytofTransform)
#'
#' quantiles <- getQuantiles(files = file.path(dir, files),
#'                           channels = channels,
#'                           transformList = transformList)
#'
#' pheatmap::pheatmap(quantiles[[1]],
#'                    cluster_rows = FALSE,
#'                    cluster_cols = FALSE,
#'                    labels_col =
#'                      paste0(FlowSOM::GetMarkers(ff, colnames(quantiles[[1]])),
#'                             " (", colnames(quantiles[[1]]), ")"),
#'                    main = files[1])
#'
#' @export
getQuantiles <- function(files,
                         channels,
                         nQ = 101,
                         minCells = 50,
                         quantileValues = NULL,
                         transformList = NULL,
                         labels = NULL,
                         selection = NULL,
                         verbose = FALSE,
                         plot = FALSE,
                         ...){

    if (is.null(labels)) labels <- files

    # Compute quantiles for each label
    if (verbose) message("Computing Quantiles")

    quantiles <- list()

    if (is.null(quantileValues)) {
        quantileValues <- c(0, (1:(nQ-1))/(nQ-1))
    } else {
        nQ <- length(quantileValues)
    }

    for(label in unique(labels)){
        ids <- which(labels == label & file.exists(files))
        if(verbose) message("  ", label, " (FileID ", paste(ids, collapse = "," ), ")")

        # Read the file(s) and transform if necessary
        if (length(ids) > 1) {
            ff <- FlowSOM::AggregateFlowFrames(files[ids], 1e12,
                                               keepOrder = TRUE,
                                               channels = channels,
                                               ...)
        } else if(length(ids) == 1) {
            o <- capture.output(ff <- flowCore::read.FCS(files[ids],...))
            if (verbose) message(o)
        } else {
            ff <- NULL
        }

        if (!is.null(ff) && !is.null(transformList)) {
            ff <- flowCore::transform(ff, transformList)
        }

        if (!is.null(ff) & !is.null(selection)){
            ff <- ff[selection[[file]], ]
        }


        # Compute quantiles for all channels to normalize
        if (!is.null(ff) && flowCore::nrow(ff) > minCells) {
            quantiles[[label]] <- apply(flowCore::exprs(ff)[, channels, drop = FALSE],
                                        2,
                                        function(x){
                                            stats::quantile(x,
                                                            quantileValues)
                                        })

            if(plot){
                textPlot(label)
                for(channel in channels){
                    dens <- stats::density(flowCore::exprs(ff)[, channel],
                                           bw = 0.1)
                    graphics::plot(dens,
                                   bty = "n", xaxt = "n", yaxt = "n",
                                   xlab = "", ylab = "", main = "",
                                   xlim = c(0,
                                            max(flowCore::exprs(ff)[, channel],
                                                8)))
                    graphics::abline(v = quantiles[[label]][, channel],
                                     col = "grey")
                    graphics::lines(dens, lwd = 2)
                }
            }
        } else {
            if(is.null(ff)){
                message("  Could not find ", paste(files[ids], collapse = ", "))
            } else {
                message("  Less then ", minCells, " cells in ", label,
                        " (", flowCore::nrow(ff), "). No quantiles computed.")
            }
            quantiles[[label]] <- matrix(NA,
                                         nrow = nQ,
                                         ncol = length(channels),
                                         dimnames =
                                             list(as.character(quantileValues),
                                                  channels))
            if(plot){
                textPlot(label)
                for(channel in channels){
                   textPlot("NA")
                }
            }
        }
    }

    return(quantiles)
}

#' QuantileNorm.train
#'
#' Learn the batch effects from control samples.
#' This function computes quantiles to describe the distribution of the data,
#' and infers spline functions to equalize these distributions over the files.
#' Typically, you will use the function \code{\link{QuantileNorm.normalize}}
#' after this function to reverse the batch effects on other files.
#'
#' @param files       Full paths of to the fcs files of the control samples
#' @param labels      A label for every file, indicating to which batch it
#'                    belongs, e.g. the plate ID.
#' @param channels    Names of the channels to normalize
#' @param transformList   Transformation list to pass to the flowCore
#'                    \code{transform} function
#' @param nQ          Number of quantiles to use. Default = 101, which results in
#'                    quantiles for every percent of the data.
#' @param quantileValues If specified, it should be a vector of length nQ with
#'                        values between 0 and 1, giving the percentages at
#'                        which the quantiles should be computed. If NULL
#'                        (default), the quantiles will be evenly distributed,
#'                        including 0 and 1.
#' @param limit       These values will be modelled to map onto themselves by the
#'                    spline
#' @param goal        Goal distribution. Default "mean", can also be nQ numeric
#'                    values or one of the batch labels.
#' @param plot        If TRUE, a plot is generated (using the \code{layout}
#'                    function) showing all quantiles. Default = FALSE.
#' @param plotTitle   Title to use in the plot. Default = "Quantiles".
#' @param verbose     If TRUE, progress updates are printed. Default = FALSE.
#' @param ...         Additional arguments to pass to read.FCS
#'
#' @return A list containing all the splines and quantile information. This can
#'         be used as input for the \code{\link{QuantileNorm.normalize}} function.
#' @seealso   \code{\link{CytoNorm.train}}, \code{\link{QuantileNorm.normalize}}
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
#' png("nQ101.png",
#'     width = length(channels) * 300,
#'     height = (nrow(train_data) * 2 + 1) * 300)
#' model_nQ_101 <- QuantileNorm.train(
#'   files = train_data$Path,
#'   labels = train_data$Batch,
#'   channels = channels,
#'   transformList = transformList,
#'   nQ = 101,
#'   plot = TRUE)
#' dev.off()
#'
#' png("nQ101_limited.png",
#'     width = length(channels) * 300,
#'     height = (nrow(train_data) * 2 + 1) * 300)
#' model_nQ_101 <- QuantileNorm.train(
#'   files = train_data$Path,
#'   labels = train_data$Batch,
#'   channels = channels,
#'   transformList = transformList,
#'   nQ = 101,
#'   limit = c(0,8),
#'   plot = TRUE)
#' dev.off()
#'
#' png("nQ_2.png",
#'     width = length(channels) * 300,
#'     height = (nrow(train_data) * 2 + 1) * 300)
#' model_nQ_2 <- QuantileNorm.train(
#'   files = train_data$Path,
#'   labels = train_data$Batch,
#'   channels = channels,
#'   transformList = transformList,
#'   nQ = 2,
#'   quantileValues = c(0.001, 0.999),
#'   plot = TRUE)
#' dev.off()
#'
#' model_goal_mean <- QuantileNorm.train(
#'   files = train_data$Path,
#'   labels = train_data$Batch,
#'   channels = channels,
#'   transformList = transformList)
#'
#' model_goal_batch1 <- QuantileNorm.train(
#'   files = train_data$Path,
#'   labels = train_data$Batch,
#'   channels = channels,
#'   transformList = transformList,
#'   goal = "PTLG021")
#'
#' model_goal_fixed <- QuantileNorm.train(
#'   files = train_data$Path,
#'   labels = train_data$Batch,
#'   channels = channels,
#'   transformList = transformList,
#'   nQ = 21,
#'   goal = seq(0, 1, by = 0.05))
#'
#' @export
QuantileNorm.train <- function(files,
                               labels,
                               channels,
                               transformList,
                               nQ = 101,
                               limit = NULL,
                               quantileValues = NULL,
                               goal = "mean",
                               verbose = FALSE,
                               plot = FALSE,
                               plotTitle = "Quantiles",
                               ...){

    if(length(labels) != length(files)){
        stop("Input parameters 'labels' and 'files'",
             " should have the same length")
    }

    labels <- as.character(labels)

    if(plot){
        xdim = 1 + length(channels)
        ydim = 2 + 2*length(unique(labels))
        graphics::layout(matrix(1:(xdim*ydim), ncol = xdim, byrow = TRUE))
        graphics::par(mar=c(1, 1, 0, 0))
        textPlot(plotTitle)
        ff_tmp <- flowCore::read.FCS(files[file.exists(files)][1],...)

        for(channel in channels){
            marker <- flowCore::getChannelMarker(ff_tmp, channel)[,2]
            if(is.na(marker)){
                textPlot(channel)
            } else {
                textPlot(paste0(marker, " (", channel, ")"))
            }
        }
    }

    quantiles <- getQuantiles(files = files,
                              labels = labels,
                              channels = channels,
                              transformList = transformList,
                              nQ = nQ,
                              quantileValues = quantileValues,
                              verbose = verbose,
                              plot = plot,
                              ...)

    if(!is.null(limit)){
        quantiles <- lapply(quantiles, function(quantile_matrix){
            rbind(quantile_matrix,
                  matrix(limit,
                         nrow = length(limit),
                         ncol = ncol(quantile_matrix),
                         byrow = FALSE,
                         dimnames = list(paste0("Fixed_",limit), NULL)))
        })
    }

    # Get the goal distributions
    if(is.character(goal) && goal == "mean"){
        refQuantiles <- matrix(apply(matrix(unlist(quantiles),
                                            ncol = length(unique(labels))),
                                     1, mean, na.rm = TRUE),
                               nrow = nQ+length(limit),
                               dimnames = list(rownames(quantiles[[1]]),
                                               channels))
    } else if (is.numeric(goal)) {
        if(length(goal) == nQ){
            refQuantiles <- matrix(goal,
                                   nrow = nQ,
                                   ncol = length(channels),
                                   dimnames = list(quantileValues,
                                                   channels))
        } else if ((nrow(goal) == nQ) & (ncol(goal) == length(channels))){
            refQuantiles <- goal
        } else {
            stop("Goal should be 'mean', a batch label, ",
                   "a numeric vector of length nQ, or a matrix with nQ rows
                   and length(channels) columns.")
        }
    } else if (goal %in% unique(labels)) {
        refQuantiles <- quantiles[[goal]]
    } else {
        stop("Goal should be 'mean', a batch label, ",
             "a numeric vector of length nQ, or a matrix with nQ rows
             and length(channels) columns.")
    }

    if(plot){
        textPlot("Goal distribution")
        for(channel in channels){
            graphics::plot(0, type = "n", xlim = c(0, 8),
                           bty = "n", xaxt = "n", yaxt = "n",
                           xlab = "", ylab = "", main = "")
            graphics::abline(v = refQuantiles[,channel])
        }
    }

    # Compute splines for each label
    splines <- list()
    if(verbose) message("Computing Splines")

    for(label in unique(labels)){
        if(verbose) message("  ",label)
        if(plot){ textPlot(label) }
        splines[[label]] <- list()
        if(!is.null(quantiles[[label]]) & !any(is.na(quantiles[[label]])) & !any(is.na(refQuantiles))){
            for(channel in channels){

                refQ <- refQuantiles[, channel]
                labelQ <- quantiles[[label]][, channel]

                if(length(unique(labelQ)) > 1){
                    suppressWarnings(spl <- stats::splinefun(labelQ,
                                                             refQ,
                                                             method="monoH.FC"))
                } else {
                    spl <- identityFunction
                    warning("Not enough unique quantiles  for ", label, " in ",
                            channel, ". The identity function will be used.")
                }
                splines[[label]][[channel]] <- spl

                if(plot){
                    graphics::plot(labelQ, refQ, xlim = c(0, 8), ylim = c(0, 8),
                                   pch = 19, bty = "n", xaxt = "n", yaxt = "n",
                                   xlab = "", ylab = "", main = "")
                    graphics::lines(c(-0.5, 8), c(-0.5, 8), col="#999999")
                    x <- seq(-0.5, 8, 0.1)
                    graphics::lines(x,
                                    splines[[label]][[channel]](x),
                                    col = "#b30000")
                }

            }
        } else {
            warning("Not enough cells for ", label,
                    "\nThe identity function will be used.")
            for(channel in channels){
                splines[[label]][[channel]] <- identityFunction
                if(plot){
                    graphics::plot(c(0,8), c(0,8), col = "#999999", type = "l",
                                   xlim = c(0,8), ylim = c(0,8),
                                   pch = 19, bty = "n", xaxt = "n", yaxt = "n",
                                   xlab = "", ylab = "", main = "")
                }
            }
        }
    }
    named.list(channels, splines, quantiles, quantileValues, refQuantiles)
}


#' QuantileNorm.normalize
#'
#' Normalize data, given the batch effects learned from control samples per
#' cell type/cluster (output from \code{\link{QuantileNorm.train}}). New fcs
#' files are written to the given output directory.
#'
#' @param model       Model of the batch effercts, as computed by
#'                    \code{\link{QuantileNorm.train}}
#' @param files       Full paths of to the fcs files of the samples to normalize
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
#' @param removeOriginal Should the original fcs be removed? Default = FALSE.
#' @param verbose     If TRUE, progress updates are printed. Default = FALSE.
#' @param ...         Additional arguments to pass to read.FCS
#'
#' @return Nothing is returned, but the new FCS files are written to the output
#'         directory
#' @seealso \code{\link{QuantileNorm.train}}
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
#' model_nQ_101 <- QuantileNorm.train(
#'   files = train_data$Path,
#'   labels = train_data$Batch,
#'   channels = channels,
#'   transformList = transformList,
#'   nQ = 101)
#'
#' QuantileNorm.normalize(model_nQ_101,
#'                        validation_data$Path,
#'                        validation_data$Batch,
#'                        transformList = transformList,
#'                        transformList.reverse = transformList.reverse)
#' @export
QuantileNorm.normalize <- function(model,
                                   files,
                                   labels,
                                   transformList,
                                   transformList.reverse,
                                   outputDir = ".",
                                   prefix = "Norm_",
                                   removeOriginal = FALSE,
                                   verbose = FALSE,
                                   ...){

    if(is.null(model$channels) |
       is.null(model$splines) |
       is.null(model$quantiles)){
        stop("The 'model' parameter should be the result of using the ",
             "QuantileNorm.train function.")
    }
    if(length(labels) != length(files)){
        stop("Input parameters 'labels' and 'files' ",
             "should have the same length")
    }

    # Create temporary directory
    if (!dir.exists(outputDir)) dir.create(outputDir)

    labels <- labels
    channels <- model$channels

    # Normalize each file based on label
    for(i in seq_along(files)){
        file <- files[i]

        if(file.exists(file)){
            label <- as.character(labels[i])
            if(verbose) message("  ",file," (",label,")")

            if(label %in% names(model$splines)){
                # Read the file
                ff <- flowCore::read.FCS(file,...)

                # Transform if necessary
                if(!is.null(transformList)){
                    #description_original <- ff@description
                    #parameters_original <- ff@parameters
                    ff <- flowCore::transform(ff, transformList)
                }

                # Overwrite the values with the normalized values
                if (verbose) message("Normalizing ",label)
                for (channel in channels) {
                    flowCore::exprs(ff)[, channel] <-
                        model$splines[[label]][[channel]](
                            flowCore::exprs(ff[, channel]))
                }

                if (!is.null(transformList.reverse)) {
                    ff <- flowCore::transform(ff, transformList.reverse)
                } else if (!is.null(transformList)) {
                    warning("Please provide a reverse transformation list if
                            you want the files to be saved in the original
                            untransformed space.")
                }

                if (removeOriginal) {
                    file.remove(file)
                }

                suppressWarnings(flowCore::write.FCS(ff,
                                                     filename = file.path(outputDir,
                                                                          paste0(prefix,
                                                                                 gsub(".*/","",file)))))

            } else {
                warning("The model was not trained for ", label, ".")
            }
        } else {
            warning(file, " does not exist, skipped.")
        }
    }
}

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
#' @param nQ          Number of quantiles to use. Default = 21, which results in
#'                    quantiles for every 5 percent of the data.
#' @param spar        Smoothing parameter, typically between (0,1]. Default=0.5
#'                    See \code{smooth.spline} from the \code{stats} package
#' @param plot        If TRUE, a plot is generated (using the \code{layout}
#'                    function) showing all quantiles. Default = FALSE.
#' @param plotTitle   Title to use in the plot. Default = "Quantiles".
#' @param verbose     If TRUE, progress updates are printed. Default = FALSE.
#'
#' @return A list containing all the splines and quantile information. This can
#'         be used as input for the \code{\link{QuantileNorm.normalize}} function.
#' @seealso   \code{\link{CytofNorm.train}}, \code{\link{QuantileNorm.normalize}}
#'
#' @examples
#'
#' @export
QuantileNorm.train <- function(files,
                               labels,
                               channels,
                               transformList,
                               nQ = 21,
                               spar = 0.5,
                               verbose = FALSE,
                               plot = FALSE,
                               plotTitle = "Quantiles"){

    if(length(labels) != length(files)){
        stop("Input parameters 'labels' and 'files'",
             " should have the same length")
    }

    labels <- as.character(labels)

    # Compute quantiles for each label
    if(verbose) message("Computing Quantiles")
    quantiles <- list()

    if(plot){
        xdim = 1 + length(channels)
        ydim = 2 + 2*length(unique(labels))
        graphics::layout(matrix(1:(xdim*ydim), ncol = xdim, byrow = TRUE))
        graphics::par(mar=c(1, 1, 0, 0))
        textPlot(plotTitle)
        ff_tmp <- flowCore::read.FCS(files[1])

        for(channel in channels){
            marker <- flowCore::getChannelMarker(ff_tmp, channel)[,2]
            if(is.na(marker)){
                textPlot(channel)
            } else {
                textPlot(paste0(marker, " (", channel, ")"))
            }
        }
    }

    for(label in unique(labels)){
        ids <- which(labels == label & file.exists(files))
        if(verbose) message("  ", label, "(", paste(ids, collapse = "," ), ")")

        # Read the file(s) and transform if necessary
        if (length(ids) > 1) {
            ff <- FlowSOM::AggregateFlowFrames(files[ids], 1e12, keepOrder = TRUE)
        } else if(length(ids) == 1) {
            ff <- flowCore::read.FCS(files[ids])
        } else {
            ff <- NULL
        }

        if (!is.null(ff) && !is.null(transformList)) {
            ff <- flowCore::transform(ff, transformList)
        }

        # Compute quantiles for all channels to normalize
        if(!is.null(ff) && flowCore::nrow(ff) > nQ){
            quantiles[[label]] <- apply(flowCore::exprs(ff)[, channels],
                                            2,
                                            function(x){
                                                stats::quantile(x,
                                                                c(0,(1:(nQ-1))/
                                                                      (nQ-1)))
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
                                  max(flowCore::exprs(ff)[, channel], 7)))
                    graphics::abline(v = quantiles[[label]][, channel],
                                     col = "grey")
                    graphics::lines(dens, lwd = 2)
                }
            }
        } else {
            warning("Not enough cells in file ",files[ids],
                    "\nThe identity function will be used.")
            quantiles[[label]] <- matrix(NA,
                                         nrow = nQ,
                                         ncol = length(channels),
                                         dimnames = list(c(0,(1:(nQ-1))/(nQ-1)),
                                                    channels))
        }
    }

    # Get the goal distributions
    refQuantiles <- matrix(apply(matrix(unlist(quantiles),
                                        ncol = length(unique(labels))),
                                 1, mean, na.rm = TRUE),
                           nrow = nQ,
                           dimnames = list(c(0,(1:(nQ-1))/(nQ-1)),
                                           channels))

    if(plot){
        textPlot("Goal distribution")
        for(channel in channels){
            graphics::plot(0, type = "n", xlim = c(0, 7),
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
        if(!is.null(quantiles[[label]]) | any(is.na(quantiles[[label]]))){
            for(channel in channels){

                refQ <- refQuantiles[, channel]
                labelQ <- quantiles[[label]][, channel]

                # Remove leading zeros and NA values to avoid computational problems

                leadingZeros <- max(sum(cumsum(refQ)<1e-5,na.rm=TRUE) -1,
                                    sum(cumsum(labelQ)<1e-5,na.rm=TRUE)-1)
                if(leadingZeros > 0){
                    leadingZeros <- seq(leadingZeros)
                } else {
                    leadingZeros <- c()
                }
                naValues <- which(is.na(labelQ)|is.na(refQ))
                toRemove <- c(leadingZeros,naValues)

                if(length(toRemove) > 0){
                    labelQ <- labelQ[-toRemove]
                    refQ <- refQ[-toRemove]
                }

                if(length(labelQ) > 3){
                    spl <- stats::smooth.spline(labelQ, refQ, spar = spar)
                } else {
                    spl <- stats::smooth.spline(1:5, 1:5, spar=spar)
                    warning("Too many zeros or missing values in the quantiles",
                            " for ", label, " in ", channel)
                }
                splines[[label]][[channel]] <- spl

                if(plot){
                    graphics::plot(labelQ, refQ, xlim = c(0, 7), ylim = c(0, 7),
                                   pch = 19, bty = "n", xaxt = "n", yaxt = "n",
                                   xlab = "", ylab = "", main = "")
                    graphics::lines(c(0, 7), c(0, 7), col="#999999")
                    graphics::lines(seq(0, 7, 0.1),
                                    stats::predict(splines[[label]][[channel]],
                                            seq(0, 7, 0.1))$y)
                }

            }
        } else {
            for(channel in channels){
                splines[[label]][[channel]] <- stats::splinefun(1:5, 1:5,
                                                                method="hyman")
                if(plot){
                    graphics::plot(c(0, 7), c(0, 7), col = "#999999", type = "l",
                                   xlim = c(0, 7), ylim = c(0, 7),
                                   pch = 19, bty = "n", xaxt = "n", yaxt = "n",
                                   xlab = "", ylab = "", main = "")
                }
            }
        }
    }
    named.list(channels, splines, quantiles, refQuantiles)
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
#'
#' @return Nothing is returned, but the new FCS files are written to the output
#'         directory
#' @seealso \code{\link{QuantileNorm.train}}
#'
#' @examples
#'
#' @export
QuantileNorm.normalize <- function(model,
                                   files,
                                   labels,
                                   transformList,
                                   transformList.reverse,
                                   outputDir = ".",
                                   prefix = "Norm_",
                                   removeOriginal = FALSE,
                                   verbose = FALSE){

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
                ff <- flowCore::read.FCS(file)

                # Transform if necessary
                if(!is.null(transformList)){
                    #description_original <- ff@description
                    #parameters_original <- ff@parameters
                    ff <- flowCore::transform(ff, transformList)
                }

                # Overwrite the values with the normalized values
                if (verbose) message("Normalizing ",label)
                for (channel in channels) {
                    flowCore::exprs(ff)[, channel] <- stats::predict(
                        model$splines[[label]][[channel]],
                        flowCore::exprs(ff[, channel]))$y
                }

                if (!is.null(transformList.reverse)) {
                    ff <- flowCore::transform(ff, transformList.reverse)
                } else if (!is.null(transformList)) {
                    warning("Please provide a reverse transformation list if
                            you want the files to be saved in the original
                            untransformed space.")
                }

                suppressWarnings(flowCore::write.FCS(ff,
                                                     filename = file.path(outputDir,
                                                                          paste0(prefix,
                                                                                 gsub(".*/","",file)))))
                if (removeOriginal) {
                    file.remove(file)
                }
            } else {
                warning("The model was not trained for ", label, ".")
            }
        } else {
            warning(file, " does not exist, skipped.")
        }
    }
}

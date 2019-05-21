#' MinMaxNorm.train
#'
#' Learn the batch effects from control samples.
#' This function computes min and max to describe the distribution of the data,
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
#' @param minmax_q    Quantiles to approximate minimum and maximum. Use c(0,1)
#'                    for exact minimum and maximum. Default = c(0.001, 0.999)
#'
#' @return A list containing the min max information. This can
#'         be used as input for the \code{\link{MinMaxNorm.normalize}} function.
#' @seealso   \code{\link{MinMaxNorm.normalize}}
#'
#' @examples
#'
#' @export
MinMaxNorm.train <- function(files,
                             labels,
                             channels,
                             transformList,
                             minmax_q = c(0.001, 0.999),
                             verbose = FALSE,
                             plot = FALSE,
                             plotTitle = "Min Max"){

    if(length(labels) != length(files)){
        stop("Input parameters 'labels' and 'files'",
             " should have the same length")
    }

    if(length(minmax_q) != 2){
        stop("minmax_q should contain two numeric values.")
    }

    labels <- as.character(labels)

    # Compute quantiles for each label
    if(verbose) message("Computing min max values")
    minmax <- list()

    for(label in unique(labels)){
        if(verbose) message("  ",label)

        # Read the file(s) and transform if necessary
        ids <- which(labels == label)
        if(length(ids) > 1){
            ff <- FlowSOM::AggregateFlowFrames(files[ids], 1e12, keepOrder = TRUE)
        } else {
            ff <- flowCore::read.FCS(files[ids])
        }
        if(!is.null(transformList)) ff <- flowCore::transform(ff, transformList)

        minmax[[label]] <- apply(flowCore::exprs(ff)[, channels],
                                 2,
                                 function(x){
                                     stats::quantile(x,
                                                     minmax_q)
                                 })
    }

    minmax[["goal"]] <- Reduce("+", minmax) / length(minmax)

    named.list(channels, minmax)
}


#' MinMaxNorm.normalize
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
#' @seealso \code{\link{MinMaxNorm.train}}
#'
#' @examples
#'
#' @export
MinMaxNorm.normalize <- function(model,
                                 files,
                                 labels,
                                 transformList,
                                 transformList.reverse,
                                 outputDir = ".",
                                 prefix = "Norm_",
                                 removeOriginal = FALSE,
                                 verbose = FALSE){

    if(is.null(model$minmax)){
        stop("The 'model' parameter should be the result of using the ",
             "MinMaxNorm.train function.")
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
        label <- as.character(labels[i])
        if(verbose) message("  ",file," (",label,")")

        if(label %in% names(model$minmax)){
            # Read the file
            ff <- flowCore::read.FCS(file)

            # Transform if necessary
            if(!is.null(transformList)){
                ff <- flowCore::transform(ff, transformList)
            }

            # Overwrite the values with the normalized values
            if (verbose) message("Normalizing ",label)
            for (channel in channels) {
                flowCore::exprs(ff)[, channel] <-
                    ((flowCore::exprs(ff[, channel]) - model$minmax[[label]][1, channel])/
                         (model$minmax[[label]][2, channel] - model$minmax[[label]][1, channel])) *
                    (model$minmax[["goal"]][2, channel] - model$minmax[["goal"]][1, channel]) +
                    model$minmax[["goal"]][2, channel]
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
    }
}

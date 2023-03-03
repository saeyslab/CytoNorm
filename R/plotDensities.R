
#' Plot densities per batch
#'
#' @param input    Named list containing paths to fcs files, for both batch and
#'                 batch_norm. For a batch, the input can be either one or more
#'                 paths to fcs files or a flowFrame.
#' @param channels Channels to plot
#' @param colors Vector with a color for each batch
#' @param transformList In case the data from the fcs files still needs to be
#'                      transformed. Default NULL, in which case no transformation
#'                      happens.
#' @param suffix Suffixes used to distinguish the original from the normalized files.
#'               Default is c("original" = "", "normalized" = "_norm")
#'
#' @return  List with 2 plots per channel
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
#'                                               xdim = 15,
#'                                               ydim = 15,
#'                                               nClus = 10,
#'                                               scale = FALSE),
#'                         normParams = list(nQ = 101),
#'                         seed = 1,
#'                         verbose = TRUE)
#'
#' CytoNorm.normalize(model = model,
#'                    files = validation_data$Path,
#'                    labels = validation_data$Batch,
#'                    transformList = transformList,
#'                    transformList.reverse = transformList.reverse,
#'                    outputDir = "Normalized",
#'                    verbose = TRUE)
#'
#' original <- list("PTLG021" = validation_data$Path[1],
#'                  "PTLG028" = validation_data$Path[2],
#'                  "PTLG034" = validation_data$Path[3])
#' normalized <- list("PTLG021_norm" = paste0("Normalized/Norm_",validation_data$File[1]),
#'                    "PTLG028_norm" = paste0("Normalized/Norm_",validation_data$File[2]),
#'                    "PTLG034_norm" = paste0("Normalized/Norm_",validation_data$File[3]))
#'
#' channels_to_plot <- c("Er170Di", "La139Di")
#' plots <- plotDensities(input = c(original, normalized),
#'                        channels = channels_to_plot,
#'                        colors = c("blue", "red", "green"),
#'                        transformList = transformList)
#'
#' p <- ggpubr::ggarrange(ggpubr::ggarrange(plotlist = plots[1:(2*length(channels_to_plot))],
#'                                          ncol = 2,
#'                                          nrow = length(channels_to_plot)),
#'                        ggpubr::ggarrange(ggpubr::as_ggplot(plots[[length(plots)]])),
#'                        ncol = 1, nrow = 2, heights = c(10,0.5))
#' print(ggpubr::annotate_figure(p,
#'                       top = paste0("Before normalization",
#'                                    "                              ",
#'                                    "After normalization")))
#' @importFrom methods is
#' @importFrom flowCore read.FCS transform
#' @importFrom FlowSOM AggregateFlowFrames
#' @importFrom dplyr sym
#' @importFrom ggplot2 ggplot aes stat_density scale_color_manual xlab
#'                     theme_minimal theme .data
#' @importFrom ggpubr get_legend
#'
#' @export
plotDensities <- function(input, # list with 4 elements B1, B2, B1_norm, B2_norm
                          channels,
                          colors,
                          transformList = NULL,
                          suffix = c("original" = "",
                                     "normalized" = "_norm")){

    batch_names <- unique(gsub(suffix["original"], "",
                               gsub(suffix["normalized"], "", names(input))))

    data <- list("original" = list(),
                 "normalized" = list())

    for(batch in batch_names){
        for(type in c("original", "normalized")){
            i <- paste0(batch, suffix[type])
            if(is.character(input[[i]])){
                if(length(input[[i]] > 1)){
                    set.seed(2023)
                    data[[type]][[batch]] <- FlowSOM::AggregateFlowFrames(fileNames = input[[i]],
                                                              cTotal = length(input[[i]])*10000)
                } else {
                    data[[type]][[batch]] <- flowCore::read.FCS(input[[i]])
                }
            } else if(methods::is(input[[i]], "flowFrame")) {
                data[[type]][[batch]] <- input[[i]]
            }
            if(!is.null(transformList)){
                data[[type]][[batch]] <- flowCore::transform(data[[type]][[batch]],
                                                             transformList)
            }
        }
    }

    # Make dfs
    dfs <- list()
    for(type in c("original", "normalized")){
        dfs[[type]] <- data.frame(do.call(rbind,
                                          lapply(data[[type]],
                                                 function(x) flowCore::exprs(x))),
                                  check.names = FALSE)
        dfs[[type]]$Batch <- unlist(lapply(batch_names,
                                           function(i) rep(x = i,
                                                           times = nrow(data[[type]][[i]]))))

        dfs[[type]] <- dfs[[type]][,c(channels, "File", "Batch")]
    }

    # Plot
    plotlist <- list()
    for (channel in channels){
        for(type in c("original", "normalized")){
            df <- dfs[[type]]
            p <- ggplot2::ggplot(df, ggplot2::aes(x=!!dplyr::sym(channel), color = .data$Batch)) +
                ggplot2::stat_density(ggplot2::aes(group = paste(.data$Batch, .data$File)),
                                      geom = "line", position = "identity",
                                      alpha = 0.2)+
                ggplot2::stat_density(geom = "line", position = "identity",
                                      linewidth = 0.7) +
                ggplot2::scale_color_manual(values = colors)+
                ggplot2::xlab(channel) +
                ggplot2::theme_minimal()
            leg <- ggpubr::get_legend(p)
            p <- p + ggplot2::theme(legend.position = "none")
            plotlist[[paste(channel,type)]] <- p
        }
    }

    plotlist[["legend"]] <- leg
    return(plotlist)
}

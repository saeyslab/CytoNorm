
#' Plot densities per batch
#'
#' @param input    Named list containing paths to fcs files, for both batch and
#'                 batch_norm. For a batch, the input can be either one or more
#'                 paths to fcs files or a flowFrame.
#' @param channels Channels to plot
#' @param colors Optional vector with a color for each batch
#' @param transformList In case the data from the fcs files still needs to be
#'                      transformed. Default NULL, in which case no transformation
#'                      happens.
#' @param suffix Suffixes used to distinguish the original from the normalized files.
#'               Default is c("original" = "", "normalized" = "_norm")
#' @param model CytoNorm model to split out the normalized density plots per 
#'              FlowSOM meta-cluster. Optional. 
#' @param show_goal If a 'model' is provided, it is possible to show the model's 
#'                  goal distribution.
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
#'                        model = model,
#'                        channels = channels_to_plot,
#'                        colors = c("blue", "red", "green"),
#'                        transformList = transformList)
#'
#' p <- ggpubr::ggarrange(ggpubr::ggarrange(plotlist = plots[1:(length(plots)-1)],
#'                                          ncol = 2,
#'                                          nrow = 2*length(channels_to_plot)),
#'                        ggpubr::ggarrange(ggpubr::as_ggplot(plots[[length(plots)]])),
#'                        ncol = 1, nrow = 2, heights = c(10,1))
#'
#' @importFrom methods is
#' @importFrom flowCore read.FCS transform flowFrame
#' @importFrom FlowSOM AggregateFlowFrames NewData GetMetaclusters
#' @importFrom dplyr sym
#' @importFrom ggplot2 ggplot aes stat_density scale_color_manual xlab
#'                     theme_minimal theme .data xlim
#' @importFrom ggpubr get_legend
#'
#' @export
plotDensities <- function(input, # list with 4 elements B1, B2, B1_norm, B2_norm
                          channels,
                          colors = NULL,
                          model = NULL,
                          transformList = NULL,
                          show_goal = FALSE,
                          suffix = c("original" = "",
                                     "normalized" = "_norm")){

    batch_names <- unique(gsub(suffix["original"], "",
                               gsub(suffix["normalized"], "", names(input))))

    data <- list("original" = list(),
                 "normalized" = list())

    for(batch in batch_names){
        for(type in c("original", "normalized")){
            i <- paste0(batch, suffix[type])
            if(i %in% names(input)){
                if(is.character(input[[i]])){
                    if(length(input[[i]] > 1)){
                        set.seed(2023)
                        data[[type]][[batch]] <- FlowSOM::AggregateFlowFrames(fileNames = input[[i]],
                                                                              cTotal = length(input[[i]])*10000)
                    } else {
                        data[[type]][[batch]] <- flowCore::read.FCS(input[[i]], truncate_max_range = FALSE)
                    }
                } else if(methods::is(input[[i]], "flowFrame")) {
                    data[[type]][[batch]] <- input[[i]]
                }
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
        print(type)
        dfs[[type]] <- data.frame(do.call(rbind,
                                          lapply(data[[type]],
                                                 function(x) flowCore::exprs(x))),
                                  check.names = FALSE)
        dfs[[type]]$Batch <- unlist(lapply(batch_names,
                                           function(i) rep(x = i,
                                                           times = nrow(data[[type]][[i]]))))

        if(!is.null(model) & type == "original"){ #!is.null(dfs[[type]])){
            mapped <- FlowSOM::NewData(model$fsom,
                                       as.matrix(dfs[[type]][model$fsom$map$colsUsed]))
            dfs[[type]][,"Cluster"] <- FlowSOM::GetMetaclusters(mapped)
            dfs[[type]] <- dfs[[type]][,c(channels, "File", "Batch", "Cluster")]
        } else {
            dfs[[type]] <- dfs[[type]][,c(channels, "File", "Batch")]
        }

    }

    if(!is.null(model)){
        dfs[["normalized"]][,"Cluster"] <- dfs[["original"]][,"Cluster"]
    }

    if(show_goal){
        if(is.null(model)){
            stop("To show the goal, a model should be provided")
        }
        quantiles <- getCytoNormQuantiles(model)

    }

    # Plot
    plotlist <- list()
    for (channel in channels){
        x_range <- c(min(sapply(dfs, function(df) min(df[[channel]], na.rm = TRUE))),
                     max(sapply(dfs, function(df) max(df[[channel]], na.rm = TRUE))))
        for(type in c("original", "normalized")){
            df <- dfs[[type]]
            p <- ggplot2::ggplot(df, ggplot2::aes(x=!!dplyr::sym(channel),
                                                  color = .data$Batch)) +
                ggplot2::stat_density(ggplot2::aes(group = paste(.data$Batch,
                                                                 .data$File)),
                                      geom = "line", position = "identity",
                                      alpha = 0.2)+
                ggplot2::stat_density(geom = "line", position = "identity",
                                      linewidth = 0.7) +
                ggplot2::xlab(paste0(FlowSOM::GetMarkers(data[["original"]][[1]], channel),
                                     " <", channel, ">")) +
                ggplot2::ylab(type) +
                ggplot2::theme_minimal() +
                ggplot2::xlim(x_range)
            
            if (!is.null(colors)){
              p <- p +
                ggplot2::scale_color_manual(values = colors)
            }

            leg <- ggpubr::get_legend(p)
            p <- p + ggplot2::theme(legend.position = "none")
            plotlist[[paste(channel,type)]] <- p
            if(!is.null(model)){
                plotlist[[paste(channel,type,"per cluster")]] <-
                    ggplot2::ggplot(df, ggplot2::aes(x=!!dplyr::sym(channel),
                                                     color = .data$Batch))

                if(show_goal){
                    quantiles_df <- do.call(rbind,
                                            lapply(seq_along(quantiles),
                                                   function(x){
                                                       df_tmp <- data.frame(Value = 1/2 * (quantiles[[x]][-1,channel] + quantiles[[x]][-nrow(quantiles[[x]]),channel]),
                                                                            Density = 1 / nrow(quantiles[[x]])/ diff(quantiles[[x]][,channel]),
                                                                            Cluster = x)
                                                       df_tmp$Density_smooth <- stats::splinefun(df_tmp$Value, df_tmp$Density, method = "monoH.FC")(df_tmp$Value)
                                                       df_tmp
                                                   }))
                    plotlist[[paste(channel,type,"per cluster")]] <-
                        plotlist[[paste(channel,type,"per cluster")]] +
                        ggplot2::geom_line(aes(x = Value, y = Density_smooth),
                                           color = "black",
                                           data = quantiles_df)
                }


                plotlist[[paste(channel,type,"per cluster")]] <-
                    plotlist[[paste(channel,type,"per cluster")]] +
                    ggplot2::stat_density(ggplot2::aes(group = paste(.data$Batch,
                                                                     .data$File)),
                                          geom = "line", position = "identity",
                                          alpha = 0.2)+
                    ggplot2::stat_density(geom = "line", position = "identity",
                                          linewidth = 0.7) +
                    ggplot2::xlab(paste0(FlowSOM::GetMarkers(data[["original"]][[1]], channel),
                                         " <", channel, ">")) +
                    ggplot2::ylab(type) +
                    ggplot2::theme_minimal() +
                    ggplot2::xlim(x_range) +
                    ggplot2::facet_grid(~ .data$Cluster) + 
                    ggplot2::theme(legend.position = "none")
                
                if (!is.null(colors)){
                  plotlist[[paste(channel,type,"per cluster")]] <- plotlist[[paste(channel,type,"per cluster")]] +
                    ggplot2::scale_color_manual(values = colors)
                }
                
                
            }
        }


    }

    plotlist[["legend"]] <- leg
    return(plotlist)
}

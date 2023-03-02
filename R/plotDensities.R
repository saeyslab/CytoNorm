
#' Plot densities per batch
#'
#' @param input    Named list containing paths to fcs files, for both batch and batch_norm
#' @param channels Channels to plot
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
#'                    verbose = TRUE)
#'
#' plots <- plotDensities(input = list("PTLG021" = validation_data$Path[1],
#'                            "PTLG028" = validation_data$Path[2],
#'                            "PTLG034" = validation_data$Path[3],
#'                            "PTLG021_norm" = paste0("Normalized/Norm_",validation_data$File[1]),
#'                            "PTLG028_norm" = paste0("Normalized/Norm_",validation_data$File[2]),
#'                            "PTLG034_norm" = paste0("Normalized/Norm_",validation_data$File[3])),
#'              channels = c("Er170Di", "La139Di"),
#'              colors = c("blue", "red", "green"),
#'              transformList = transformList)
#' @importFrom methods is
#' @importFrom flowCore sampleNames flowSet
#'
#' @export
plotDensities <- function(input, # list with 4 elements B1, B2, B1_norm, B2_norm
                          channels,
                          colors,
                          transformList = NULL){

    batch_names <- unique(gsub("_norm", "", names(input)))

    data <- list()
    data_norm <- list()

    for(i in batch_names){
        if(is.character(input[[i]])){
            if(length(input[[i]] > 1)){
                set.seed(2023)
                data[[i]] <- FlowSOM::AggregateFlowFrames(fileNames = input[[i]],
                                                          cTotal = length(input[[i]])*10000)
            } else {
                data[[i]] <- read.FCS(input[[i]])
            }
        } else if(is.flowFrame(input[[i]])) {
            data[[i]] <- input[[i]]
        }
        if(!is.null(transformList)) data[[i]] <- flowCore::transform(data[[i]], transformList)

        i <- paste0(i,"_norm")
        if(is.character(input[[i]])){
            if(length(input[[i]] > 1)){
                set.seed(2023)
                data_norm[[i]] <- FlowSOM::AggregateFlowFrames(fileNames = input[[i]],
                                                          cTotal = length(input[[i]])*10000)
            } else {
                data_norm[[i]] <- read.FCS(input[[i]])
            }
        } else if(is.flowFrame(input[[i]])) {
            data_norm[[i]] <- input[[i]]
        }
        if(!is.null(transformList)) data_norm[[i]] <- flowCore::transform(data_norm[[i]], transformList)

    }

    # Make dfs
    df <- data.frame(do.call(rbind, lapply(data, function(x)flowCore::exprs(x))),
                     check.names = FALSE)
    df$Batch <- unlist(lapply(batch_names,
                              function(i) rep(x = i, times = nrow(data[[i]]))))

    df <- df[,c(channels, "File", "Batch")]

    df_norm <- data.frame(do.call(rbind, lapply(data_norm, function(x)flowCore::exprs(x))),
                          check.names = FALSE)
    df_norm$Batch <- unlist(lapply(batch_names,
                                   function(i) rep(x = i, times = nrow(data_norm[[paste0(i, "_norm")]]))))
    df_norm <- df_norm[,c(channels, "File", "Batch")]

    # Plot
    plotlist <- list()
    for (channel in channels){
        p <- ggplot2::ggplot(df, ggplot2::aes(x=!!dplyr::sym(channel), color = Batch)) +
            ggplot2::stat_density(ggplot2::aes(group = paste(Batch, File)),
                                  geom = "line", position = "identity",
                                  alpha = 0.2)+
            ggplot2::stat_density(geom = "line", position = "identity",
                                  linewidth = 0.7) +
            ggplot2::scale_color_manual(values = colors)+
            ggplot2::xlab(channel) +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = "none")
        plotlist[[channel]] <- p

        p <- ggplot2::ggplot(df_norm, ggplot2::aes(x=!!dplyr::sym(channel), color = Batch)) +
            ggplot2::stat_density(ggplot2::aes(group = paste(Batch, File)),
                                               geom = "line", position = "identity", alpha = 0.2)+
            ggplot2::stat_density(geom = "line", position = "identity", size = 0.7) +
            #stat_density(geom = "line", position = "identity", linewidth = 0.8) +
            ggplot2::scale_color_manual(values = colors)+
            ggplot2::xlab(channel) +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.direction = "horizontal")
        leg <- ggpubr::get_legend(p)
        p <- p + ggplot2::theme(legend.position = "none")
        plotlist[[paste0(channel,"_norm")]] <- p
    }

    plotlist[["legend"]] <- leg
    return(plotlist)
}

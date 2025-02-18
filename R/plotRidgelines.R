#' Plot ridgeline plots
#'
#' @param input    Paths to fcs files.
#' @param batch Optional vector with a batch label for each file
#' @param channels Channels to plot
#' @param colors Optional vector with a color for each batch
#' @param transformList In case the data from the fcs files still needs to be
#'                      transformed. Default NULL, in which case no transformation
#'                      happens.
#' @param quantiles Optional vector specifying the quantiles that need to be
#'                  highlighted and connected.
#'
#' @return  Named list with a plot per channel
#'
#' @examples
#'
#' dir <- system.file("extdata", package = "CytoNorm")
#' files <- list.files(dir, pattern = "fcs$")
#' ff <- flowCore::read.FCS(file.path(dir, files[1]))
#' channels <- grep("Di$", flowCore::colnames(ff), value = TRUE)
#' transformList <- flowCore::transformList(channels,
#'                                          cytofTransform)
#' beforeNorm <- plotRidgelines(input = file.path(dir, files), 
#'                              batch = stringr::str_match(files, "PTLG[0-9]*")[,1],
#'                              channels = c("Er170Di", "La139Di"),
#'                              colors = c("PTLG021" = "#d8e2dc", 
#'                                         "PTLG028" = "#ffe5d9",
#'                                         "PTLG034" = "#ffcad4"),
#'                              transformList = transformList,
#'                              quantiles = c(0.01, 0.25, 0.5, 0.75, 0.99))
#'                              
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
#'                    outputDir = "Normalized",
#'                    verbose = TRUE)
#'
#' files <- list.files("Normalized", pattern = "fcs$")
#' afterNorm <- plotRidgelines(input = file.path("Normalized", files), 
#'                             batch = stringr::str_match(files, "PTLG[0-9]*")[,1],
#'                             channels = c("Er170Di", "La139Di"),
#'                             colors = c("PTLG021" = "#d8e2dc", 
#'                                        "PTLG028" = "#ffe5d9",
#'                                        "PTLG034" = "#ffcad4"),
#'                             transformList = transformList,
#'                             quantiles = c(0.01, 0.25, 0.5, 0.75, 0.99))
#'
#' p <- list()
#' for (i in 1:length(beforeNorm)){
#'   p[[length(p)+1]] <- ggpubr::ggarrange(beforeNorm[[i]], afterNorm[[i]], nrow =1)
#' }
#'
#' @importFrom methods is
#' @importFrom flowCore transform exprs
#' @importFrom FlowSOM AggregateFlowFrames
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr arrange
#' @importFrom stats quantile
#' @importFrom ggplot2 ggplot aes scale_fill_manual ggplot_build geom_path
#'                     geom_point xlab
#'                     theme_minimal theme_minimal

#' @export
#' 
plotRidgelines <- function(input, 
                           batch = NULL,
                           channels,
                           colors = NULL,
                           transformList = NULL,
                           quantiles = c(seq(0,1,0.1))){
  
  
  data <- FlowSOM::AggregateFlowFrames(fileNames = input,
                                       cTotal = length(input)*10000)
  if(!is.null(transformList)){
    data <- flowCore::transform(data, transformList)
  }
  
  df <- data.frame(flowCore::exprs(data)[,c(channels, "File")], check.names = FALSE)
  df$FileName <- sub("(?i).*/(.*).fcs", "\\1", input)[df$File]
  df$Batch <- batch[df$File]
  
  df_l <- data.frame(tidyr::pivot_longer(df, 
                                         cols = -c("File", "FileName", "Batch"), 
                                         names_to = "Channel", 
                                         values_to = "Intensity"),
                     check.names = FALSE)
  
  plotlist <- list()
  
  for (channel in channels){
    fileQuantiles <- matrix(data = NA,
                            nrow = length(input), ncol = length(quantiles)+1,
                            dimnames = list(unique(df$FileName), c(paste0("q", quantiles), "File")))
    for (filename in unique(df$FileName)){
      fileQuantiles[filename,1:length(quantiles)] <- quantile(df[df$FileName == filename, channel], quantiles)
    }
    fileQuantiles[,"File"] <- 1:nrow(fileQuantiles)
    fileQuantiles_l <- data.frame(tidyr::pivot_longer(data.frame(fileQuantiles), 
                                                      cols = colnames(fileQuantiles)[startsWith(colnames(fileQuantiles), "q")], 
                                                      names_to = "Quantile", values_to = "x"))
    fileQuantiles_l$y <- NA
    
    p <- ggplot2::ggplot(df_l[df_l$Channel == channel,], ggplot2::aes(x = Intensity, y = FileName, fill = Batch)) +
      ggridges::geom_density_ridges() +
      ggplot2::scale_fill_manual(values = colors) +
      ggridges::theme_ridges() 
    p_build <- ggplot2::ggplot_build(p)$data[[1]] # Get plot coordinates
    for(row in 1:nrow(fileQuantiles_l)){
      file <- fileQuantiles_l$File[row]
      p_build_sub <- p_build[p_build$y == file,]
      i <- which(abs(p_build_sub$x - fileQuantiles_l[row, "x"]) == 
                   min(abs(p_build_sub$x - fileQuantiles_l[row, "x"])))
      fileQuantiles_l[row,"y"] <- p_build_sub$ymax[i]
    }
    fileQuantiles_l[,"FileName"] <- sub("(?i).*/(.*).fcs", "\\1", input)[fileQuantiles_l[,"File"]]
    fileQuantiles_l$Quantile <- as.factor(fileQuantiles_l$Quantile)
    fileQuantiles_l <- dplyr::arrange(fileQuantiles_l, Quantile, y)
    
    plotlist[[channel]] <- ggplot2::ggplot(df_l[df_l$Channel == channel,], ggplot2::aes(x = Intensity, y = FileName)) +
      ggridges::geom_density_ridges(ggplot2::aes(fill = Batch)) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::geom_path(data = fileQuantiles_l, 
                         position = "identity",
                         ggplot2::aes(x=x, y=y, group = Quantile, color = Quantile)) +
      ggplot2::geom_point(data = fileQuantiles_l,
                          ggplot2::aes(x = x, y = y, color = Quantile)) +
      ggplot2::theme_minimal()
    
  }
  return(plotlist)
}

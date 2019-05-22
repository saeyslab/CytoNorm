
#' emdEvaluation
#'
#' Evaluate how much files differ by computing the maximum Earth Movers Distance
#' for all markers and cellTypes.
#'
#' @param files     Full paths of to the fcs files of the control samples.
#' @param transformList Transformation list to pass to the flowCore
#'                  \code{transform} function
#' @param markers   Markers to evaluate (corresponding with the column
#'                  names of the flow frame)
#' @param manual    A list which contains for every file a factor array. These
#'                  arrays contain a cell label for every cell in the files. All
#'                  arrays should have the same levels. Default = NULL, all
#'                  cells are evaluated together.
#' @param binSize   Binsize to approximate distribution. Default = 0.1.
#' @param return_all If TRUE, distributions and pairwise distances are returned
#'                   as well. Default = FALSE.
#'
#' @return A matrix in which the rows represent the cell types, the columns
#' reprents the markers and the values represent the maximal earth movers
#' distances for the distributions between all files
#'
#' @examples
#'    # Describe file names
#'    dir <- system.file("extdata",package="CytoNorm")
#'    fileNames <- c("Gates_PTLG021_Unstim_Control_1.fcs",
#'                    "Gates_PTLG026_Unstim_Control_1.fcs")
#'    labels <- c("PTLG021","PTLG026")
#'    ff <- flowCore::read.FCS(file.path(dir,fileNames[1]))
#'    markersToNormalize <- flowCore::colnames(ff)[c(10, 11, 14, 16:35, 37, 39:51)]
#'
#'    # Build transform list
#'    transformList <- flowCore::transformList(markersToNormalize,
#'                                         flowCore::arcsinhTransform(
#'                                            transformationId="cytofTransform",
#'                                            a=0,b=(1/5),c=0))
#'    emdEvaluation(file.path(dir,fileNames),
#'                  transformList,
#'                  markersToNormalize)
#'
#' @export
emdEvaluation <- function(files,
                          transformList,
                          markers,
                          manual = NULL,
                          binSize = 0.1,
                          return_all = FALSE){

    if(is.null(manual)){
        cellTypes <- c("AllCells")
    } else {
        cellTypes <- levels(manual[[1]])
    }

    distr <- list()
    for(file in files){
        print(file)

        distr[[file]] <- list()

        ff <- flowCore::read.FCS(file)
        ff <- flowCore::transform(ff,transformList)
        for(cellType in cellTypes){
            if(is.null(manual)){
                selection <- seq_len(flowCore::nrow(ff))
            } else {
                selection <- manual[[gsub("^Norm_","",gsub(".*/","",file))]]==cellType
            }
            distr[[file]][[cellType]] <-
                apply(flowCore::exprs(ff)[selection,
                                          markers],
                      2,
                      function(x){
                          graphics::hist(x,
                                         breaks = seq(-100,100,by=binSize),
                                         plot = FALSE)$counts
                      })
            any(distr[[file]][[cellType]] != 0)
        }
    }

    distances <- list()
    for(cellType in cellTypes){
        distances[[cellType]] <- list()
        for(marker in markers){
            distances[[cellType]][[marker]] <- matrix(NA,
                                                      nrow=length(files),
                                                      ncol=length(files),
                                                      dimnames=list(files,
                                                                    files))

            for(i in seq_along(files)[-length(files)]){
                file1 <- files[i]
                for(j in seq(i+1,length(files))){
                    file2 <- files[j]
                    distances[[cellType]][[marker]][file1,file2] <-
                        emdist::emd2d(
                            matrix(distr[[file1]][[cellType]][,marker]),
                            matrix(distr[[file2]][[cellType]][,marker]))
                }
            }
        }
    }

    comparison <- matrix(NA,
                         nrow=length(cellTypes),
                         ncol=length(markers),
                         dimnames = list(cellTypes, markers))

    for(cellType in cellTypes){
        for(marker in markers){
            comparison[cellType,marker] <- max(distances[[cellType]][[marker]],
                                               na.rm = TRUE)
        }
    }

    if(return_all){
        return(named.list(distr, distances, comparison))
    } else {
        return(comparison)
    }
}

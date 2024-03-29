% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluation.R
\name{madEvaluation}
\alias{madEvaluation}
\title{madEvaluation}
\usage{
madEvaluation(
  files,
  channels,
  transformList = NULL,
  prefix = "^Norm_",
  manual = NULL,
  return_all = FALSE
)
}
\arguments{
\item{files}{Full paths of to the fcs files of the control samples.}

\item{channels}{Channels to evaluate (corresponding with the column
names of the flow frame)}

\item{transformList}{Transformation list to pass to the flowCore
\code{transform} function. Default NULL}

\item{prefix}{Prefix present in the files, which will be removed to match
the manual list.}

\item{manual}{A list which contains for every file a factor array. These
arrays contain a cell label for every cell in the files. All
arrays should have the same levels. Default = NULL, all
cells are evaluated together.}

\item{return_all}{If TRUE, individual MAD values are returned
as well. Default = FALSE.}
}
\value{
A matrix in which the rows represent the cell types, the columns
reprents the markers and the values represent the median MAD values for
the distributions of all files
}
\description{
Evaluate whether you lose biological information by checking whether the
MAD stays similar before and after normalization.
}
\examples{
   # Describe file names
   dir <- system.file("extdata",package="CytoNorm")
   fileNames <- c("Gates_PTLG021_Unstim_Control_1.fcs",
                   "Gates_PTLG028_Unstim_Control_1.fcs")
   labels <- c("PTLG021","PTLG028")
   ff <- flowCore::read.FCS(file.path(dir,fileNames[1]))
   channelsToNormalize <- flowCore::colnames(ff)[c(10, 11, 14, 16:35, 37, 39:51)]

   # Build transform list
   transformList <- flowCore::transformList(channelsToNormalize,
                                            cytofTransform)
   res <- madEvaluation(file.path(dir,fileNames),
                        transformList = transformList,
                        channelsToNormalize)

}

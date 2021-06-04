#' cytofTransform
#'
#' Arcsinh transformation with cofactor 5
#'
#' @param x Value to transform
#' @export
cytofTransform <- flowCore::arcsinhTransform(transformationId="cytofTransform",
                                             a=0, b=(1/5), c=0)

#' cytofTransform.reverse
#'
#' @param x Value to tranform
#' Function to reverse arcsinh transformation with cofactor 5
#'
#' @export
cytofTransform.reverse <- function(x){
    return(sinh(x)/(1/5))
}

#' named.list
#'
#' Creates a list, in which the names of the items correspond to the names
#' of the variables.
#'
#' @param ... variables to add in the list
#' @export
named.list <- function(...) {
    l <- list(...)
    names(l) <- sapply(substitute(list(...)), deparse)[-1]
    l
}

#' textPlot
#'
#' Creates an empty plot with a label in it
#'
#' @param text Text to show in the plot
#' @export
textPlot <- function(text){
    graphics::plot(0,
                   xlim = c(-1,1), ylim = c(-1,1),
                   type = "n",
                   xaxt = "n", yaxt = "n",
                   bty = "n", xlab = "", ylab = "")
    graphics::text(0, 0, text, adj = 0.5)
}

#' identityFunction
#'
#' Returns the input value
#' @param x The value to return
#' @export
identityFunction <- function(x){
    x
}

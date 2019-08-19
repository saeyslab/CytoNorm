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

#' applyStaticGating
#'
#' Parse a flowjo workspace and apply the first gating on a different file
#' @param wsp_path Path of the flowjo workspace
#' @param fcs_path Path of the fcs file
#'
#' @return A matrix of TRUE/FALSE values, with one row for every cell and one
#'         column for every gate.
#' @export
applyStaticGating <- function(wsp_path,
                              fcs_path){
    wsp <- CytoML::open_flowjo_xml(wsp_path)
    gs <- CytoML::flowjo_to_gatingset(wsp,
                                      sampNloc = "sampleNode",
                                      name = "All Samples",
                                      execute = FALSE)

    gs_new <-  flowWorkspace::gh_apply_to_new_fcs(gs[[1]],
                                                  fcs_path)
    gate_names <- flowWorkspace::gs_get_pop_paths(gs[[1]], path = "auto")[-1]
    gatingMatrix <-
        flowWorkspace::gh_pop_get_indices_mat(gs_new[[1]],
                                              paste(gate_names, collapse = "|"))


    CytoML::flowjo_ws_close(wsp)

    return(gatingMatrix)
}


#' jaccard
#'
#' Computes the jaccard index given two logical vectors of equal length:
#' intersection / union
#'
#' @param a logical vector
#' @param b logical vector
#'
#' @examples
#' a <- c(TRUE, TRUE, FALSE, FALSE, TRUE)
#' b <- c(TRUE, FALSE, TRUE, FALSE, TRUE)
#' jaccard(a, b)
#'
#' @export
jaccard <- function(a, b){
    sum(a & b) / sum(a | b)
}

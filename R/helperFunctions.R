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

    CytoML::flowjo_ws_close(wsp)

    gs_new <-  flowWorkspace::gh_apply_to_new_fcs(gs[[1]],
                                                  fcs_path)
    gate_names <- flowWorkspace::gs_get_pop_paths(gs[[1]], path = "auto")[-1]
    gatingMatrix <-
        flowWorkspace::gh_pop_get_indices_mat(gs_new[[1]],
                                              paste(gate_names, collapse = "|"))

    return(gatingMatrix)

    ############################################################################

        # wsp <- CytoML::open_flowjo_xml(wsp_path)
        # o <- capture.output(
        #     gs <- suppressMessages(
        #         CytoML::flowjo_to_gatingset(wsp,
        #                                     sampNloc = "sampleNode",
        #                                     name = "All Samples")))
        # CytoML::flowjo_ws_close(wsp)
        #
        # gate_names <- flowWorkspace::gs_get_pop_paths(gs[[1]], path = "full")[-1]
        #
        # ff <- flowCore::read.FCS(fcs_path)
        # ff <- flowCore::transform(ff,
        #                           transformList(names(flowWorkspace::gh_get_transformations(gs[[1]])),
        #                                         flowWorkspace::gh_get_transformations(gs[[1]])))
        #
        #
        # # Why is this line needed???
        # ff@exprs[,grep("Di", flowCore::colnames(ff))] <- ff@exprs[,grep("Di", flowCore::colnames(ff))] / 32.61719
        #
        # gs_new <- flowWorkspace::GatingSet(flowSet(ff))
        #
        # for(gate_name in gate_names){
        #     gate <- flowWorkspace::gh_pop_get_gate(gs[[1]], gate_name)
        #     parent <- gsub("/[^/]*$", "", gate_name)
        #     if(parent == "") parent <- "root"
        #     flowWorkspace::gs_pop_add(gs = gs_new,
        #                               gate = gate,
        #                               parent = parent)
        # }
        # flowWorkspace::recompute(gs_new)
        #
        # gatingMatrix <- matrix(NA,
        #                        nrow = nrow(ff),
        #                        ncol = length(gate_names),
        #                        dimnames = list(NULL,
        #                                        gate_names))
        # for(gate_name in gate_names){
        #     gatingMatrix[, gate_name] <- flowWorkspace::gh_pop_get_indices(
        #         gs_new[[1]],
        #         gate_name)
        #
        #     gate <- flowWorkspace::gh_pop_get_gate(gs[[1]], gate_name)
        #     # if("boundaries" %in% slotNames(gate)){
        #     #     plot(ff@exprs[,colnames(gate@boundaries)], pch = ".",
        #     #          main = sum(gatingMatrix[,gate_name]))
        #     #     lines(gate@boundaries, col = "red")
        #     #     print(gate@boundaries)
        #     # } else if("min" %in% slotNames(gate)) {
        #     #     plot(ff@exprs[1:10000,names(gate@min)], pch = ".",
        #     #          main = sum(gatingMatrix[,gate_name]))
        #     #     lines(c(gate@min[1], gate@min[1], gate@max[1], gate@max[1], gate@min[1]),
        #     #           c(gate@min[2], gate@max[2], gate@max[2], gate@min[2], gate@min[2]),
        #     #           col = "red")
        #     # }
        # }
        # colnames(gatingMatrix) <- flowWorkspace::gs_get_pop_paths(gs[[1]],
        #                                                           path = "auto")[-1]
        # return(gatingMatrix)
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

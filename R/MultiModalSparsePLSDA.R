#' Multi-modal supervised integration of MultiAssayExperiment assays
#'
#' @inheritParams MultiModalSparsePLS
#' @param formula A formula of form LHS ~ RHS where LHS is a name from colData(data) 
#' and RHS is a list of assays (e.g. assay_1 + assay_2 + ...). 
#' @param ... Other arguments passed to mixOmics::block.splsda
#'
#' @return
#' @export
#'
#' @example ./examples/MultiModalSparsePLSDA-examples.R
MultiModalSparsePLSDA <- function(data,
                                  formula,
                                  keep_features = NULL,
                                  ncomp = 2,
                                  design='full',
                                  scale = TRUE,
                                  ...) {
    mc <- mget(names(formals()), sys.frame(sys.nframe()))
    if(any(!MultiAssayExperiment::complete.cases(data))) {
        stop('Some samples do not have complete observations.',
             ' Use MultiAssayExperiment::intersectColumns to match the observations', call. = FALSE)
    }
    mc <- .get_xy(mc = mc, DA = TRUE, block = TRUE)
    mc$data <- mc$formula <- NULL 
    
    if (!is.null(keep_features) & (names(keep_features) %!=% names(mc$X))) {
        stop("'keep_features' names do not match the assays in the RHS of formula.")
    }
    
    mc$X <- .get_primary_names(mae = data, X = mc$X)
    
    result <-
        block.splsda(
            X = mc$X,
            Y = mc$Y,
            ncomp = ncomp,
            keepX = mc$keep_features,
            design = design,
            scale = mc$scale,
            ...
        )
    
    result$call <- match.call()
    class(result) <- c(class(result), 'multimodal')
    return(result)
    
}
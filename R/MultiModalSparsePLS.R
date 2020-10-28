#' Multi-modal integration of MultiAssayExperiment assays
#'
#' @inheritParams block.spls
#' @param data A \code{MultiAssayExperiment} object with matching cells/samples
#' across all experiments.
#' @param formula A formula of form LHS ~ RHS where LHS is an assay and RHS is 
#' a list of assays (e.g. assay_1 + assay_2 + ...)
#' @param keep_features Named list of number of features to keep on each component
#' and assay. See examples.
#' @param design Matrix or character. Connects assays during integration.
#' Each value is a continuous number b/w 0 (not connected) and 1 (fully connected).
#' 'null' will connect the LHS to all RHS only and 'full' will fully connect all
#' assays.
#' @param ... Other arguments passed to mixOmics::block.spls
#'
#' @return A 'multimodal' object
#' @export
#' @example ./examples/MultiModalSparsePLS-examples.R
#' 
MultiModalSparsePLS <- function(data,
                                formula,
                                keep_features = NULL,
                                ncomp = 2,
                                design='full',
                                scale = TRUE,
                                ...) {
    require(MultiAssayExperiment)
    mc <- mget(names(formals()), sys.frame(sys.nframe()))
    .check_data_rows_and_cols(data = data)
    
    mc <- .get_xy(mc = mc, DA = FALSE, block = TRUE)
    mc$data <- mc$formula <- NULL 
    
    mc$X <- c(list('Y'=mc$Y), mc$X)
    names(mc$X)[1] <- as.character(formula)[2]
    mc$Y <- NULL
    
    if (!is.null(keep_features) & (names(keep_features) %!=% names(mc$X))) {
        stop("'keep_features' names do not match the assays in formula.")
        }
    
    mc$X <- .get_primary_names(mae = data, X = mc$X)
    
    result <-
        block.spls(
            X = mc$X,
            indY = 1,
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
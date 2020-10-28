#' Tuning function for multi-modal supervised integration of MultiAssayExperiment assays
#'
#' @inheritParams MultiModalSparsePLSDA
#' @param test_keep_features Named list of integer vectors specifying the search grid
#' @param ... Other arguments passed ro tune.block.splsda
#' @return A list of class'tune'
#' @example ./examples/TuneMultiModalSparsePLSDA-examples.R
#' @export
TuneMultiModalSparsePLSDA <- function(data,
                                      formula,
                                      test_keep_features = NULL,
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
    
    if (names(test_keep_features) %!=%  names(mc$X)) {
        stop("'test_keep_features' names do not match the assays in formula.")
    }
    
    mc$X <- .get_primary_names(mae = data, X = mc$X)
    
    result <-
        tune.block.splsda(
            X = mc$X,
            Y = mc$Y,
            ncomp = ncomp,
            test.keepX = test_keep_features,
            design = design,
            scale = mc$scale,
            ...
        )
    
    result$call <- match.call()
    return(result)
}

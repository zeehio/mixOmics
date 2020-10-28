#' create design matrix for multi-modal pls
#'
#' @param X Named list of datasets. When
#' type = 'null', for non-DA analyses the first one is 
#' taken to be the response matrix which is fully connected to others.
#' @param type Character, one of c('full', 'null'). If 'full' all modalities will
#' be connected, otherwise only the first modality is connected. Default is 'null'.
#'
#' @return A matrix whose dimension names match the names of \code{keep_features}
#' @noRd
#' @examples
#' .create_design(list(rna = c(20, 30), met = c(40, 20), acc = c(20, 40)), type = 'full')
.create_design <- function(X, type = 'null', DA = FALSE) {
    modals <- names(X)
    
    design <-  matrix(1, nrow = length(modals), ncol = length(modals) ,dimnames = rep(list(modals), 2))
    diag(design) <- 0
    ## null - only weights on Y
    if (type != 'full') {
        design <- matrix(0, nrow = length(modals), ncol = length(modals) ,dimnames = rep(list(modals), 2))
        if (isFALSE(DA)) {
            design[1,-1] <- design[-1,1] <- 1
        }
    }
    return(design)
}
## ------------------------------------------------------------------------ ##
## Match the names of observations/cells/samples for all modalities
.get_primary_names <- function(mae, X) {
    sm <- MultiAssayExperiment::sampleMap(mae)
    
    res <- mapply(assay_value=X, assay_type=names(X), function(assay_value, assay_type) {
        sm <- sm[(sm$assay == assay_type), ]
        sm <- data.frame(names = sm$primary, row.names = sm$colname)
        rownames(assay_value) <- sm[rownames(assay_value),,drop=FALSE]$names
        assay_value
    }, SIMPLIFY = FALSE)
    res
}
## ------------------------------------------------------------------------ ##
.check_data_rows_and_cols <- function(data) {
    if(any(!MultiAssayExperiment::complete.cases(data))) {
        stop('Some samples do not have complete observations.',
             ' Use MultiAssayExperiment::intersectColumns to match the observations', call. = FALSE)
    }
    
    rn<- lapply(MultiAssayExperiment::assays(data), rownames)
    factor_rn <- factor(names(rn), ordered = TRUE)
    
    for (x in  factor_rn) {
        for (y in factor_rn) {
            if ( x < y) {
                commons <- intersect(rn[[x]], rn[[y]])
                
                if (length(commons) > 1) {
                    warning("assays ", x, " and ", y, " have common feature names. ",
                            "This creates issues in some visualisations. ", call. = FALSE)
                    # warning(sprintf("some assays have common feature names. ",
                    #                 "This might create issues in some visualisations. ",
                    #                 "Consider creating unique feature names for all assays."), call. = FALSE)
                } 
            }
        }
    }
}

## ------------------------------------------------------------------------ ##

#' Make unique feature names for assays
#'
#' @param data A \code{MultiAssayExperiment} object
#'
#' @return A \code{MultiAssayExperiment} object with unique feature names
#'  in all assays
#'
#' @examples
#' \dontrun{
#' library(MultiAssayExperiment)
#' data("miniACC")
#' miniACC <- MatchedAssayExperiment(miniACC)
#' make_unique_feature_names(miniACC)
#' }
#' @export
make_unique_feature_names <- function(data) {

    
    has_duplicate_feats <- .check_duplicate_feature_names(data)
    if (any(has_duplicate_feats)) {
        
        experiments <- experiments(data)
        col_data <- colData(data)
        sample_map <- sampleMap(data)
        meta_data <- metadata(data)
        
        for (j in seq_along(has_duplicate_feats)) {
            if (isTRUE(has_duplicate_feats[j])) {
                message('Making unique feature names by adding assay name to feature name for: ', names(experiments)[j])
                rownames(experiments[[j]]) <- paste0(rownames(experiments[[j]]), '_', names(experiments)[j])
            }
        }
        
        res <- MultiAssayExperiment::MultiAssayExperiment(experiments = experiments, 
                                                          colData = col_data, sampleMap = sample_map, 
                                                          metadata = meta_data)
        return(res)
    } else {
        message("No duplicate feature names")
        return(data)
    }
    
}



## ------------------------------------------------------------------------ ##
.check_duplicate_feature_names <- function(data) {
    experiments <- experiments(data)
    
    has_duplicate_feats <- vector(length = length(experiments))
    
    for (i in  seq_along(experiments)) {
        for (j in seq_along(experiments)) {
            if ( j > i) {
                rn<- lapply(experiments, rownames)
                commons <- intersect(rn[[i]], rn[[j]])
                
                if (length(commons) > 1) {
                    has_duplicate_feats[j] <- TRUE
                } 
            }
        }
    }
    
    has_duplicate_feats
}
#' @noRd
'standardise_to_-1_and_1' <- function(vec) 
{
    2*(vec - min(vec, na.rm = TRUE))/ (max(vec, na.rm = TRUE) - min(vec, na.rm = TRUE)) - 1
}

#' @noRd
'standardise_columns_to_-1_and_1' <- function(x) 
{
    apply(x, 2, 'standardise_to_-1_and_1')
}

#' Get block alignment in latent space
#'
#'
#' @param object A pls-derived object
#'
#' @param X.block X block name
#' @param Y.block Reference block name
#' @param weighted Logical, whether the distances should be weighted by the 
#' explained variance of Y components.
#'
#' @noRd
#' @examples 
#' \dontrun{
#' }
.get_pls_alignment <-
    function(X.block = 'X',
             object,
             Y.block = 'Y',
             weighted = TRUE
    )
{
    weighted_manhattan <- function(x, w=NULL) {
        if (is.null(w))
            w <- rep(1, ncol(x))
        
        ## w to have a mean of 1
        w <- w/mean(w)
        dist(x*w, method = 'manhattan')
    }
    
    variates <-  object$variates
    variates <-  lapply(variates, 'standardise_columns_to_-1_and_1')
    x <- variates[[X.block]]
    y <- variates[[Y.block]]
    w <- object$explained_variance[[Y.block]]
    ## median weighted manhattan distance between x and y variates
    dev <- weighted_manhattan(x = abs(x-y), w = if (weighted) w else NULL)
    ## median weighted manhattan distance between y samples
    y.dist <- weighted_manhattan(x = y, w = if (weighted) w else NULL)
    median.y.dist <- median(as.vector(y.dist))
    
    ## deviation to median distance
    as.vector(dev/median.y.dist)
    }

get_pls_alignment <-
    function(object,
             Y.block = 'Y',
             weighted = TRUE
    )
    {
        blocks <- names(object$variates)[-which(names(object$variates) == Y.block)]
        lapply(.name_list(blocks), function(x){
            .get_pls_alignment(
                X.block = x,
                object = object,
                Y.block = Y.block,
                weighted = weighted
            )
        })
    }

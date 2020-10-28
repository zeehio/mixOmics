## ----------- .formula_checker -----------
#Check formula's class and ensure non-NULL X and Y are not provided with it
#
#'Gets a call which includes a $formula entry and expects it to be of form
#\code{Y~X} or \code{Y~X1+X2+...} (for block), it also checks that $X and $Y are NULL.
#
#@param mc A call list
#
#@return Exception handler
#
#' @noRd
.formula_checker <- function(mc, block=FALSE){
  fl <- vapply(as.list(mc$formula)[-1L], length, 1) ## formula list lengths
  
  if (isTRUE(block)) {
    ## check formula is Y~X1+X2+...
    if ( length(fl) != 2 || fl[1] != 1 || fl[2] < 2)
      .inv_bformula()
  } else {
    ## check formula is Y~X
    if ( length(fl) != 2 || fl[1] != 1 || fl[2] != 1)
      .inv_sformula()
  }
  ## X and Y must be NULL
  if (!all(vapply(mc[c("X", "Y")], is.null, TRUE)))
    .inv_signature()
}


## ----------- .call_return ----------- 
## TODO keep only one returner function
## function to add the call to the result if asked
## The aim is to reduce the methods code
## Handling it in here to avoid unnecessary copy of the fully evaluated call
## unless user asks for it. No full copy is made unless ret.call=TRUE.
.call_return <-
  function(result = list(), ## result from internal
           mcr, ## always match.call()
           fun.name = 'pca', ## name of the function
           pframe = 3L) {  ## how many frames above to evaluate in?
    
    if ( isTRUE(mcr$ret.call) ) {
      {
        mcr[[1]] = as.name(fun.name)
        mcr[-1L] <- lapply(mcr[-1L], function(x) eval.parent(x, n = pframe))
      }
      ## makes expect_identical() easy to have call as first arg and just
      ## drop it using pca.res[-1]
      
      result <- structure(c(list(call = mcr), result), class = class(result))
    }
    return(result)
  }

## ----------- .matched_samples ----------- 
#' Keep onlymatched samples in MAE object
#'
#' Creating this so I would not have to import every time
#' 
#' 
# @param mae 
#' @importFrom MultiAssayExperiment MatchedAssayExperiment
#' @noRd
.matched_samples <- function(mae) {
  message("Keeping matching samples in MultiAssayExperiment object ...")
  mae <- MatchedAssayExperiment(mae)
}

## ------------------------------------------------ ##
#' @title Get matrix from assay name
#'
#' @param assay.name Name of assay
#' @param mae_data MAE data
#' @importFrom MultiAssayExperiment assay
#' 
#' @return tansposed assay data
#' @noRd
.getMatrix <- function(assay.name, mae_data) {
  as.matrix(t(assay(mae_data, assay.name)))
}

## ------------------------------------------------ ##
#' @title Given data in mc, get X for PCA family
#'
#' @param mc A match.call
#' @return Adjusted mc with \code{(N x P)}  \code{X} to pass to internal
#' @importFrom MultiAssayExperiment colData assays
#' 
#' @noRd
.get_x <- function(mc) {
  
  ## ensure X and Y are character
  if (!(is(mc$X, "character") & length(mc$X) == 1 )) {
    .stop("'X' must be one valid assay name from 'data'", 
          .subclass = "inv_XY")
  }
  
  ## ensure X is a valid assay
  invalid_X <- !(mc$X %in% names(assays(mc$data)))
  if (any(invalid_X)) {
    .stop(.subclass = "inv_XY",
          message = sprintf("%s is not a valid assay name from 'data'", sQuote(mc$X)))
  }
  
  ## get X
  mc$X <- .getMatrix(assay.name = mc$X, mc$data)
  mc
}
## ------------------------------------------------ ##
#' @title Given a formula, get assay names for .get_xy
#'
#' @param mc A match.call
#'
#' @return Adjusted mc to pass to .get_xy
#' @noRd
.get_xy_names_from_formula <- function(mc, block=FALSE) {
  .formula_checker(mc, block = block) ## check formula validity
  ## esnure X and Y are NULL
  if ( !(is.null(mc$X) && is.null(mc$Y)) )
    .stop(message = "Where 'data' and 'formula' are provided 'X' and 'Y' should be NULL.", 
          .subclass = "inv_signature")
  
  assay.names <- rownames(attr(terms(mc$formula), "factors"))
  mc$Y <- assay.names[1]
  mc$X <- assay.names[-1]
  mc
}

## ------------------------------------------------ ##
#' @title Given a formula of matrices, get X and Y
#'
#' @param mc A match.call
#' @param DA Logical, TRUE if DA analysis
#' @param block Logical, TRUE if block analysis
#'
#' @return Adjusted mc to pass to internal
#' @noRd
.get_xy_from_formula <- function(mc, DA=FALSE, block=FALSE) {

  .formula_checker(mc = mc, block = block) ## check formula validity
  mf <- stats::model.frame(mc$formula) ## THANK YOU stats::model.frame *cries*
  mc$Y <- mf[[1]]
  if (isFALSE(DA)) mc$Y <- as.matrix(mc$Y)
  mc$X <- as.matrix(mf[[2]])
  mc$data <- mc$formula <- NULL 
}
## ------------------------------------------------ ##
## Given data in mc, get X and Y from assay/colData for PLS and PLSDA family
## First, we'll define a function that assumes X and Y are given
#' @title Get X and Y to pass to internal PLS(DA) functions
#'
#' @param mc A match.call
#' @param DA Logical, TRUE if discriminatory analysis
#' @param block Logical, TRUE if block.pls family
#' @importFrom MultiAssayExperiment colData assays
#' 
#' @return Adjusted mc with X and Y ready to pass to internal
#' @noRd
.get_xy <- function(mc, DA=FALSE, block=FALSE) {
  
  ## if formula given, check and get names
  if (!is.null(mc$formula)) {
    mc <- .get_xy_names_from_formula(mc= mc, block = block)
  }
  
  ## if a X is list, unlist it
  if (is(mc$X, "list")) { 
    mc$X <- unlist(mc$X)
  }
  if (!is.null(mc$indY)) { ## if indY provided
    if (!is.null(mc$Y)) {
      .stop("Only one of 'Y' and 'indY' should be provided", .subclass = "inv_signature")
    }
    mc$Y <- mc$X[indY]
    mc$X <- mc$X[-indY]
  }
  ## ensure X and Y are character
  if (!is(mc$X, "character")) {
    .stop("'X' must be a valid assay name", 
          .subclass = "inv_XY")
  }
  
  if (!is(mc$Y, "character")) {
    
    if (isTRUE(DA)) {
      .stop("'Y' must be a valid assay name", 
            .subclass = "inv_XY")
    } else {
      .stop("'Y' must be a valid name from colData(data)", 
            .subclass = "inv_XY")
    }
  }
  
  ## ensure X is valid assay(s)
  invalid_X <- mc$X[!mc$X %in% names(assays(mc$data))]
  if (length(invalid_X) != 0) {
    .stop(.subclass = "inv_XY",
          message = paste0("Not valid assay name(s): ", 
                           paste0(invalid_X, collapse = ", ")))
  } else if (isTRUE(block) & (length(mc$X) < 2)) {
    .stop(.subclass = "inv_XY",
          message = paste0("For block methods, more than 1 matrices must be provided to X"))
  }
  
  ## ensure Y is valid
  ## if !DA could be assay or colData - should turn into numeric matrix with checks
  ## if DA must only be colData - internal will check its sanity
  col.Y <- mc$Y %in% names(colData(mc$data))
  
  if (isTRUE(DA) & isFALSE(col.Y)) 
    .stop("'Y' must be a character from a valid colData(data)", .subclass = "inv_XY")
  assay.Y <- mc$Y %in% names(assays(mc$data))
  if (!(col.Y | assay.Y)) 
    .stop("'Y' must be a character from assays(data) or colData(data)", .subclass = "inv_XY")
  
  ## get Y
  if (isTRUE(DA) || col.Y) {
    mc$Y <- colData(mc$data)[, mc$Y]
  } else {
    mc$Y <- .getMatrix(assay.name = mc$Y,  mc$data)
  }
  
  ## get X - named list first
  tmp <- mc$X
  mc$X <- as.list(mc$X)
  names(mc$X) <- tmp
  mc$X <- lapply(mc$X, function(x) .getMatrix(assay.name = x, mc$data))
  ## unlist if not block
  if (isFALSE(block)) {
    mc$X <- mc$X[[1]]
  }
  mc
}

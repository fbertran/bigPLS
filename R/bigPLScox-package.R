#' bigPLS-package
#'
#' Provides Partial least squares Regression for regular, generalized linear and Cox models for big data. It allows for missing data in the explanatory variables. Repeated k-fold cross-validation of such models using various criteria. Bootstrap confidence intervals constructions are also available.
#'
#' @aliases bigPLS-package bigPLS NULL
#' 
#' @references TODO
#' 
#' @seealso [big_pls_cox()] and [big_pls_cox_gd()]
#' 
#' @examples
#' set.seed(314)
#' library(bigPLS)
#' data(sim_data)
#' head(sim_data)
#' 
"_PACKAGE"


#' @import bigSurvSGD
# #' @importFrom bigSurvSGD bigSurvSGD lambdaMaxC oneChunkC oneObsPlugingC
# #' @importClassesFrom bigSurvSGD bigSurvSGD
# #' @importMethodsFrom bigSurvSGD plot.bigSurvSGD print.bigSurvSGD
#' @importFrom graphics abline axis layout legend segments
#' @importFrom grDevices dev.new
#' @importFrom stats as.formula is.empty.model model.matrix model.response model.weights rexp runif var
# #' @importFrom utils head
#' @import bigmemory
#' @useDynLib bigPLS, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import bigalgebra
NULL

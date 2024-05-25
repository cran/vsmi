#'vsmi: Variable selection for multiple imputed data
#'
#'This is a package to implementation
#'penalized weighted least-squares estimate for variable selection on correlated multiply imputed data
#'and penalized estimating equations for generalized linear models with multiple imputation.
#'
#'@section Functions:
#'\code{\link[vsmi]{PEE}}:Penalized estimating equations for generalized linear models with multiple imputation
#'
#'\code{\link[vsmi]{PWLS}} : Penalized weighted least-squares estimate for variable selection on correlated multiply imputed data
#'
#'\code{\link[vsmi]{generate_pwls_missing_data}} : Generate example missing data for  PWLS
#'
#'\code{\link[vsmi]{generate_pee_missing_data}} : Generate example missing data for PEE
#'
#'@name vsmi
#'@docType package
#'@import Matrix
#'@import MASS
#'@import mice
#'@import qif
#'@importFrom stats binomial
#'@importFrom stats poisson
#'@importFrom stats rbinom
#'@importFrom stats rpois
#'@importFrom stats uniroot
#'@importFrom stats cor


NULL

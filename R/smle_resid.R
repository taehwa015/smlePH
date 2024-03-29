#' @importFrom stats as.formula lm predict sd quantile
#' @importFrom splines2 mSpline
#' @importFrom MASS ginv
NULL
#' Extract residuals of the full likelihood proportional hazards model
#'
#' This function extracts residuals of the full likelihood proportional hazards model estimated by the sieve estimation. Deviance-type and score-type residuals are available.
#'
#'
#' @param y survival time (> 0).
#' @param d right-censoring indicator, \code{1}: observed; \code{0}: right-censored.
#' @param x p-dimensional covariates matrix.
#' @param fit an object comes from the function \code{smle_ph}.
#' @param type type of residual, either \code{deviance} or \code{score}.
#'
#' @return \code{smle_resid} returns a numeric vector (if \code{type = "deviance"}) or a matrix (if \code{type = "score"}) of residuals extracted from the \code{object}.
#'
#' @details
#' see Halabi et al., (2024+) for detailed method explanation.
#'
#' @references
#' Halabi et al., (2024+) Sieve maximum full likelihood estimation for the proportional hazards model
#'
#'
#' @examples
#' \dontshow{
#' library(smlePH)
#' set.seed(111)
#' n = 200
#' beta = c(1, -1, 0.5, -0.5, 1)
#' p = length(beta)
#' beta = matrix(beta, ncol = 1)
#' R = matrix(c(rep(0, p^2)), ncol = p)
#' diag(R) = 1
#' mu = rep(0, p)
#' SD = rep(1, p)
#' S = R * (SD %*% t(SD))
#' x = MASS::mvrnorm(n, mu, S)
#' T = (-log(runif(n)) / (2 * exp(x %*% beta)))^(1/2)
#' C = runif(n, min = 0, max = 2.9)
#' y = apply(cbind(T,C), 1, min)
#' d = (T <= C)+0
#' ord = order(y)
#' y = y[ord]; x = x[ord,]; d = d[ord]
#' fit = smle_ph(y = y, d = d, x = x)
#' }
#' library(smlePH)
#' # The 'fit' comes from an example description of smle_ph()
#' smle_resid(y = y, d = d, x = x, fit = fit, type = "deviance")
#' smle_resid(y = y, d = d, x = x, fit = fit, type = "score")
#' @export


smle_resid = function(y,
                      d,
                      x,
                      fit,
                      type = c("score","deviance"))
{
  ord = order(y)
  ut = y = y[ord]
  d = d[ord]
  x = as.matrix(as.matrix(x)[ord,])
  beta = fit$Coef[,1]; Haz = fit$Cum.hazard[,1]
  xbeta = drop(x%*%beta)
  if (type == "score") res = ((x*(d-Haz*exp(xbeta))))
  if (type == "deviance")
  {
    id0 = which(d==0)
    id1 = which(d==1)
    dev = NULL
    tmp = Haz[id0]*exp(xbeta[id0])
    dev[id0] = sign(tmp)*sqrt(2)*sqrt(abs(tmp))
    tmp = -xbeta[id1] + log(pmax(1e-4,Haz[id1])) + Haz[id1]*exp(xbeta[id1]) -1
    dev[id1] = sign(tmp)*sqrt(2)*sqrt(abs(tmp))
    res = dev
  }
  res
}

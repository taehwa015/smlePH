#' @importFrom stats as.formula lm predict sd quantile optim pnorm
#' @importFrom splines2 mSpline
#' @importFrom MASS ginv
NULL
#' Fit the full likelihood proportional hazards model
#'
#' Fit the proportional hazards model with maximum full likelihood estimation. Sieve estimation is used for estimating the baseline hazard function.
#'
#'
#' @param y n-vector of survival time (> 0).
#' @param d n-vector of right-censoring indicator, \code{1}: observed; \code{0}: right-censored.
#' @param x p-dimensional matrix of covariates. 
#'
#' @return \code{smle_ph} returns a list containing the following components:
#' \itemize{
#'   \item \code{Coef}: regression estimator and its inferential results.
#'   \item \code{Cum.hazard}: baseline cumulative hazard function estimates.
#' }
#'
#' @details
#' see Choi et al., (2026+) for detailed method explanation.
#'
#' @references
#' Choi et al., (2026+) Residual-Based Sieve Maximum Full Likelihood Estimation for the Proportional Hazards Model
#'
#'
#' @examples
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
#' smle_ph(y = y, d = d, x = x)
#' @export

smle_ph = function(y,
                   d,
                   x)
{
  phfunc_sieve = function(y,
                          d,
                          x,
                          degree,
                          nknots)
  {
    ord = order(y)
    ut = y = y[ord]
    d = d[ord]
    x = as.matrix(as.matrix(x)[ord,])
    if (nknots == 0) {
      knots = numeric(0) 
    } else {
      qq = seq(0,1,by=1/(nknots+1))[-c(1,nknots+2)]
      knots = quantile(ut,qq)
    }
    # isp=iSpline(ut,degree=degree,knots=knots,intercept = FALSE)
    msp=mSpline(ut,degree=degree,knots=knots,intercept = FALSE)
    
    n = nrow(x)
    p = ncol(x)
    q = ncol(msp)
    slike = function(theta) {
      beta = theta[1:p]
      gamma = theta[-(1:p)]
      haz = pmax(as.vector(msp%*%gamma),1e-5)
      Haz = cumsum(haz*diff(c(0,ut)))
      xbeta=drop(x%*%beta)
      sum((d*(log(haz)+xbeta) - Haz*exp(xbeta)))
    }
    theta0 = c(rep(0,p), (rep(1,q)))
    ineqA = matrix(0,q,q)
    for (i in 1:q) {
      if (i == 1) ineqA[i,i] = 1
      else ineqA[i,c(i-1,i)] = c(-1,1)
    }
    ineqA = cbind(matrix(0, q, p), ineqA)
    ineqB = matrix(0, nrow=q)
    ineqA = cbind(matrix(0,q,p),diag(1,q))
    ineqB = matrix(0,nrow=q)
    theta = optim(theta0, slike,control = list(fnscale=-1))$par
    coef = theta[1:p]
    gamma = theta[-(1:p)]
    haz = pmax(as.vector(msp%*%gamma),1e-5)
    Haz = cumsum(haz*diff(c(0,ut)))
    list(coef = coef, like = slike(theta) ,msp = msp, isp=isp,scoef = gamma,
         Haz = Haz, haz = haz, time = y)
  }
  if (any(y<0)) y = exp(y)
  ord = order(y)
  ut = y = y[ord]
  d = d[ord]
  x = as.matrix(as.matrix(x)[ord,])
  grids = expand.grid(3, 0:(floor(nrow(x)^(1/3))))
  like_val = apply(grids, 1, function(a) phfunc_sieve(y, d, x, a[1], a[2])$like)
  val = -2*like_val + 2*rowSums(grids)*log(log(nrow(x)))
  opt = as.numeric(grids[which.min(val),])
  tmp = phfunc_sieve(y, d, x, opt[1], opt[2])
  tmp$opt = opt
  msp = tmp$msp
  exbeta=exp(drop(x%*%tmp$coef))
  a = -crossprod(x,x*tmp$Haz*exbeta)
  isp = apply(msp, 2, function(aa) cumsum(aa*diff(c(0,y))))
  g1=d*x-x*tmp$Haz*exp(drop(x%*%tmp$coef))
  g2=d*(tmp$msp)/pmax(tmp$haz,min(tmp$haz[tmp$haz>0])) - isp*exp(drop(x%*%tmp$coef))
  A11=crossprod(g1)
  A12=crossprod(g1,g2)
  A21=t(A12)
  A22=crossprod(g2)
  tmp$se = sqrt(diag(solve((A11-A12%*%solve(A22)%*%A21))))
  
  coef.smr = data.frame("Coefficients" = tmp$coef,
                        "Std. Error"= tmp$se,
                        "t statistic" = abs(tmp$coef/tmp$se),
                        "Two-sided p-value" = 2*pnorm(abs(tmp$coef/tmp$se),
                                          lower.tail = FALSE))
  c.haz = data.frame("chaz" = tmp$Haz, "time" = tmp$time)
  res = list("Coef" = round(coef.smr, 3), "Cum.hazard" = c.haz)
  res
}

# roxygen2::roxygenize()
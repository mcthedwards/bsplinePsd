#' Generate a linear, quadratic, or cubic B-spline density basis
#' @description This function generates a linear, quadratic, or cubic B-spline density basis.
#' @details \link{splineDesign} is used to generate a linear, quadratic, or cubic B-spline basis.  Each B-spline is then normalised to become a B-spline density using analytical integration.  Note that the two end knots are each coincident four times for cubic B-splines, three times for quadratic B-splines, and two times for linear B-splines.
#' @importFrom splines splineDesign
#' @importFrom Rcpp evalCpp
#' @useDynLib bsplinePsd, .registration = TRUE
#' @export
#' @param x numeric vector for which the B-spline densities are to be generated
#' @param knots knots used to generate the linear, quadratic, or cubic B-spline densities
#' @param degree integer equal to 1, 2, or 3 (default), specifying the degree of the B-splines (1 = linear, 2 = quadratic, 3 = cubic)
#' @return matrix of the B-spline density basis
#' @seealso \link{splineDesign}
#' @examples 
#' \dontrun{
#' 
#' # Generate basis functions
#' x = seq(0, 1, length = 256)
#' knots = sort(c(0, runif(10), 1))
#' basis = dbspline(x, knots)
#' 
#' # Plot basis functions
#' plot(x, basis[1, ], type = "l", ylim = c(min(basis), max(basis)))
#' for (i in 2:nrow(basis)) lines(x, basis[i, ], col = i)
#' }
dbspline = function (x, knots, degree = 3)  {
  
  if (!(degree %in% c(1, 2, 3))) stop("degree must be 1, 2, or 3 in this version of the package")
  
  knots.mult <- c(rep(knots[1], degree), knots, rep(knots[length(knots)], degree))
  
  B <- splines::splineDesign(knots.mult, x, ord = degree + 1, outer.ok = TRUE)
  
  if (degree == 3) {  # Cubic B-spline normaliser
    bs_int <- rep(NA, length = length(knots.mult) - 4)
    for (ii in 1:(length(knots.mult) - 4)) {
      part <- knots.mult[ii:(ii + 4)]
      bs_int[ii] <- AnIn1(part) + AnIn2(part) + AnIn3(part) + 
        AnIn4(part) + AnIn5(part) + AnIn6(part) + AnIn7(part) + 
        AnIn8(part)
      if (bs_int[ii] == 0) {
        bs_int[ii] <- Inf
      }
    }
  }
  
  if (degree == 2) {  # Quadratic B-spline normaliser
    bs_int <- rep(NA, length = length(knots.mult) - 3)
    for (ii in 1:(length(knots.mult) - 3)) {
      part <- knots.mult[ii:(ii + 3)]
      bs_int[ii] <- AnInQ1(part) + AnInQ2(part) + AnInQ3(part) + AnInQ4(part) 
      if (bs_int[ii] == 0) {
        bs_int[ii] <- Inf
      }
    }
  }
  
  if (degree == 1) { # Linear B-spline normaliser
    bs_int = (knots.mult[3:length(knots.mult)] - knots.mult[1:(length(knots.mult) - 2)]) / 2  # Area of triangle (height always 1)
    if (any(bs_int == 0)) bs_int[which(bs_int == 0)] = Inf  # Makes B.norm = 0 rather than NaN
  }
  
  B.norm <- t(B) / bs_int  # Normalise
  
  return(B.norm)
  
}


#' Compute unnormalised PSD using random mixture of B-splines
#' @importFrom Rcpp evalCpp
#' @useDynLib bsplinePsd, .registration = TRUE
#' @keywords internal
qpsd <- function(omega, k, v, w, u, z, degree, recompute, db.list) {
  
  #####
  # Find weights for B-spline mixture using stick-breaking
  #####
  p <- pFromV(v)
  weight <- mixtureWeight(p, w, k)
  
  if (recompute == TRUE) {
    #####
    # Find knots for B-splines using another stick-breaking
    #####
    q <- pFromV(u)  # u here
    newk <- k - degree  # CAUTION: This is the denominator for the second DP. 
    knot.diffs <- mixtureWeight(q, z, newk)
    knots <- c(0, cumsum(knot.diffs))
    
    # B-spline density matrix
    db.list <- dbspline(omega, knots, degree) 
  }
  
  #####
  # B-spline mixture
  #####
  
  # Call stored matrix for k mixtures and matrix multiply with weights
  psd <- densityMixture(weight, db.list)
  epsilon <- 1e-20 
  psd <- pmax(psd, epsilon)
  
  return(list(psd = psd,
              knots = knots,
              db.list = db.list))
  
}

#' Unnormalised log joint prior
#' @keywords internal
lprior <- function(k, v, w, u, z, tau, k.theta, 
                   MG, G0.alpha, G0.beta, 
                   MH, H0.alpha, H0.beta, 
                   tau.alpha, tau.beta) {
  
  logprior <- (MG - 1) * sum(log(1 - v)) +  # log prior for V's - beta(1, MG)
    sum((G0.alpha - 1) * log(w) + (G0.beta - 1) * log(1 - w)) +  # log prior for W's - beta(a, b)
    (MH - 1) * sum(log(1 - u)) +  # log prior for U's - beta(1, MH)
    sum((H0.alpha - 1) * log(z) + (H0.beta - 1) * log(1 - z)) -  # log prior for Z's - beta(a, b)
    k.theta * k ^ 2 -   # log prior for k
    (tau.alpha + 1) * log(tau) - tau.beta / tau  # log prior for tau (Inverse Gamma)
  
  return(logprior)
  
}

#' log Whittle likelihood
#' @keywords internal
llike <- function(omega, FZ, k, v, w, u, z, tau, pdgrm, degree, recompute, db.list) {
  
  # Calculates Whittle log-likelihood for Gaussian errors
  
  n <- length(FZ)
  ###m <- n - 2  # Hard coded for even length time series
  
  # Which boundary frequencies to remove from likelihood computation
  if (n %% 2) {  # Odd length time series
    bFreq <- 1  # Remove first
  } 
  else {  # Even length time series
    bFreq <- c(1, n)  # Remove first and last
  }
  
  # Un-normalised PSD (defined on [0, 1])
  qq.psd <- qpsd(omega, k, v, w, u, z, degree, recompute, db.list)
  q = unrollPsd(qq.psd$psd, n)  # Unrolls the unnormalised PSD to length n
  
  # Normalised PSD (defined on [0, pi])
  f <- tau * q
  
  # Whittle log-likelihood
  llike <- -sum(log(f[-bFreq]) + pdgrm[-bFreq] / (f[-bFreq] * 2 * pi)) / 2
  
  return(list(llike = llike,
              db.list = qq.psd$db.list))  
  
}

#' Unnormalised log posterior
#' @keywords internal
lpost <- function(omega, FZ, k, v, w, u, z, tau, k.theta, 
                  MG, G0.alpha, G0.beta, 
                  MH, H0.alpha, H0.beta, 
                  tau.alpha, tau.beta,
                  pdgrm, degree, recompute, db.list) {
  
  ll <- llike(omega, FZ, k, v, w, u, z, tau, pdgrm, degree, recompute, db.list)
  
  # Unnormalised log posterior
  lp <- ll$llike + 
    lprior(k, v, w, u, z, tau, k.theta, 
           MG, G0.alpha, G0.beta, 
           MH, H0.alpha, H0.beta, 
           tau.alpha, tau.beta)
  
  return(list(lp = lp,
              db.list = ll$db.list))
  
}

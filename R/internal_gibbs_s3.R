#' Plot method for psd class
#' @description This function plots the log periodogram, log posterior median PSD, and log 90\% credible region PSD.  The x-axis uses angular frequency and the y-axis is plotted on the log scale.  The PSD at the zero frequency is removed from the plot.
#' @export
#' @param x an object of class psd
#' @param ... other graphical parameters from the plot.default function
#' @return plot of the estimate of the log PSD
#' @seealso \link{gibbs_bspline}
#' @examples 
#' \dontrun{
#' 
#' # Simulate AR(4) data
#' n = 2 ^ 7
#' ar.ex = c(0.9, -0.9, 0.9, -0.9)
#' data = arima.sim(n, model = list(ar = ar.ex))
#' data = data - mean(data)
#' 
#' # Run MCMC with linear B-spline prior
#' mcmc = gibbs_bspline(data, 4000, 2000, 1, degree = 1)
#' 
#' # Plot result
#' plot(mcmc)
#' 
#' # Plot result with title
#' plot(mcmc, main = "log PSD estimate using the linear B-spline prior")
#' }
plot.psd = function(x, ...) {  # Plot method for "psd" class
  freq = seq(0, pi, length = length(x$pdgrm))
  graphics::plot.default(freq[-1], log(x$pdgrm[-1]), type = "l", col = "grey",
                         xlab = "Frequency", ylab = "log PSD", ...)
  graphics::lines(freq[-1], log(x$psd.median[-1]), lwd = 2)
  graphics::lines(freq[-1], log(x$psd.p05[-1]), lwd = 2, lty = 2, col = 4)
  graphics::lines(freq[-1], log(x$psd.p95[-1]), lwd = 2, lty = 2, col = 4)
  graphics::legend("topright", legend = c("periodogram", "posterior median", "90% credible region"), 
                   col = c("grey", "black", "blue"), lwd = c(1, 2, 2), lty = c(1, 1, 2))
}
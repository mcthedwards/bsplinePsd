## bsplinePsd 0.4.1

Previous versions only allowed the user to use cubic B-spline densities.  This version allows the user to choose between linear, quadratic, or cubic B-spline densities.  The user can specify degree = 1 for linear B-spline densities, degree = 2 for quadratic B-spline densities, and degree = 3 (default) for cubic B-spline densities.

The function gibbs_bspline can now handle odd length time series.

An S3 plot method has been included so the user can easily plot their PSD estimate.

## bsplinePsd 0.2.0

Added an argument called k1 in the gibbs_bspline function.  This allows the user to specify the starting value for parameter k.  If well-chosen, this can speed up convergence significantly.  The default is set to 20, which works well on all of the cases I have come across.  If missing (NA), then a random integer between 5 and kmax will be selected as the starting value for k.

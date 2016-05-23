# Armijo Rule (or Sufficient Decrease Condition)
#
# @param phi_0 Value of phi at alpha = 0
# @param dphi_0 the directional derivative at alpha = 0
# @param alpha the step length.
# @param phi_a the value of phi at \code{alpha}.
# @param c1 the sufficient decrease constant.
# @return \code{TRUE} if the step \code{alpha} represents a sufficient decrease.
armijo_ok <- function(phi_0, dphi_0, alpha, phi_a, c1) {
  phi_a <= phi_0 + c1 * alpha * dphi_0
}

# Curvature Condition
#
# Line search test.
#
# Ensures that the directional derivative of the line search direction at a
# candidate step size is greater than a specified fraction of the slope of the
# line at the starting point of the search. This condition is used to stop step
# sizes being too small.
#
# In combination with the sufficient decrease condition \code{\link{armijo_ok}}
# these conditions make up the Wolfe conditions.
#
# @param d0 Directional derivative at starting point.
# @param da Directrional derivative at step alpha.
# @param c2 Curvature condition constant. Should take a value between \code{c1}
#  (the constant used in the sufficient decrease condition check) and 1.
# @return \code{TRUE} if the curvature condition is met.
curvature_ok <- function(d0, da, c2) {
  da > c2 * d0
}


curvature_oks <- function(step0, step_a, c2) {
  curvature_ok(step0$d, step_a$d, c2)
}

# Strong Curvature Condition
#
# Line search test.
#
# Ensures that the value of the directional derivative of the line search
# direction at a candidate step size is equal to or greater than a specified
# fraction of the slope of the line at the starting point of the search, while
# having the same direction. This condition is used to make the step size lie
# close to a stationary point. Unlike the normal curvature condition, a step
# size where the sign of the gradient changed (e.g. the minimum had been
# skipped) would not be acceptable for the strong curvature condition.
#
# In combination with the sufficient decrease condition \code{\link{armijo_ok}}
# these conditions make up the Strong Wolfe conditions.
#
# @param d0 Directional derivative at starting point.
# @param da Directrional derivative at step alpha.
# @param c2 Curvature condition constant. Should take a value between \code{c1}
#  (the constant used in the sufficient decrease condition check) and 1.
# @return \code{TRUE} if the strong curvature condition is met.
strong_curvature_ok <- function(d0, da, c2) {
  abs(da) <= -c2 * d0
}

strong_wolfe_ok <- function(f0, d0, alpha, fa, da, c1, c2) {
  armijo_ok(f0, d0, alpha, fa, c1) &&
    strong_curvature_ok(d0, da, c2)
}

strong_wolfe_oks <- function(step0, step, c1, c2) {
  armijo_ok(step0$f, step0$d, step$alpha, step$f, c1) &&
    strong_curvature_ok(step0$d, step$d, c2)
}

armijo_oks <- function(step0, step_a, c1) {
  armijo_ok(step0$f, step0$d, step_a$alpha, step_a$f, c1)
}

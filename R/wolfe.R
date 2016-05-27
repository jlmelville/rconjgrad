# Armijo Rule (or Sufficient Decrease Condition)
#
# Line search test.
#
# The sufficient decrease condition is met if the line search step length yields
# a decrease in the function value that is sufficiently large (relative to the
# size of the step).
#
# This test prevents large step sizes that, while representing a function value
# decrease, don't reduce it by very much, which could indicate that the
# function minimum has been stepped over and you're now going back up the slope.
# Also, this condition can always be met by taking a sufficiently small step,
# so line searches involving this condition can always terminate. The downside
# is that you can end up taking very small steps, so it's usual to combine this
# condition with one that encourages larger step sizes.
#
# @param f0 Value of function at starting point of line search.
# @param d0 Directional derivative at starting point of line search.
# @param alpha the step length.
# @param fa Value of function at alpha.
# @param c1 the sufficient decrease constant. Should take a value between 0 and
#   1.
# @return \code{TRUE} if the step \code{alpha} represents a sufficient decrease.
armijo_ok <- function(f0, d0, alpha, fa, c1) {
  fa <= f0 + c1 * alpha * d0
}

# Armijo Rule (or Sufficient Decrease Condition)
#
# Line search test.
#
# The sufficient decrease condition is met if the line search step length yields
# a decrease in the function value that is sufficiently large (relative to the
# size of the step).
#
# This test prevents large step sizes that, while representing a function value
# decrease, don't reduce it by very much, which could indicate that the
# function minimum has been stepped over and you're now going back up the slope.
# Also, this condition can always be met by taking a sufficiently small step,
# so line searches involving this condition can always terminate. The downside
# is that you can end up taking very small steps, so it's usual to combine this
# condition with one that encourages larger step sizes.
#
# @param step0 Line search values at starting point.
# @param step Line search value at a step along the line.
# @param c1 the sufficient decrease constant. Should take a value between 0 and
#   1.
# @return \code{TRUE} if the step represents a sufficient decrease.
armijo_ok_step <- function(step0, step, c1) {
  armijo_ok(step0$f, step0$d, step$alpha, step$f, c1)
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
# @param da Directional derivative at step alpha.
# @param c2 Curvature condition constant. Should take a value between \code{c1}
#  (the constant used in the sufficient decrease condition check) and 1.
# @return \code{TRUE} if the curvature condition is met.
curvature_ok <- function(d0, da, c2) {
  da > c2 * d0
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
# @param step0 Line search values at starting point.
# @param step Line search value at a step along the line.
# @param c2 Curvature condition constant. Should take a value between \code{c1}
#  (the constant used in the sufficient decrease condition check) and 1.
# @return \code{TRUE} if the curvature condition is met.
curvature_ok_step <- function(step0, step, c2) {
  curvature_ok(step0$d, step$d, c2)
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
# @param step0 Line search values at starting point.
# @param step Line search value at a step along the line.
# @param c2 Curvature condition constant. Should take a value between \code{c1}
#  (the constant used in the sufficient decrease condition check) and 1.
# @return \code{TRUE} if the curvature condition is met.
strong_curvature_ok_step <- function(step0, step, c2) {
  strong_curvature_ok(step0$d, step$d, c2)
}

# Are the Strong Wolfe Conditions Met?
#
# Step size check.
#
# Returns true if the Strong Wolfe conditions are met, consisting of the
# sufficient decrease conditions and the strong curvature condition.
#
# @param f0 Function value at starting point.
# @param d0 Directional derivative value at starting point.
# @param alpha Step length.
# @param fa Function value at alpha.
# @param da Directional derivative at alpha.
# @param c1 Constant used in sufficient decrease condition. Should take a value
#   between 0 and 1.
# @param c2 Constant used in curvature condition. Should take a value between
#   c1 and 1.
# @return TRUE if the Strong Wolfe condition is met by the step size.
strong_wolfe_ok <- function(f0, d0, alpha, fa, da, c1, c2) {
  armijo_ok(f0, d0, alpha, fa, c1) &&
    strong_curvature_ok(d0, da, c2)
}

# Are the Strong Wolfe Conditions Met?
#
# Line search test.
#
# Returns true if the candidate step size meets the Strong Wolfe conditions,
# consisting of the sufficient decrease conditions and the strong curvature
# condition.
#
# @param step0 Line search values at starting point of line search.
# @param step Line search value at candiate step size.
# @param c1 Constant used in sufficient decrease condition. Should take a value
#   between 0 and 1.
# @param c2 Constant used in curvature condition. Should take a value between
#   c1 and 1.
# @return TRUE if the Strong Wolfe condition is met by the step size.
strong_wolfe_ok_step <- function(step0, step, c1, c2) {
  armijo_ok_step(step0, step, c1) && strong_curvature_ok_step(step0, step, c2)
}


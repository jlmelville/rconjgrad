# Dot Product
#
# Returns the dot product of two vectors.
#
# @param a Numeric vector.
# @param b Numeric vector of the same dimensions as a.
# @return Dot product of \code{a} and \code{b}.
dot <- function(a, b) {
  sum(a * b)
}

# Prints out a vector. For debugging purposes.
format_vec <- function(vec) {
  paste(formatC(vec), collapse = ' ')
}

# Construct Line Function
#
# Factory function for line function. Returns a 1D function of alpha, the step
# size, the return value of which is a list containing the function value,
# directional derivative, and other values calculated for that distance
# along the search direction.
#
# @param fn Function.
# @param gr Gradient.
# @param par Vector of parameters.
# @param debug If TRUE, log information about the line search values.
# @param ... Other parameters to pass to fn and gr when they are invoked.
# @return Line function.
make_phi <- function(fn, gr, par, pv, debug = FALSE, ...) {
  function(alpha) {
    xa <- par + alpha * pv
    f <- fn(xa, ...)
    df <- gr(xa, ...)
    d <- dot(df, pv)
    if (debug) {
      message(
        " p = ", format_vec(pv),
        " alpha = " , formatC(alpha),
        " f = ", formatC(f),
        " df = ", format_vec(df),
        " d = ", formatC(d))
    }
    list(f = f, df = df, d = d, alpha = alpha, par = xa)
  }
}

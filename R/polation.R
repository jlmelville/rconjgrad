cubic_extrapolate <- function(x1, f1, g1, x2, f2, g2) {
  A <- 6 * (f1 - f2) + 3 * (g2 + g1) * (x2 - x1)
  B <- 3 * (f2 - f1) - (2 * g1 + g2) * (x2 - x1)
  suppressWarnings(
    x1 - g1 * (x2 - x1) ^ 2 / (B + sqrt(B * B - A * g1 * (x2 - x1)))
  )
}

cubic_extrapolate_step <- function(step1, step2) {
  cubic_extrapolate(step1$alpha, step1$f, step1$d, step2$alpha, step2$f, step2$d)
}

cubic_interpolate <- function(x1, f1, g1, x2, f2, g2) {
  nwc(x1, f1, g1, x2, f2, g2)
  # A <- 6 * (f1 - f2) / (x2 - x1) + 3 * (g2 + g1)
  # B <- 3 * (f2 - f1) - (2 * g1 + g2) * (x2 - x1)
  # # num. error possible, ok!
  # suppressWarnings(
  #   x1 + (sqrt(B * B - A * g1 * (x2 - x1) ^ 2) -  B) / A
  # )
#  A <- 6 * (f1 - f2) + 3 * (g2 + g1) * (x2 - x1)
#  B <- 3 * (f2 - f1) - (2 * g1 + g2) * (x2 - x1)
#  x1 - g1 * (x2 - x1) ^ 2 / (B + sqrt(B * B - A * g1 * (x2 - x1)))
}

nwc <- function(x0, f0, df0, x1, f1, df1, ignoreWarnings = FALSE) {
  d1 <- df0 + df1 - 3 * ((f0 - f1) / (x0 - x1))

  if (ignoreWarnings) {
    suppressWarnings(
      d2 <- sign(x1 - x0) * sqrt(d1 * d1 - df0 * df1)
    )
  }
  else {
    d2 <- sign(x1 - x0) * sqrt(d1 * d1 - df0 * df1)
  }
  x1 - (x1 - x0) * ((df1 + d2 - d1) / (df1 - df0 + 2 * d2))

}

cubic_interpolate_step <- function(step1, step2) {
  cubic_interpolate(step1$alpha, step1$f, step1$d,
                    step2$alpha, step2$f, step2$d)
}

quadratic_interpolate <- function(x1, f1, g1, x2, f2) {
  x1 - (0.5 * g1 * (x2 - x1) ^ 2) / (f2 - f1 - g1 * (x2 - x1))
}

quadratic_interpolate_step <- function(step1, step2) {
  quadratic_interpolate(step1$alpha, step1$f, step1$d,
                        step2$alpha, step2$f)
}

mtq <- function(stx, fx, dx, stp, fp) {
  stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2) * (stp - stx)
}

mts <- function(stx, dx, stp, dp) {
  stp + (dp / (dp - dx)) * (stx - stp)
}

quadratic_interpolateg <- function(x1, g1, x2, g2) {
  #x1 - (0.5 * g1 * (x2 - x1) ^ 2) / (f2 - f1 - g1 * (x2 - x1))
  x2 + (x1 - x2) * g2 / (g2 - g1)
}

#' Tweak Extrapolated Point
#'
#' Prevents the extrapolated point from being too far away from or to close to
#' the points used in the extrapolation.
#'
#' @param xnew 1D position of the new point.
#' @param x1 1D position of the first points used in the extrapolation.
#' @param x2 1D position of the second point used in the extrapolation.
#' @param EXT Maximum multiple of \code{x2} that \code{xnew} is allowed to be
#'  extrapolated to.
#' @param INT Given the distance between \code{x1} and \code{x2}, specified what
#'  multiple of that distance is the minimum allowed distance for \code{xnew}
#'  from \code{x2}.
#' @return A value of \code{xnew} that obeys the minimum and maximum distance
#'  constraints from \code{x2}.
tweak_extrapolation <- function(xnew, x1, x2, EXT, INT) {
  # num prob | wrong sign?
  if (!is.double(xnew) || is.nan(xnew) || is.infinite(xnew) || xnew < 0) {
    # extrapolate maximum amount
    xnew <- x2 * EXT
  } else if (xnew > (x2 * EXT)) {
    # new point beyond extrapolation limit?
    # extrapolate maximum amount
    xnew <- x2 * EXT
  } else if (xnew < (x2 + INT * (x2 - x1))) {
    # new point too close to previous point?
    xnew <- x2 + INT * (x2 - x1)
  }
  xnew
}

#' Tweak Interpolated Point
#'
#' Prevents interpolated point from getting too close to either of the
#' points used for the interpolation. If the point is not a number or infinite,
#' then it is set to the bisection of the position of the two interpolating
#' points before the check for a too-close approach is carried out.
#'
#' @param xnew Position of the interpolated point.
#' @param x1 Position of the first point used for interpolation.
#' @param x2 Position of the second point used for interpolation.
#' @param INT Given the distance between \code{x1} and \code{x2}, specifies what
#'  multiple of that distance is the minimum allowed distance for \code{xnew}
#'  from \code{x1} or \code{x2}.
#' @return Tweaked position of \code{xnew} such that it is not too close to
#' either \code{x1} or \code{x2}.
tweak_interpolation <- function(xnew, x1, x2, INT) {
  if (is.nan(xnew) || is.infinite(xnew)) {
    # if we had a numerical problem then bisect
    xnew <- (x1 + x2) / 2
  }
  # don't accept too close
  max(min(xnew, x2 - INT * (x2 - x1)), x1 + INT * (x2 - x1))
}

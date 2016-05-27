dot <- function(a, b) {
  sum(t(a) %*% b)
}

format_vec <- function(vec) {
  paste(formatC(vec), collapse = ' ')
}

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
    list(f = f, df = df, d = d, alpha = alpha, x = xa)
  }
}

make_step0 <- function(fn, gr, x, pv, f = fn(x), df = gr(x)) {
  list(
    x = x,
    alpha = 0,
    f = f,
    df = df,
    d = dot(pv, df)
  )
}


# max_alpha <- function(fcn, am, mu) {
#   v <- fcn(0)
#   (v$f - fcn(am)$f)/(-mu *  v$g)
# }

# fd <- function(f, d = 1e-6) {
#   function(par) {
#     (f(par + d) - f(par - d)) / (2 * d)
#   }
# }

# mvec <- function(vec) {
#   paste(vec, collapse = '; ')
# }

#' Conjugate Gradient Minimization
#'
#' Optimization method.
#
#' Conjugate gradient optimimzation routine. A translation and slight
#' modification of the Matlab code \code{minimize.m} by
#' \href{http://learning.eng.cam.ac.uk/carl/code/minimize/}{Carl Edward Rasmussen}.
#'
#' @details
#' The \code{line_search} parameter takes one of two string values indicating
#' the type of line search function desired: either the original line search
#' from the Matlab code, or the modified More'-Thuente line search algorithm
#' implemented in MINPACK.
#'
#' Alternatively, you may provide your own line search function. Doing so
#' requires some knowledge of internal interfaces, which is described in the
#' package help documentation which can be accessed via the following help
#' command:
#'
#' \code{package?rconjgrad}
#'
#' Factory functions that generate a suitable line search function for the
#' existing line search methods are \code{\link{more_thuente}} and
#' \code{\link{rasmussen}}.
#'
#' @param par Vector of initial numeric parameters.
#' @param fn Function to optimize. Should take a vector of the length of
#'  \code{par} and return a scalar numeric value.
#' @param gr Gradient of \code{fn}. Should take a vector of the length of
#'  \code{par} and return a vector of the gradients with respect to each
#'  parameter.
#' @param c1 Constant used in sufficient decrease condition. Should take a value
#'   between 0 and 1.
#' @param c2 Constant used in curvature condition. Should take a value between
#'   \code{c1} and 1.
#' @param max_iter Maximum number of iterations to carry out the optimization.
#' @param max_line_fn Maximum number of evaluations of \code{fn} to carry out
#'  during each line search.
#' @param red Scalar to determine initial step size guess.
#' @param max_alpha_mult Maximum scale factor to use when guessing the initial
#'  step size for the next iteration.
#' @param abstol Absolute tolerance. If not \code{NULL}, then optimization will
#'  stop early if the return value of \code{fn} falls below this value.
#' @param reltol Relative tolerance. If not \code{NULL}, then optimization will
#'  stop early if the relative decrease in the value of \code{fn} on successive
#'  iterations falls below this value.
#' @param line_search Line search type. Can be one of
#'  \itemize{
#'    \item \code{"r"} The Rasmussen line search as originally implemented in
#'    the original Matlab code.
#'    \item \code{"mt"} More'-Thuente line search as originally implemented in
#'    MINPACK.
#'  }
#' You may also assign a line search function directly to this parameter. See
#' 'Details' for more information.
#' @param ortho_restart If \code{TRUE}, then if successive conjugate gradient
#'  directions are not sufficiently orthogonal, reset the search direction to
#'  steepest descent.
#' @param nu If the dot product of the old and new conjugate gradient direction
#'  (normalized with respect to inner product of the new direction) exceeds
#'  this value, then the two directions are considered non-orthogonal and the
#'  search direction is reset to steepest descent. Only used if
#'  \code{ortho_restart} is \code{TRUE}.
#' @param prplus If \code{TRUE} then the 'PR+' variant of the Polak-Ribiere
#'  update will be used: when the beta scale factor used to calculate the
#'  new direction is negative, the search direction will be reset to
#'  steepest descent.
#' @param eps Epsilon for avoiding numerical issues.
#' @param verbose If \code{TRUE} log information about the status of the
#'  optimization at each iteration.
#' @param debug If \code{TRUE} logs \emph{lots} of information about the status
#'  of the optimization.
#' @param ... Other parameters to pass to the \code{fn} and \code{gr} functions.
#' @return List containing:
#' \itemize{
#'  \item \code{par} Optimized parameters.
#'  \item \code{value} Return value of \code{fn} for the optimized parameters.
#'  \item \code{values} Vector of return values of \code{fn} at each iteration
#'    of the optimization.
#'  \item \code{iter} Number of iterations optimization took place over.
#'  \item \code{counts} Sublist of two values, giving the number of evaluations
#'    of \code{fn} and \code{gr}, respectively.
#'  \item \code{nresets} Number of times the optimizer was reset to steepest
#'    descent.
#' }
#' @export
#' @examples
#' # The venerable Rosenbrock Banana function
#' rosenbrock_banana <- list(
#' fr = function(x) {
#'  x1 <- x[1]
#'  x2 <- x[2]
#'  100 * (x2 - x1 * x1) ^ 2 + (1 - x1) ^ 2
#' },
#' grr = function(x) {
#'  x1 <- x[1]
#'  x2 <- x[2]
#' c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
#'    200 *      (x2 - x1 * x1))
#' })
#'
#' # Default is to use Rasmussen line search with c1 = 0.05 and c2 = 0.1
#' res <- conj_grad(par = c(-1.2, 1),
#'                  fn = rosenbrock_banana$fr, gr = rosenbrock_banana$grr)
#'
#' # Turning down eps, abstol and reltol to compare with Matlab result at
#' # http://learning.eng.cam.ac.uk/carl/code/minimize/
#' # But must also put an upper limit on the number of line function evaluations
#' # because numerical errors prevent convergence (so don't do this normally!)
#' res <- conj_grad(par = c(0, 0),
#'                  fn = rosenbrock_banana$fr, gr = rosenbrock_banana$grr,
#'                  eps = .Machine$double.xmin, reltol = .Machine$double.xmin,
#'                  abstol = .Machine$double.xmin, max_line_fn = 20)
#' res$par # c(1, 1)
#' res$value # 1.232595e-32
#' res$iter # 19 iterations
#' res$counts # c(79, 79) 79 fn and 79 gr evaluations
#'
#' # Use More'-Thuente line search with typical CG Wolfe parameters mentioned
#' # in Nocedal and Wright's book on numerical optimization.
#' res <- conj_grad(par = c(-1.2, 1),
#'                  fn = rosenbrock_banana$fr, gr = rosenbrock_banana$grr,
#'                  line_search = "mt",
#'                  c1 = 1e-4, c2 = 0.1)
#' \dontrun{
#' # Can pass a function to line_search if you want to write your own
#' # line search function. This example is the same as the previous one, but
#' # uses the More-Thuente factory function.
#' # Yes, you do have to specify c1 and c2 in two separate places. Sorry.
#' res <- conj_grad(par = c(-1.2, 1),
#'                  fn = rosenbrock_banana$fr, gr = rosenbrock_banana$grr,
#'                  line_search = more_thuente(c1 = 1e-4, c2 = 0.1),
#'                  c1 = 1e-4, c2 = 0.1)
#' }
conj_grad <- function(par, fn, gr,
                      c1 = c2 / 2,
                      c2 = 0.1,
                      max_iter = 1000,
                      max_line_fn = Inf,
                      red = 1,
                      max_alpha_mult = 10,
                      abstol = .Machine$double.eps,
                      reltol = sqrt(.Machine$double.eps),
                      line_search = "r",
                      ortho_restart = FALSE, nu = 0.1,
                      prplus = FALSE,
                      eps = .Machine$double.eps,
                      verbose = FALSE, debug = FALSE,
                      ...) {

  if (class(line_search) == "character") {
    line_search <- tolower(line_search)
    if (line_search == "mt") {
      line_search <- more_thuente(c1 = c1, c2 = c2, max_fn = max_line_fn)
    }
    else if (line_search == "r") {
      line_search <- rasmussen(c1 = c1, c2 = c2, max_fn = max_line_fn)
    }
    else {
      stop("Unknown line_search type '", line_search, "'")
    }
  }
  else if (class(line_search) != "function") {
    stop("line_search parameter must be either a valid string or a function")
  }
  nresets <- 0
  nfn <- 0
  ngr <- 0
  # calculate function value and gradient at initial location
  f0 <- fn(par, ...)
  if (is.nan(f0)) {
    f0 <- Inf
  }
  df0 <- gr(par, ...)
  nfn <- nfn + 1
  ngr <- ngr + 1
  # pv is the descent direction (steepest initially)
  # d0 is the directional derivative (phi')
  pv <- -df0
  d0 <- dot(df0, pv)
  step0 <- list(alpha = 0, f = f0, d = d0, df = df0)

  fX <- step0$f

  # initial step is red/(||pv||+1)
  alpha <- red / (1 - d0)

  iter <- 0
  while (iter < max_iter) {
    iter <- iter + 1

    phi <- make_phi(fn, gr, par, pv, debug = debug, ...)

    ls_result <- line_search(phi, step0, alpha)

    step <- ls_result$step
    nfn <- nfn + ls_result$nfn
    ngr <- ngr + ls_result$ngr

    if (!strong_wolfe_ok_step(step0, step, c1, c2)) {
      if (verbose) {
        message("Could not satisfy the Strong Wolfe Conditions, exiting")
      }
      break
    }

    # update variables
    par <- step$par
    # update costs
    fX <- c(fX, step$f)

    if (!is.null(abstol) && step$f < abstol) {
      if (verbose) {
        message("Reached absolute tolerance, exiting")
      }
      break
    }

    rtol <- reltol * (abs(step$f) + reltol)
    reduction <- step$f / (step0$f + eps)
    if (reduction < rtol) {
      if (verbose) {
        message("Reached relative tolerance, exiting")
      }
      break
    }

    # Polack-Ribiere CG direction
    beta <- pr_update(step0, step, eps)
    pv <- (beta * pv) - step$df
    ortho_test <- abs(dot(step$df, step0$df)) / dot(step$df, step$df)

    # update old values
    d_old <- step0$d
    step0 <- step
    step0$alpha <- 0
    step0$d <- dot(step0$df, pv)

    # reset CG if new direction is not a descent direction
    reset_CG <- FALSE
    if (step0$d > 0) {
      if (verbose) {
        message("New CG direction is not a descent direction, ",
                  "resetting to steepest descent")
      }
      reset_CG <- TRUE
    }
    else if (ortho_restart && ortho_test >= nu) {
      if (verbose) {
        message("New CG direction is not sufficiently orthogonal, ",
                "resetting to steepest descent")
      }
      reset_CG <- TRUE
    }
    else if (prplus && beta < 0) {
      if (verbose) {
        message("PR beta is -ve, resetting to steepest descent")
      }
      reset_CG <- TRUE
    }
    if (reset_CG) {
      nresets <- nresets + 1
      pv <- -step0$df
      step0$d <- dot(step0$df, pv)
    }

    # new initial step size
    slope_ratio <- d_old / (step0$d - eps)
    old_alpha <- alpha
    alpha <- step$alpha * min(max_alpha_mult, slope_ratio)

    if (verbose) {
      message("iter: ", iter,
              " phi(a) = ", formatC(step$f),
              " phi'(a) = ", formatC(step$d),
              " alpha_init = ", formatC(old_alpha),
              " alpha = ", formatC(step$alpha),
              " slope_ratio = ", formatC(slope_ratio), " ",
                                 formatC(d_old), " / ", formatC(step0$d),
              " ortho = ", formatC(ortho_test),
              " beta = ", formatC(beta)
      )
    }
  }

  list(par = par, value = fX[length(fX)], values = fX, iter = iter,
       counts = c(nfn, ngr), nresets = nresets)
}

# Polak-Ribiere Update
#
# Scale factor applied to the previous conjugate gradient direction.
#
# Uses the formula of Polak and Ribiere. Generally considered a better choice
# than the Fletcher-Reeves method.
#
# @param step0 Line search values at starting point of line search.
# @param step Line search value at current step size.
# @param eps Epsilon for avoiding numerical issues.
# @return Beta parameter for conjugate gradient update.
pr_update <- function(step0, step, eps = .Machine$double.eps) {
  (dot(step$df, step$df) - dot(step$df, step0$df)) /
    (dot(step0$df, step0$df) + eps)
}


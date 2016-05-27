conj_grad <- function(par, fn, gr,
                      c1 = c2 / 2,
                      c2 = 0.1,
                      max_iter = 100,
                      red = 1,
                      max_alpha_mult = 10,
                      eps = .Machine$double.eps,
                      abstol = NULL,
                      reltol = sqrt(.Machine$double.eps),
                      verbose = FALSE, debug = FALSE,
                      line_search = rasmussen(c1 = c1, c2 = c2),
                      ortho_restart = FALSE, nu = 0.1,
                      prplus = FALSE,
                      ...) {
  nfn <- 0

  # calculate function value and gradient at initial location
  f0 <- fn(par, ...)
  if (is.nan(f0)) {
    f0 <- Inf
  }
  df0 <- gr(par, ...)
  nfn <- nfn + 1
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

    if (!strong_wolfe_ok_step(step0, step, c1, c2)) {
      if (verbose) {
        message("Could not satisfy the Strong Wolfe Conditions, exiting")
      }
      break
    }

    # update variables
    par <- par + step$alpha * pv
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
    #f_old <- step0$f
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
       counts = c(nfn, nfn))
}

pr_update <- function(step0, step, eps = .Machine$double.eps) {
  (dot(step$df, step$df) - dot(step$df, step0$df)) / (dot(step0$df, step0$df) + eps)
}


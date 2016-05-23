find_sane <- function(phi_alpha, alpha, min_alpha = 0, max_fn = 20) {
  nfn <- 0
  while (nfn < max_fn) {
    step <- phi_alpha(alpha)
    if (is.nan(step$f)) {
      step$f <- Inf
    }
    nfn <- nfn + 1

    if (is.infinite(step$f) || any(is.nan(step$df) + is.infinite(step$df))) {
      # bisect and try again
      alpha <- (min_alpha + alpha) / 2
    } else {
      break
    }
  }
  list(step = step, nfn = nfn)
}

extrapolate_step_size <- function(phi_alpha, alpha, step0,
                                  max_fn = 20,
                                  c1, c2, EXT, INT) {
  step <- list(alpha = alpha)

  # keep a count of number of function evaluations
  nfn <- 0
  while (1) {
    sane_result <- find_sane(phi_alpha, step$alpha, max_fn, min_alpha = 0)
    nfn <- nfn + sane_result$nfn
    max_fn <- max_fn - sane_result$nfn
    step <- sane_result$step

    # are we done extrapolating?
    if (extrapolation_ok(step0, step, c1, c2) || max_fn == 0) {
      break
    }
    step$alpha <- tweaked_extrapolation(step0, step, EXT, INT)
    #    message("ext_alpha = ", formatC(step$alpha))
  }

  list(step = step, nfn = nfn)
}

extrapolation_ok <- function(step0, step_a, c1, c2) {
  #  message("curvature ok?")
  #  if (curvature_oks(step0, step_a, c2)) {
  #   message("Curvature is ok, stopping extrapolation")
  #  }
  #  message("armijo ok?")
  # if (!armijo_oks(step0, step_a, c1)) {
  #   message("Armijo is not ok, stopping extrapolation")
  # }
  curvature_oks(step0, step_a, c2) || !armijo_oks(step0, step_a, c1)
}

tweaked_extrapolation <- function(step0, step3, EXT, INT) {
  ext_alpha <- cubic_extrapolates(step0, step3)
  tweak_extrapolation(ext_alpha, step0$alpha, step3$alpha, EXT, INT)
}

interpolate_step_size <- function(phi_alpha, step0, step3,
                                  max_fn = 20,
                                  c1, c2, INT) {
  step2 <- step0
  nfn <- 0
  # keep interpolating
  while (!strong_wolfe_oks(step0, step3, c1, c2) && nfn < max_fn) {
    # choose subinterval
    if (step3$d > 0 || !armijo_oks(step0, step3, c1)) {
      #      message("step4<-3")
      step4 <- step3
    } else {
      #      message("step2<-3")
      step2 <- step3
    }

    #    message("s2 f = ", formatC(step2$f), " d = ", formatC(step2$d),
    #            " s4 f = ", formatC(step4$f), " d = ", formatC(step4$d))
    if (step4$f > step0$f) {
      #      message("quadratic interpolate")
      step3$alpha <- quadratic_interpolates(step2, step4)
    } else {
      #      message("cubic interpolate")
      step3$alpha <- cubic_interpolates(step2, step4)
    }

    #    message("pre-tweak = ", formatC(step3$alpha),
    #            " a2 = ", formatC(step2$alpha),
    #            " a4 = ", formatC(step4$alpha))
    step3$alpha <- tweak_interpolation(step3$alpha, step2$alpha, step4$alpha, INT)
    #   message("post-tweak = ", formatC(step3$alpha))

    step3 <- phi_alpha(step3$alpha)
    nfn <- nfn + 1
  }
  list(step = step3, nfn = nfn)
}


ras_ls <- function(phi_alpha, alpha, step0, max_fn = 20,
                   c1 = 0.1, c2 = 0.1 / 2, EXT = 3.0, INT = 0.1) {
  #  message("initial alpha = ", formatC(alpha))

  if (c2 < c1) {
    stop("Rasmussen line search: c2 < c1")
  }

  #message("Extrapolating")
  # extrapolate from initial alpha until either curvature condition is met
  # or the armijo condition is NOT met
  ex_result <- extrapolate_step_size(phi_alpha, alpha, step0,
                                     max_fn,
                                     c1, c2, EXT, INT)

  step3 <- ex_result$step
  nfn <- ex_result$nfn
  max_fn <- max_fn - nfn
  if (max_fn <= 0) {
    return(ex_result)
  }

  #  message("alpha ext = ", formatC(step3$alpha))
  # interpolate until the Strong Wolfe conditions are met
  int_result <- interpolate_step_size(phi_alpha, step0, step3, max_fn,
                                      c1, c2, INT)
  int_result$nfn <- int_result$nfn + nfn
  int_result
}


rasmussen <- function(c1 = c2 / 2, c2 = 0.1,
                      INT = 0.1, EXT = 3.0, max_fn = 20) {
  if (c2 < c1) {
    stop("rasmussen line search: c2 < c1")
  }
  function(phi, step0, alpha) {
    ras_ls(phi, alpha, step0, max_fn, c1, c2, EXT, INT)
  }
}


library(rcgmin)
context("CG with More'-Thuente")

fxs <- c(2.420000e+01, 4.128256e+00, 3.837450e+00, 3.366653e+00, 1.495371e+00,
         1.461521e+00, 1.199143e+00, 9.773584e-01, 7.542285e-01, 5.737405e-01,
         5.163744e-01, 3.001563e-01, 2.117800e-01, 1.957472e-01, 4.448084e-02,
         3.848279e-02, 1.955651e-02, 6.787355e-03, 2.485170e-03, 2.934491e-04,
         5.966681e-05, 1.974353e-07, 3.267785e-09, 7.230753e-11, 4.037467e-11,
         2.588282e-11, 2.110310e-20, 2.045993e-20, 1.817153e-27, 1.874777e-29)

pr_result <- conj_grad(par = c(-1.2, 1), fn = rosenbrock_banana$fr,
                       gr = rosenbrock_banana$grr, eps = .Machine$double.xmin,
                       reltol = .Machine$double.xmin,
                       abstol = .Machine$double.xmin,
                       line_search = more_thuente(c1 = 0.05, c2 = 0.1))
expect_equal(pr_result$value, 1.874777e-29, tolerance = 1e-7)
expect_equal(pr_result$values, fxs, tolerance = 1e-7)
expect_equal(pr_result$iter, 30)
expect_equal(pr_result$counts, c(150, 150))

fxs00 <-
  c(1.000000e+00, 7.713628e-01, 5.946036e-01, 4.247294e-01, 3.067932e-01,
    2.746692e-01, 1.454919e-01, 9.938971e-02, 8.945326e-02, 6.868146e-03,
    5.527433e-03, 5.856635e-04, 1.341593e-04, 3.954782e-06, 2.753610e-07,
    2.641489e-08, 3.203393e-11, 3.594079e-13, 3.323960e-15, 7.501299e-21,
    2.320237e-28, 3.081488e-31)


pr_result00 <- conj_grad(par = c(0, 0), fn = rosenbrock_banana$fr,
                         gr = rosenbrock_banana$grr, eps = .Machine$double.xmin,
                         line_search = more_thuente(c1 = 0.05, c2 = 0.1),
                         reltol = .Machine$double.xmin,
                         abstol = .Machine$double.xmin)
expect_equal(pr_result00$value, 3.081488e-31, tolerance = 1e-7)
expect_equal(pr_result00$values, fxs00, tolerance = 1e-7)
expect_equal(pr_result00$iter, 22)
expect_equal(pr_result00$counts, c(122, 122))

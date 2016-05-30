library(rcgmin)
context("CG with Rasmussen")

fxs <- c(2.420000e+01, 4.128200e+00, 3.840252e+00, 3.517264e+00, 1.748701e+00,
         1.677741e+00, 1.390890e+00, 1.169249e+00, 8.912086e-01, 6.928948e-01,
         5.983546e-01, 1.050386e-01, 1.004605e-01, 5.656690e-02, 1.195812e-03,
         9.640880e-04, 3.146628e-05, 1.154634e-06, 1.182478e-08, 1.297227e-12,
         2.699643e-15, 1.755510e-15, 9.997385e-16, 8.598584e-29, 7.753024e-30,
         1.836567e-30)
pr_result <- conj_grad(par = c(-1.2, 1), fn = rosenbrock_banana$fr,
                      gr = rosenbrock_banana$grr, eps = .Machine$double.xmin,
                      reltol = .Machine$double.xmin,
                      abstol = .Machine$double.xmin, max_line_fn = 20)
expect_equal(pr_result$value, 4.979684e-30, tolerance = 1e-7)
expect_equal(pr_result$values, fxs, tolerance = 1e-7)
expect_equal(pr_result$iter, 26)
expect_equal(pr_result$counts, c(105, 105))

# example from http://learning.eng.cam.ac.uk/carl/code/minimize/
fxs00 <- c(1.00000000000000, 0.77160942667725, 0.58224024884105,
           0.40492742502160, 0.32466327341368, 0.28960411147824,
           0.07623420070067, 0.06786211944378, 0.03378423679313,
           0.00108990808914, 0.00108795243321, 0.00008974308332,
           0.00000012183819, 0.00000000675602, 0.00000000000000,
           0.00000000000000, 0.00000000000000, 0.00000000000000,
           0.00000000000000)
pr_result00 <- conj_grad(par = c(0, 0), fn = rosenbrock_banana$fr,
                       gr = rosenbrock_banana$grr, eps = .Machine$double.xmin,
                       reltol = .Machine$double.xmin,
                       abstol = .Machine$double.xmin, max_line_fn = 20)
expect_equal(pr_result00$value, 4.930381e-32, tolerance = 1e-7)
expect_equal(pr_result00$values, fxs00, tolerance = 1e-7)
expect_equal(pr_result00$iter, 19)
expect_equal(pr_result00$counts, c(79, 79))

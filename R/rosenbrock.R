## Rosenbrock Banana function
rosenbrock_banana <- list(
  fr = function(x) {
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1) ^ 2 + (1 - x1) ^ 2
  },
  grr = function(x) { ## Gradient of 'fr'
    x1 <- x[1]
    x2 <- x[2]
    c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
      200 *      (x2 - x1 * x1))
  }
)

fcn <- function(x) {
  list(f = rosenbrock_banana$fr(x), g = rosenbrock_banana$grr(x))
}

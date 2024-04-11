library(rgl)
library(mvtnorm)

############## benchmark functions
rastrigin_func <- function(x, A = 10) {
    n <- length(x)
    return(A * n + sum(x^2) - A * sum(cos(2 * pi * x)))
}

gaussian_pdf_func <- function(x) {
    f1 <- dmvnorm(x, mean = c(0.25, 0.75), sigma = diag(rep(0.01, d)))
    f2 <- dmvnorm(x, mean = c(0.75, 0.25), sigma = diag(rep(0.01, d)))

    return(f1 - f2)
}

ackley_func <- function(x) {
    return(-20 * exp(-0.2 * sqrt(0.5 * sum(x^2))) - exp(0.5 * sum(cos(2 * pi * x))) + exp(1) + 20)
}

sphere_func <- function(x) {
    return(sum(x^2))
}

############# grid points
s <- 100
grid_points <- matrix(rep(0, d * s^d), nrow = s^d)
for (c in 1:d) {
    grid_points[, c] <- rep(kronecker(seq(0, 1, length.out = s), rep(1, s^(d - c))), s^(c - 1))
}
# plot(grid_points)


zvalue <- apply(grid_points, 1, rastrigin_func)
# zvalue <- apply(grid_points, 1, gaussian_pdf_func)
# zvalue <- apply(grid_points, 1, ackley_func)
# zvalue <- apply(grid_points, 1, sphere_func)

open3d()
plot3d(grid_points[, 1], grid_points[, 2], zvalue)

library("hetGP")

f <- function(x) {
    return(f1d2(sum(x)))
}

y <- function(x) {
    # set.seed(123)
    return(f(x) + rnorm(1, sd = 0.01))
}

predict_mean <- function(X0, XN, Y, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    x0 <- matrix(X0, nrow = 1, byrow = TRUE)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type)

    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1]))

    L <- t(chol(K))
    v <- solve(L, t(KX))
    u <- solve(L, Y)
    return(t(v) %*% u)
}

predict_mean2 <- function(X0, XN, Y, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    x0 <- matrix(X0, nrow = 1, byrow = TRUE)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type)

    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1]))

    return(KX %*% solve(K) %*% Y)
}


N <- 20
cov_type <- "Gaussian"
nu.hat <- 1
theta <- 0.1
noise <- 0.01

XN <- matrix(runif(N, 0, 1), nrow = N)
Y <- apply(XN, 1, y)

X0 <- runif(1, 0, 1)
p1 <- predict_mean(X0, XN, Y, cov_type, nu.hat, theta, noise)
p2 <- predict_mean2(X0, XN, Y, cov_type, nu.hat, theta, noise)

print(c(p1, p2))

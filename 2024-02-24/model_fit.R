library("hetGP")
library("Metrics")
library("etasFLP")
library("psych")

setwd("D:/workspace/projects/R/hetGP")

# real function
f <- function(x) {
    return(f1d2(sum(x)) + sin(sum(x)^2))
}

# real function with constant noise = 0.01
y <- function(x) {
    # set.seed(123)
    return(f(x) + rnorm(1, sd = 0.01))
}

likelihood <- function(par, XN, Y, cov_type = "Gaussian") {
    theta <- par[1:ncol(XN)]
    sigma2 <- par[ncol(XN) + 1]
    n <- length(Y)
    K <- cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(sigma2, n)

    L <- t(chol(K))
    v <- solve(L, Y)

    ldetK <- determinant(K, logarithm = TRUE)$modulus

    l <- (n / 2) * log(t(v) %*% v) + (1 / 2) * ldetK
    return(l)
}

likelihood_derivative_validate <- function(par, XN, Y, cov_type = "Gaussian") {
    n <- length(Y)
    theta <- par[1:ncol(XN)]
    sigma2 <- par[ncol(XN) + 1]
    K <- cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(sigma2, n)

    Ki <- solve(K)
    KiY <- Ki %*% Y

    dlltheta <- rep(0, length(theta))
    for (k in 1:length(dlltheta)) {
        dotK <- K * as.matrix(dist(XN[, k]))^2 / (theta[k]^2)
        dlltheta[k] <- n * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) -
            sum(diag(Ki %*% dotK))
    }
    dllsigma2 <- n * t(KiY) %*% KiY / (t(Y) %*% KiY) - sum(diag(Ki))
    g <- -c(dlltheta / 2, dllsigma2 / 2)
    return(g)
}

likelihood_derivative <- function(par, XN, Y, cov_type = "Gaussian") {
    n <- length(Y)
    theta <- par[1:ncol(XN)]
    sigma2 <- par[ncol(XN) + 1]
    K <- cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(sigma2, n)

    Ki <- solve(K)
    KiY <- Ki %*% Y

    dlltheta <- rep(0, length(theta))
    for (k in 1:length(dlltheta)) {
        dotK <- K * as.matrix(dist(XN[, k]))^2 / (theta[k]^2)
        dlltheta[k] <- n * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) -
            tr(Ki %*% dotK)
    }
    dllsigma2 <- n * t(KiY) %*% KiY / (t(Y) %*% KiY) - tr(Ki)
    g <- -c(dlltheta / 2, dllsigma2 / 2)
    return(g)
}

mle_fit <- function(XN, Y, cov_type = "Gaussian") {
    fit <- optim(
        c(rep(0.1), 0.01),
        likelihood, likelihood_derivative,
        method = "L-BFGS-B",
        lower = c(0.0001, 0.0001), upper = c(5, 0.1),
        XN = XN, Y = Y, cov_type = cov_type
    )

    theta <- fit$par[1:ncol(XN)]
    sigma2 <- fit$par[ncol(XN) + 1]

    K <- cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(sigma2, length(Y))
    L <- t(chol(K))
    v <- solve(L, Y)

    nu_hat <- t(v) %*% v / length(Y)

    return(list(theta = theta, sigma2 = sigma2, nu_hat = nu_hat))
}

# nLL <- function(par, X, Y) {
#     theta <- par[1:ncol(X)]
#     tau2 <- par[ncol(X) + 1]
#     n <- length(Y)
#     K <- cov_gen(X1 = X, theta = theta) + diag(tau2, n)
#     Ki <- solve(K)
#     ldetK <- determinant(K, logarithm = TRUE)$modulus
#     (n / 2) * log(t(Y) %*% Ki %*% Y) + (1 / 2) * ldetK
# }

# gnLL <- function(par, X, Y) {
#     n <- length(Y)
#     theta <- par[1:ncol(X)]
#     tau2 <- par[ncol(X) + 1]
#     K <- cov_gen(X1 = X, theta = theta) + diag(tau2, n)
#     Ki <- solve(K)
#     KiY <- Ki %*% Y
#     dlltheta <- rep(0, length(theta))
#     for (k in 1:length(dlltheta)) {
#         dotK <- K * as.matrix(dist(X[, k]))^2 / (theta[k]^2)
#         dlltheta[k] <- n * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) -
#             sum(diag(Ki %*% dotK))
#     }
#     dlltau2 <- n * t(KiY) %*% KiY / (t(Y) %*% KiY) - sum(diag(Ki))
#     -c(dlltheta / 2, dlltau2 / 2)
# }

N <- 5
XN <- matrix(seq(0, 1, length.out = N), nrow = N)
# x1 <- kronecker(rep(1, N), seq(0, 1, length = N))
# x2 <- kronecker(seq(0, 1, length = N), rep(1, N))
# XN <- matrix(c(x1, x2), nrow = N^2)
print(XN)
Y <- apply(XN, 1, y)
print(Y)

# library("lhs")
# XN <- 6 * randomLHS(40, 2) - 2
# XN <- rbind(XN, XN)
# Y <- XN[, 1] * exp(-XN[, 1]^2 - XN[, 2]^2) + rnorm(nrow(XN), sd = 0.01)
# print(XN)
# print(Y)

delta <- 1e-9
l1 <- likelihood(c(0.5, 0.5 + delta), XN, Y)
l2 <- likelihood(c(0.5, 0.5), XN, Y)
dl <- likelihood_derivative(c(0.5, 0.5), XN, Y)
print(c((l1 - l2) / delta, dl))

# params <- mle_fit(XN, Y)
# print(params)

# mod <- mleHomGP(X = XN, Z = Y, lower = c(0.0001, 0.0001), upper = c(5, 0.1), known = list(beta0 = 0), init = c(list(theta = c(0.1), g = 0.01)))


# out <- optim(c(rep(0.1, 2), 0.1 * var(Y)), nLL, gnLL,
#     method = "L-BFGS-B",
#     lower = Lwr, upper = c(rep(Upr, 2), var(Y)), X = XN, Y = Y
# )
# out$par

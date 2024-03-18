library("optimx")
library("hetGP")
library("Metrics") # compute rmse
library("psych") # compute trace
library("BB")

setwd("D:/workspace/projects/R/hetGP")

f <- function(x) {
    return(sum(cos(8 * pi * x)))
}

# real function with constant noise = 0.01
y <- function(x) {
    return(f(x) + rnorm(1, sd = 0.01))
}

h <- function(X0) {
    # return(1 / sqrt(2 * pi) * exp(-X0^2 / 2))
    return(1)
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
    return(l * 10)
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
    return(g * 10)
}

mle_fit <- function(XN, Y, cov_type = "Gaussian") {
    fit <- optim(
        c(0.1, 0.01),
        likelihood, likelihood_derivative,
        method = "L-BFGS-B",
        lower = c(1e-4, 1e-4), upper = c(1, 0.1),
        XN = XN, Y = Y, cov_type = cov_type
    )

    theta <- fit$par[1:ncol(XN)]
    sigma2 <- fit$par[ncol(XN) + 1]

    K <- cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(sigma2, length(Y))
    L <- t(chol(K))
    v <- solve(L, Y)

    nu_hat <- as.numeric(t(v) %*% v / length(Y))
    return(list(theta = theta, sigma2 = sigma2, nu_hat = nu_hat))
}

mse <- function(par, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    x0 <- matrix(par, nrow = 1, byrow = TRUE)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type)

    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1]))

    L <- t(chol(K))
    v <- solve(L, t(KX))

    mse <- nu.hat - t(v) %*% v
    return(mse)
}

diff <- function(par, XN, Xnew, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    x0 <- matrix(par, nrow = 1, byrow = TRUE)
    xN1 <- matrix(Xnew, nrow = 1, byrow = TRUE)

    sigma2 <- mse(Xnew, XN, cov_type, nu.hat, theta, noise)
    sigma2 <- sigma2 + noise

    K1 <- nu.hat * cov_gen(X1 = x0, X2 = xN1, theta = theta, type = cov_type)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type)
    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1]))
    Kn <- nu.hat * cov_gen(X1 = xN1, X2 = XN, theta = theta, type = cov_type)

    L <- t(chol(K))
    v <- solve(L, t(KX))
    u <- solve(L, t(Kn))

    # print(c(K1, t(v) %*% u, K1 - t(v) %*% u))

    return(-1 / sigma2 * (K1 - t(v) %*% u)^2)
    # return((K1 - t(v) %*% u)^2)
}

g2 <- function(par, XN, XN1, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    x0 <- matrix(par, nrow = 1, byrow = TRUE)
    xN1 <- matrix(XN1, nrow = 1, byrow = TRUE)

    K1 <- nu.hat * cov_gen(X1 = x0, X2 = xN1, theta = theta, type = cov_type)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type)
    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1]))
    Kn <- nu.hat * cov_gen(X1 = xN1, X2 = XN, theta = theta, type = cov_type)

    L <- t(chol(K))
    v <- solve(L, t(KX))
    u <- solve(L, t(Kn))

    return(-(K1 - t(v) %*% u)^2)
}

# Laplace_method_v0 <- function(par, X0, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
#     x0 <- matrix(X0, nrow = 1, byrow = TRUE) # Ex
#     xN1 <- matrix(par, nrow = 1, byrow = TRUE) # xn+1
#     K1 <- nu.hat * cov_gen(X1 = x0, X2 = xN1, theta = theta, type = cov_type) # k(x, xn+1)
#     KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type) # kn(x)
#     K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1])) # Kn
#     Kn <- nu.hat * cov_gen(X1 = xN1, X2 = XN, theta = theta, type = cov_type) # kn(xn+1)
#     L <- t(chol(K)) # Kn = LL^T
#     v1 <- solve(L, t(KX))
#     v2 <- solve(L, diag(X0 - c(XN)) %*% t(KX))
#     v3 <- solve(L, diag((X0 - c(XN))^2) %*% t(KX))
#     u <- solve(L, t(Kn))
#     sigma2 <- nu.hat - t(u) %*% u + noise
#     g <- K1 - t(v1) %*% u
#     gd1 <- -2 / theta * (x0 - xN1) * K1 + 2 / theta * t(v2) %*% u
#     gd2 <- (-2 / theta + 4 / theta^2 * (x0 - xN1)^2) * K1 + 2 / theta * t(v1) %*% u - 4 / theta^2 * t(v3) %*% u
#     I <- -1 / sigma2 * sqrt(pi) * h(X0) * abs(g)^3 / sqrt(abs(gd2 * g - gd1^2))
#     return(I)
# }

Laplace_method_v1 <- function(par, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    # out <- multistart(matrix(seq(0, 1, length = 50), nrow = 50),
    #     fn = g2,
    #     method = "L-BFGS-B",
    #     lower = c(0), upper = c(1),
    #     XN = XN, XN1 = par, cov_type = cov_type, nu.hat = nu.hat, theta = theta, noise = noise
    # )
    # X0 <- out$p1[which.min(out$value)]
    X0 <- par

    x0 <- matrix(X0, nrow = 1, byrow = TRUE) # argmin -g^2
    xN1 <- matrix(par, nrow = 1, byrow = TRUE) # xn+1

    K1 <- nu.hat * cov_gen(X1 = x0, X2 = xN1, theta = theta, type = cov_type) # k(x, xn+1)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type) # kn(x)
    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1])) # Kn
    Kn <- nu.hat * cov_gen(X1 = xN1, X2 = XN, theta = theta, type = cov_type) # kn(xn+1)

    L <- t(chol(K)) # Kn = LL^T

    v1 <- solve(L, t(KX))
    v2 <- solve(L, diag(X0 - c(XN)) %*% t(KX))
    v3 <- solve(L, diag((X0 - c(XN))^2) %*% t(KX))
    u <- solve(L, t(Kn))

    sigma2 <- nu.hat - t(u) %*% u + noise

    g <- K1 - t(v1) %*% u
    gd1 <- -2 / theta * (x0 - xN1) * K1 + 2 / theta * t(v2) %*% u
    gd2 <- (-2 / theta + 4 / theta^2 * (x0 - xN1)^2) * K1 + 2 / theta * t(v1) %*% u - 4 / theta^2 * t(v3) %*% u

    I <- -1 / sigma2 * sqrt(pi) * h(X0) * abs(g)^3 / sqrt(abs(gd2 * g - gd1^2))
    return(I)
}

Laplace_method_v2 <- function(par, X0, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    x0 <- matrix(X0, nrow = 1, byrow = TRUE) # argmax h
    xN1 <- matrix(par, nrow = 1, byrow = TRUE) # xn+1

    K1 <- nu.hat * cov_gen(X1 = x0, X2 = xN1, theta = theta, type = cov_type) # k(x, xn+1)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type) # kn(x)
    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1])) # Kn
    Kn <- nu.hat * cov_gen(X1 = xN1, X2 = XN, theta = theta, type = cov_type) # kn(xn+1)

    L <- t(chol(K)) # Kn = LL^T

    v1 <- solve(L, t(KX))
    v2 <- solve(L, diag(X0 - c(XN)) %*% t(KX))
    v3 <- solve(L, diag((X0 - c(XN))^2) %*% t(KX))
    u <- solve(L, t(Kn))

    sigma2 <- nu.hat - t(u) %*% u + noise

    g <- K1 - t(v1) %*% u
    gd1 <- -2 / theta * (x0 - xN1) * K1 + 2 / theta * t(v2) %*% u
    gd2 <- (-2 / theta + 4 / theta^2 * (x0 - xN1)^2) * K1 + 2 / theta * t(v1) %*% u - 4 / theta^2 * t(v3) %*% u

    I <- -1 / sigma2 * sqrt(pi) * h(X0) * abs(g)^3 / sqrt(abs(gd2 * g - gd1^2))
    return(I)
}

MC_integral <- function(XN1, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    xgrid <- sort(runif(3000))
    # delta <- xgrid[2] - xgrid[1]
    I <- 0
    for (i in 1:length(xgrid)) {
        I <- I + diff(xgrid[i], XN, XN1, cov_type, nu.hat, theta, noise)
    }
    return(I / length(xgrid))
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

cov_type <- "Gaussian"
N <- 30
XN <- matrix(seq(0, 1, length.out = N), nrow = N)
XN1 <- runif(1)
print(sprintf("XN1: %f", XN1))
Y <- apply(XN, 1, y)

mod <- mleHomGP(X = XN, Z = Y, covtype = cov_type)
theta <- mod$theta
nu.hat <- mod$nu_hat
noise <- mod$g
print(sprintf("theta: %f, nu.hat: %f, noise: %f", theta, nu.hat, noise))

#################### validation
# out <- multistart(matrix(seq(0, 1, length = 50), nrow = 50),
#     fn = g2,
#     method = "L-BFGS-B",
#     lower = c(0), upper = c(1),
#     XN = XN, XN1 = XN1, cov_type = cov_type, nu.hat = nu.hat, theta = theta, noise = noise
# )
# X0 <- out$p1[which.min(out$value)]
# print(sprintf("X0: %f", X0))
# xgrid <- seq(0, 1, length = 1000)
# ygrid <- c()
# for (i in 1:1000) {
#     ygrid <- c(ygrid, g2(xgrid[i], XN, XN1, cov_type, nu.hat, theta, noise))
# }
# plot(xgrid, ygrid, col = "black") +
#     abline(v = XN1, col = "blue") +
#     abline(v = X0, col = "red")
# print(c(
#     Laplace_method_v1(XN1, XN, cov_type, nu.hat, theta, noise),
#     MC_integral(XN1, XN, cov_type, nu.hat, theta, noise)
# ))

#################### plot of Laplace and Monte Carlo
# xgrid <- seq(0, 1, length = 500)
# ygrid <- c()
# zgrid <- c()
# for (i in 1:500) {
#     XN1 <- xgrid[i]
#     ygrid <- c(ygrid, Laplace_method_v1(XN1, XN, cov_type, nu.hat, theta, noise))
#     zgrid <- c(zgrid, MC_integral(XN1, XN, cov_type, nu.hat, theta, noise))
# }
# plot(xgrid, ygrid, col = "black") +
#     lines(xgrid, zgrid, col = "blue")


# #################### train
step <- 50
start <- Sys.time()
for (k in 1:step) {
    out <- multistart(matrix(seq(0, 1, length = 10), nrow = 10),
        fn = Laplace_method_v1,
        method = "L-BFGS-B",
        lower = c(0), upper = c(1),
        XN = XN, cov_type = cov_type, nu.hat = nu.hat, theta = theta, noise = noise
    )
    # print(out)
    opt <- out$p1[which.min(out$value)]

    Xnew <- matrix(rep(opt, 1), ncol = 1)
    Ynew <- apply(Xnew, 1, y)

    XN <- rbind(XN, Xnew)
    Y <- c(Y, Ynew)

    # fit <- mle_fit(XN, Y, cov_type)
    # theta <- as.numeric(fit$theta)
    # nu.hat <- fit$nu_hat
    # noise <- fit$sigma2
    mod <- update(mod, Xnew = Xnew, Znew = Ynew, ginit = mod$g * 1.01)
    theta <- mod$theta
    nu.hat <- mod$nu_hat
    noise <- mod$g

    print(sprintf("iteration: %d, Xnew: %f, theta: %f, nu_hat: %f, noise: %f", k, Xnew, theta, nu.hat, noise))
}
print(sprintf("train time: %f", Sys.time() - start))


# #################### visualization
# par(mfrow = c(1, 1))

xgrid <- seq(0, 1, length.out = 1000)
pgrid <- c()
fgrid <- c()
for (i in 1:length(xgrid)) {
    fgrid <- c(fgrid, f(xgrid[i]))
    pgrid <- c(pgrid, predict_mean(xgrid[i], XN, Y, cov_type, nu.hat, theta, noise))
}
plot(xgrid, fgrid, col = "black") +
    lines(xgrid, pgrid, col = "blue") +
    abline(v = XN[1:N], col = "yellow") +
    abline(v = XN[(N + 1):(N + step)], col = "green", lty = 3)

print(rmse(fgrid, pgrid))


# ##################### train IMSPE
# XN <- matrix(XN[1:N, ], nrow = N)
# Y <- Y[1:N]

# start <- Sys.time()
# for (i in 1:step) {
#     out <- IMSPE_optim(mod)
#     Xnew <- matrix(out$par, nrow = 1)
#     Ynew <- apply(Xnew, 1, y)
#     # print(Xnew)

#     XN <- rbind(XN, Xnew)
#     Y <- c(Y, Ynew)

#     mod <- update(mod, Xnew = Xnew, Znew = Ynew, ginit = mod$g * 1.01)
#     print(sprintf("iteration: %d, Xnew: %f, theta: %f, nu_hat: %f, noise: %f", i, Xnew, mod$theta, mod$nu_hat, mod$g))
# }
# print(sprintf("train time: %f", Sys.time() - start))

# pgrid <- c()
# for (i in 1:length(xgrid)) {
#     pgrid <- c(pgrid, predict_mean(xgrid[i], XN, Y, cov_type, nu.hat, theta, noise))
# }
# plot(xgrid, fgrid, col = "black") +
#     lines(xgrid, pgrid, col = "blue") +
#     abline(v = XN[1:N], col = "yellow") +
#     abline(v = XN[(N + 1):(N + step)], col = "green", lty = 3)
# print(rmse(fgrid, pgrid))

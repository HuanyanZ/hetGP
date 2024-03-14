# set.seed(123)
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
    # set.seed(123)
    return(f(x) + rnorm(1, sd = 0.01))
}

# MSE function
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

integral <- function(par, X0, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01, var = 1 / 12) {
    x0 <- matrix(X0, nrow = 1, byrow = TRUE) # Ex
    xN1 <- matrix(par, nrow = 1, byrow = TRUE) # xn+1

    K1 <- nu.hat * cov_gen(X1 = x0, X2 = xN1, theta = theta, type = cov_type) # k(x, xn+1)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type) # kn(x)
    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1])) # Kn
    Kn <- nu.hat * cov_gen(X1 = xN1, X2 = XN, theta = theta, type = cov_type) # kn(xn+1)

    L <- t(chol(K))

    v1 <- solve(L, t(KX))
    v2 <- solve(L, diag(X0 - c(XN)) %*% t(KX))
    v3 <- solve(L, diag((X0 - c(XN))^2) %*% t(KX))
    u <- solve(L, t(Kn))

    sigma2 <- nu.hat - t(u) %*% u + noise

    g <- K1 - t(v1) %*% u
    gd1 <- -2 / theta * (x0 - xN1) * K1 + 2 / theta * t(v2) %*% u
    gd2 <- (-2 / theta + 4 / theta^2 * (x0 - xN1)^2) * K1 + 2 / theta * t(v1) %*% u - 4 / theta^2 * t(v3) %*% u

    # print(c(g, gd1, gd2, g^2 + var * (gd1^2 + g * gd2)))
    I <- -1 / sigma2 * (g^2 + var * (gd1^2 + g * gd2))
    return(I)
}

integral_exp <- function(par, X0, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01, var = 1 / 12) {
    return(exp(integral(par, X0, XN, cov_type, nu.hat, theta, noise, var)))
}

MC_integral <- function(Xnew, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    xgrid <- seq(0, 1, length = 1000)
    delta <- xgrid[2] - xgrid[1]
    I <- 0
    for (i in 1:length(xgrid)) {
        I <- I + delta * diff(xgrid[i], XN, Xnew, cov_type, nu.hat, theta, noise)
    }
    return(I / length(xgrid))
}

IMSE <- function(XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    xgrid <- seq(0, 1, length = 1000)
    delta <- xgrid[2] - xgrid[1]
    I <- 0
    for (i in 1:length(xgrid)) {
        I <- I + delta * mse(xgrid[i], XN, cov_type, nu.hat, theta, noise)
    }
    return(I / length(xgrid))
}

cov_type <- "Gaussian"
N <- 10
XN <- matrix(seq(0, 1, length.out = N), nrow = N)
# XN <- matrix(sort(runif(N, 0, 1)), nrow = N)

# Xnew <- matrix(0.5)
# XN1 <- rbind(XN, Xnew)

# IN <- IMSE(XN)
# IN1 <- IMSE(XN1)
# print(c(IN1 - IN, MC_integral(Xnew, XN)))

Y <- apply(XN, 1, y)

mod <- mleHomGP(X = XN, Z = Y, covtype = cov_type)

step <- 6
par(mfrow = c(6, 1))

xgrid <- seq(0, 1, length = 500)


##################### train Taylor expansion IMSE
for (k in 1:step) {
    ygrid <- c()
    for (i in 1:500) {
        ygrid <- c(ygrid, integral(xgrid[i], 0.5, XN))
    }

    xnew <- xgrid[which.min(ygrid)]
    plot(xgrid, ygrid, col = "black") +
        abline(v = xnew, col = "red")
    Ynew <- apply(matrix(xnew), 1, y)

    XN <- rbind(XN, matrix(xnew))
    mod <- update(mod, Xnew = matrix(xnew), Znew = Ynew, ginit = mod$g * 1.01)
}






##################### train Monte Carlo IMSE
# XN <- matrix(XN[1:N], nrow = N)
# for (k in 1:step) {
#     ygrid <- c()
#     for (i in 1:500) {
#         ygrid <- c(ygrid, MC_integral(xgrid[i], XN))
#     }

#     xnew <- xgrid[which.min(ygrid)]
#     plot(xgrid, ygrid, col = "black") +
#         abline(v = xnew, col = "red")
#     Ynew <- apply(matrix(xnew), 1, y)

#     XN <- rbind(XN, matrix(xnew))
#     mod <- update(mod, Xnew = matrix(xnew), Znew = Ynew, ginit = mod$g * 1.01)
# }


##################### train IMSPE
# XN <- matrix(XN[1:N], nrow = N)
# for (k in 1:step) {
#     ygrid <- c()
#     for (i in 1:500) {
#         ygrid <- c(ygrid, crit_IMSPE(xgrid[i], mod))
#     }

#     xnew <- IMSPE_optim(mod, h = 0)$par
#     plot(xgrid, ygrid, col = "green") +
#         abline(v = xnew, col = "red")
#     Ynew <- apply(matrix(xnew), 1, y)

#     XN <- rbind(XN, matrix(xnew))
#     mod <- update(mod, Xnew = matrix(xnew), Znew = Ynew, ginit = mod$g * 1.01)
# }

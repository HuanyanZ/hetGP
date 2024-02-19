library("hetGP")
library("DEoptim")
library("Metrics")
library("optimx")

setwd("D:/workspace/projects/R/hetGP")

# real function
f <- function(x) {
    return(f1d2(sum(x)) + exp(-sum(x)) + sin(sum(x)^2))
}

# real function with constant noise = 0.01
y <- function(x) {
    # set.seed(123)
    return(f(x) + rnorm(1, sd = 0.01))
}

g <- function(par, X0, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    x0 <- matrix(X0, nrow = 1, byrow = TRUE) # Ex
    xN1 <- matrix(par, nrow = 1, byrow = TRUE) # xn+1

    K1 <- nu.hat * cov_gen(X1 = x0, X2 = xN1, theta = theta, type = cov_type) # k(x, xn+1)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type) # kn(x)
    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1])) # Kn
    Kn <- nu.hat * cov_gen(X1 = xN1, X2 = XN, theta = theta, type = cov_type) # kn(xn+1)

    L <- t(chol(K))
    v <- solve(L, t(KX))
    u <- solve(L, t(Kn))

    g <- K1 - t(v) %*% u
    return(g)
}

g_derivative_1 <- function(par, X0, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    x0 <- matrix(X0, nrow = 1, byrow = TRUE) # Ex
    xN1 <- matrix(par, nrow = 1, byrow = TRUE) # xn+1

    K1 <- nu.hat * cov_gen(X1 = x0, X2 = xN1, theta = theta, type = cov_type) # k(x, xn+1)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type) # kn(x)
    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1])) # Kn
    Kn <- nu.hat * cov_gen(X1 = xN1, X2 = XN, theta = theta, type = cov_type) # kn(xn+1)

    L <- t(chol(K))
    v <- solve(L, diag(X0 - c(XN)) %*% t(KX))
    u <- solve(L, t(Kn))

    gd1 <- -2 / theta * (x0 - xN1) * K1 + 2 / theta * t(v) %*% u
    return(gd1)
}

g_derivative_2 <- function(par, X0, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    x0 <- matrix(X0, nrow = 1, byrow = TRUE) # Ex
    xN1 <- matrix(par, nrow = 1, byrow = TRUE) # xn+1

    K1 <- nu.hat * cov_gen(X1 = x0, X2 = xN1, theta = theta, type = cov_type) # k(x, xn+1)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type) # kn(x)
    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1])) # Kn
    Kn <- nu.hat * cov_gen(X1 = xN1, X2 = XN, theta = theta, type = cov_type) # kn(xn+1)

    L <- t(chol(K))
    v1 <- solve(L, t(KX))
    v2 <- solve(L, diag((X0 - c(XN))^2) %*% t(KX))
    u <- solve(L, t(Kn))

    gd2 <- (-2 / theta + 4 / theta^2 * (x0 - xN1)^2) * K1 + 2 / theta * t(v1) %*% u - 4 / theta^2 * t(v2) %*% u
    return(gd2)
}

# g_validate <- function(par, X0, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
#     K <- c()
#     for (i in 1:length(XN)) {
#         for (j in 1:length(XN)) {
#             K <- c(K, -(XN[i, ] - XN[j, ])^2 / theta)
#         }
#     }
#     K <- matrix(nu.hat * exp(K), nrow = length(XN)) + diag(rep(noise, length(XN)))

#     Kn <- c()
#     for (i in 1:length(XN)) {
#         Kn <- c(Kn, -(par - XN[i, ])^2 / theta)
#     }
#     Kn <- matrix(nu.hat * exp(Kn), ncol = 1)

#     KX <- c()
#     for (i in 1:length(XN)) {
#         KX <- c(KX, -(X0 - XN[i, ])^2 / theta)
#     }
#     KX <- matrix(nu.hat * exp(KX), ncol = 1)

#     K1 <- nu.hat * exp(-(X0 - par)^2 / theta)

#     return(K1 - t(KX) %*% solve(K) %*% Kn)
# }

# g_derivative_1_validate <- function(par, X0, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
#     x0 <- matrix(X0, nrow = 1, byrow = TRUE) # Ex
#     xN1 <- matrix(par, nrow = 1, byrow = TRUE) # xn+1

#     K1 <- nu.hat * cov_gen(X1 = x0, X2 = xN1, theta = theta, type = cov_type) # k(x, xn+1)
#     KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type) # kn(x)
#     K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1])) # Kn
#     Kn <- nu.hat * cov_gen(X1 = xN1, X2 = XN, theta = theta, type = cov_type) # kn(xn+1)

#     return(-2 / theta * (x0 - xN1) * K1 + 2 / theta * Kn %*% solve(K) %*% diag(X0 - c(XN)) %*% t(KX))
# }


# g_derivative_2_validate <- function(par, X0, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
#     x0 <- matrix(X0, nrow = 1, byrow = TRUE) # Ex
#     xN1 <- matrix(par, nrow = 1, byrow = TRUE) # xn+1

#     K1 <- nu.hat * cov_gen(X1 = x0, X2 = xN1, theta = theta, type = cov_type) # k(x, xn+1)
#     KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type) # kn(x)
#     K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1])) # Kn
#     Kn <- nu.hat * cov_gen(X1 = xN1, X2 = XN, theta = theta, type = cov_type) # kn(xn+1)

#     return(-2 / theta * K1 + 4 / (theta^2) * (x0 - xN1)^2 * K1 + 2 / theta * KX %*% solve(K) %*% t(Kn) - 4 / (theta^2) * Kn %*% solve(K) %*% diag(X0 - c(XN)) %*% diag(X0 - c(XN)) %*% t(KX))
# }

# mse <- function(par, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
#     x0 <- matrix(par, nrow = 1, byrow = TRUE)
#     KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type)

#     K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1]))

#     L <- t(chol(K))
#     v <- solve(L, t(KX))

#     mse <- nu.hat - t(v) %*% v
#     return(mse)
# }

# integral_validate <- function(par, X0, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01, var = 1 / 12) {
#     sigma2 <- mse(par, XN, cov_type, nu.hat, theta, noise)
#     sigma2 <- sigma2 + noise

#     g <- g(par, X0, XN, cov_type, nu.hat, theta, noise)
#     gd1 <- g_derivative_1(par, X0, XN, cov_type, nu.hat, theta, noise)
#     gd2 <- g_derivative_2(par, X0, XN, cov_type, nu.hat, theta, noise)

#     return(-1 / sigma2 * (g^2 + var * (gd1^2 + g * gd2)))
# }

# integral_MC <- function(par, xgrid, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01, var = 1 / 12) {
#     sigma2 <- mse(par, XN, cov_type, nu.hat, theta, noise)
#     sigma2 <- sigma2 + noise

#     I <- 0
#     for (i in 1:length(xgrid)) {
#         g <- g(par, xgrid[i], XN, cov_type, nu.hat, theta, noise)
#         I <- I + g^2
#     }
#     I <- -1 / sigma2 * I / length(xgrid)

#     return(I)
# }

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
    return(-1 / sigma2 * (g^2 + var * (gd1^2 + g * gd2)))
}

integral_derivative <- function(par, X0, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01, var = 1 / 12) {
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

    u1 <- solve(L, t(Kn))
    u2 <- solve(L, diag(par - c(XN)) %*% t(Kn))

    sigma2 <- nu.hat - t(u1) %*% u1 + noise

    g <- K1 - t(v1) %*% u1
    g1 <- -2 / theta * (x0 - xN1) * K1 + 2 / theta * t(v2) %*% u1
    g2 <- (-2 / theta + 4 / theta^2 * (x0 - xN1)^2) * K1 + 2 / theta * t(v1) %*% u1 - 4 / theta^2 * t(v3) %*% u1

    E <- g^2 + var * (g1^2 + g * g2)

    gd <- -2 / theta * (par - X0) * K1 + 2 / theta * t(v1) %*% u2
    g1d <- 2 / theta * K1 - 4 / theta^2 * (par - X0)^2 * K1 - 4 / theta^2 * t(v2) %*% u2
    g2d <- 12 / theta^2 * (par - X0) * K1 - 8 / theta^3 * (par - X0)^3 * K1 - 4 / theta^2 * t(v1) %*% u2 + 8 / theta^3 * t(v3) %*% u2

    Ed <- 2 * g * gd + var * (2 * g1 * g1d + gd * g2 + g * g2d)

    Id <- -1 * (2 / sigma2^2 * E * (-2 / theta) * t(u1) %*% u2 + 1 / sigma2 * Ed)
    return(Id)
}

# predict mean of posterior distribution
predict_mean <- function(X0, XN, Y, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    x0 <- matrix(X0, nrow = 1, byrow = TRUE)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type)

    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1]))

    L <- t(chol(K))
    v <- solve(L, t(KX))
    u <- solve(L, Y)
    return(t(v) %*% u)
}

N <- 10
cov_type <- "Gaussian"
noise <- 0.01

# xgrid <- sort(runif(1000, 0, 1))
xgrid <- sort(rnorm(1000, 0.85, 1.2))
# xgrid <- sort(rt(1000, 50))
# xgrid <- sort(rchisq(1000, 50))


xgrid <- (xgrid - min(xgrid)) / (max(xgrid) - min(xgrid))
X0 <- mean(xgrid)
var <- var(xgrid)
print(c(X0, var))

XN <- sample(xgrid, N)
XN <- matrix(XN, nrow = N)
# print(XN)

Y <- apply(XN, 1, y)
mod <- mleHetGP(X = XN, Z = Y, covtype = cov_type)
theta <- mod$theta
nu.hat <- mod$nu_hat
print(c(theta, nu.hat))

# Xnew <- sample(xgrid, 1)
# print(Xnew)

# g1 <- g(Xnew, 0.5, XN, cov_type = "Gaussian", nu.hat, theta, noise)
# g2 <- g_validate(Xnew, 0.5, XN, cov_type = "Gaussian", nu.hat, theta, noise)
# print(c(g1, g2))

# delta <- 1e-9
# g1 <- g(Xnew, 0.5, XN, cov_type = "Gaussian", nu.hat, theta, noise)
# g2 <- g(Xnew, 0.5 - delta, XN, cov_type = "Gaussian", nu.hat, theta, noise)
# gd1 <- g_derivative_1(Xnew, 0.5, XN, cov_type = "Gaussian", nu.hat, theta, noise)
# print(c(gd1, (g1 - g2) / delta))

# delta <- 1e-9
# g1 <- g_derivative_1(Xnew, 0.5, XN, cov_type = "Gaussian", nu.hat, theta, noise)
# g2 <- g_derivative_1(Xnew, 0.5 - delta, XN, cov_type = "Gaussian", nu.hat, theta, noise)
# gd2 <- g_derivative_2(Xnew, 0.5, XN, cov_type = "Gaussian", nu.hat, theta, noise)
# print(c(gd2, (g1 - g2) / delta))



# igrid <- c()
# for (i in 1:length(xgrid)) {
#     igrid <- c(igrid, integral(xgrid[i], X0, XN, cov_type = "Gaussian", nu.hat, theta, noise, var))
# }
# plot(xgrid, igrid)
# print(integral(Xnew, X0, XN, cov_type = "Gaussian", nu.hat, theta, noise, var))
# print(integral_MC(Xnew, xgrid, XN, cov_type = "Gaussian", nu.hat, theta, noise))


# delta <- 1e-9
# i1 <- integral(Xnew + delta, 0.5, XN, cov_type = "Gaussian", nu.hat, theta, noise)
# i2 <- integral(Xnew, 0.5, XN, cov_type = "Gaussian", nu.hat, theta, noise)
# id <- integral_derivative(Xnew, 0.5, XN, cov_type = "Gaussian", nu.hat, theta, noise)
# print(c(id, (i1 - i2) / delta))

# xgrid <- runif(2000, 0, 1)
# ygrid <- c()
# zgrid <- c()
# for (i in 1:2000) {
#     ygrid <- c(ygrid, integral(xgrid[i], 0.5, XN, cov_type = "Gaussian", nu.hat, theta, noise))
# }
# plt1 <- plot(xgrid, ygrid)



step <- 100
print(Sys.time())
for (i in 1:step) {
    # out <- DEoptim(
    #     integral,
    #     lower = c(0), upper = c(1),
    #     DEoptim.control(trace = FALSE, itermax = 100),
    #     X0 = X0, XN = XN, theta = theta, nu.hat = nu.hat, var = var,
    # )
    # opt <- out$optim$bestmem

    out <- optim(c(X0), integral, integral_derivative,
        method = "L-BFGS-B",
        lower = c(0), upper = c(1),
        X0 = X0, XN = XN, theta = theta, nu.hat = nu.hat, var = var,
    )
    opt <- out$par

    # out <- multistart(matrix(c(0, X0, 1), nrow = 3), integral, integral_derivative,
    #     method = "L-BFGS-B",
    #     lower = c(0), upper = c(1),
    #     X0 = X0, XN = XN, theta = theta, nu.hat = nu.hat, var = var,
    # )
    # opt <- out$p1[which.min(out$value)]

    XN <- rbind(XN, matrix(opt))
}
print(Sys.time())


# print(XN)


par(mfrow = c(2, 1))

Y <- apply(XN, 1, y)
pgrid <- c()
fgrid <- c()
for (i in 1:length(xgrid)) {
    fgrid <- c(fgrid, f(xgrid[i]))
    pgrid <- c(pgrid, predict_mean(xgrid[i], XN, Y, cov_type, nu.hat, theta, noise))
}
plt1 <- plot(xgrid, fgrid, col = "black") +
    lines(xgrid, pgrid, col = "blue") +
    abline(v = XN, col = "green") +
    abline(v = X0, col = "red")
# print(plt1)
print(rmse(fgrid, pgrid))



XN <- XN[1:N, ]
print(Sys.time())
for (i in 1:step) {
    out <- IMSPE_optim(mod)
    Xnew <- matrix(out$par, nrow = 1)
    Ynew <- apply(Xnew, 1, y)

    XN <- rbind(XN, Xnew)

    mod <- update(mod, Xnew = Xnew, Znew = Ynew, ginit = mod$g * 1.01)
}
print(Sys.time())

Y <- apply(XN, 1, y)

pgrid <- c()
for (i in 1:length(xgrid)) {
    pgrid <- c(pgrid, predict_mean(xgrid[i], XN, Y, cov_type, nu.hat, theta, noise))
}
plt2 <- plot(xgrid, fgrid, col = "black") +
    lines(xgrid, pgrid, col = "blue") +
    abline(v = XN, col = "green") +
    abline(v = X0, col = "red")
print(rmse(fgrid, pgrid))

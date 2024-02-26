# set.seed(123)
library("optimx")
library("hetGP")
library("Metrics") # compute rmse
library("psych") # compute trace
library("BB")

setwd("D:/workspace/projects/R/hetGP")

# real function
# f <- function(x) {
#     return(f1d2(sum(x)) + sin(sum(x)^2))
# }
f <- function(x) {
    return(sum(cos(4 * pi * x)))
}

# real function with constant noise = 0.01
y <- function(x) {
    # set.seed(123)
    return(f(x) + rnorm(1, sd = 0.01))
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
    return(-abs(I))
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

N <- 10
cov_type <- "Gaussian"
noise <- 0.01

xgrid <- sort(runif(1000, 0, 1))

X0 <- mean(xgrid)
var <- var(xgrid)
print(sprintf("X's mean: %f, X's variance: %f", X0, var))
# sprintf("Current working dir: %s", X0)

XN <- matrix(seq(0, 1, length.out = N), nrow = N)
# XN <- rbind(XN, XN)
print(sprintf("input XN: %f", XN))
Y <- apply(XN, 1, y)



# delta <- 1e-9
# i1 <- integral(X0 + delta, 0.5, XN, cov_type = "Gaussian", nu.hat, theta, noise)
# i2 <- integral(X0, 0.5, XN, cov_type = "Gaussian", nu.hat, theta, noise)
# id <- integral_derivative(X0, 0.5, XN, cov_type = "Gaussian", nu.hat, theta, noise)
# print(c(id, (i1 - i2) / delta))


mod <- mleHomGP(X = XN, Z = Y, covtype = cov_type)
theta <- mod$theta
nu.hat <- mod$nu_hat
print(sprintf("theta: %f, nu.hat: %f", theta, nu.hat))

par(mfrow = c(7, 1))

step <- 20
start <- Sys.time()
for (i in 1:step) {
    # out <- optim(c(X0), integral, integral_derivative,
    #     method = "L-BFGS-B",
    #     lower = c(0), upper = c(1),
    #     X0 = X0, XN = XN, nu.hat = nu.hat, theta = theta, noise = noise, var = var
    # )
    # opt <- out$par

    out <- multistart(matrix(seq(0, 1, length = 20), nrow = 20),
        fn = integral,
        # gr = integral_derivative,
        method = "L-BFGS-B",
        lower = c(0), upper = c(1),
        X0 = X0, XN = XN, nu.hat = nu.hat, theta = theta, noise = noise, var = var
    )
    # print(out)
    opt <- out$p1[which.min(out$value)]

    # out <- multiStart(matrix(seq(0, 1, length = 20), nrow = 20),
    #     fn = integral,
    #     action = "optimize",
    #     lower = c(0), upper = c(1),
    #     control = list(trace = FALSE),
    #     quiet = TRUE,
    #     X0 = X0, XN = XN, nu.hat = nu.hat, theta = theta, noise = noise, var = var
    # )
    # opt <- as.numeric(out$par[which.min(out$fvalue)])

    Xnew <- matrix(opt, nrow = 1)
    Ynew <- apply(Xnew, 1, y)


    if (i %in% c(1:5)) {
        ygrid <- c()
        zgrid <- c()
        for (j in 1:1000) {
            ygrid <- c(ygrid, integral(xgrid[j], X0, XN, nu.hat = nu.hat, theta = theta, noise = noise))
            zgrid <- c(zgrid, integral_derivative(xgrid[j], X0, XN, nu.hat = nu.hat, theta = theta, noise = noise))
        }
        plot(xgrid, ygrid, col = "black") +
            # lines(xgrid, zgrid, col = "blue") +
            abline(v = Xnew, col = "red")
    }


    XN <- rbind(XN, matrix(opt))
    Y <- c(Y, Ynew)

    # mod <- mleHomGP(
    #     X = XN, Z = Y,
    #     lower = c(0.0001, 0.0001), upper = c(1, 0.1),
    #     known = list(beta0 = 0),
    #     init = c(list(theta = c(0.1), g = 0.01))
    # )
    # print(c(mod$theta, mod$nu_hat, mod$g))
    # theta <- mod$theta
    # nu.hat <- mod$nu_hat

    fit <- mle_fit(XN, Y, cov_type)
    theta <- as.numeric(fit$theta)
    nu.hat <- fit$nu_hat
    noise <- fit$sigma2

    print(sprintf("iteration: %d, Xnew: %f, theta: %f, nu_hat: %f, noise: %f", i, Xnew, theta, nu.hat, noise))
}
print(Sys.time() - start)

# print(XN)

# par(mfrow = c(2, 1))

# Y <- apply(XN, 1, y)
pgrid <- c()
fgrid <- c()
for (i in 1:length(xgrid)) {
    fgrid <- c(fgrid, f(xgrid[i]))
    pgrid <- c(pgrid, predict_mean(xgrid[i], XN, Y, cov_type, nu.hat, theta, noise))
}
plt1 <- plot(xgrid, fgrid, col = "black") +
    lines(xgrid, pgrid, col = "blue") +
    abline(v = XN[1:N], col = "yellow") +
    abline(v = XN[(N + 1):(N + step)], col = "green", lty = 3) +
    abline(v = X0, col = "red")
# print(plt1)
print(rmse(fgrid, pgrid))


XN <- matrix(XN[1:N, ], nrow = N)
Y <- Y[1:N]

start <- Sys.time()
for (i in 1:step) {
    out <- IMSPE_optim(mod)
    Xnew <- matrix(out$par, nrow = 1)
    Ynew <- apply(Xnew, 1, y)
    # print(Xnew)

    XN <- rbind(XN, Xnew)
    Y <- c(Y, Ynew)

    mod <- update(mod, Xnew = Xnew, Znew = Ynew, ginit = mod$g * 1.01)
}
print(Sys.time() - start)

# Y <- apply(XN, 1, y)

pgrid <- c()
for (i in 1:length(xgrid)) {
    pgrid <- c(pgrid, predict_mean(xgrid[i], XN, Y, cov_type, nu.hat, theta, noise))
}
plt2 <- plot(xgrid, fgrid, col = "black") +
    lines(xgrid, pgrid, col = "blue") +
    abline(v = XN[1:N], col = "yellow") +
    abline(v = XN[(N + 1):(N + step)], col = "green", lty = 3) +
    abline(v = X0, col = "red")
print(rmse(fgrid, pgrid))

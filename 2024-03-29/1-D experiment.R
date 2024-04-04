library("optimx")
library("hetGP")
library("Metrics") # compute rmse
library("psych") # compute trace
library("BB")

setwd("D:/workspace/projects/R/hetGP")

f <- function(x) {
    return(sum(cos(8 * pi * x / (x + 0.2))) + sum(x))
}

# real function with constant noise = 0.01
y <- function(x) {
    return(f(x) + rnorm(1, sd = 0.01))
}

h <- function(x) {
    # return(1) # uniform(0, 1)
    # return(exp(-x^2 / 2)) # norm(0.5, 1)
    return(x^(-0.5) * (1 - x)^(-0.5)) # beta(0.5, 0.5)
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

Laplace_method <- function(par, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
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
N <- 10

# XN <- seq(0, 1, length.out = N)
# XN <- runif(N, 0, 1)
# XN <- rnorm(N, 0.5, 0.1)
XN <- rbeta(N, 0.5, 0.5)

XN <- sort(XN[0 <= XN & XN <= 1])
N <- length(XN)
print(sprintf("XN length: %d", N))
XN <- matrix(XN, nrow = N)

Y <- apply(XN, 1, y)

mod <- mleHomGP(X = XN, Z = Y, covtype = cov_type)
theta <- mod$theta
nu.hat <- mod$nu_hat
noise <- mod$g
print(sprintf("theta: %f, nu.hat: %f, noise: %f", theta, nu.hat, noise))


# #################### train Laplace
step <- 100
start <- Sys.time()
for (k in 1:step) {
    out <- optim(c(0.5),
        fn = Laplace_method,
        method = "Brent",
        lower = c(0), upper = c(1),
        XN = XN, cov_type = cov_type, nu.hat = nu.hat, theta = theta, noise = noise
    )
    # print(out)
    opt <- out$par

    Xnew <- matrix(c(opt), nrow = 1)
    Ynew <- apply(Xnew, 1, y)

    XN <- rbind(XN, Xnew)
    Y <- c(Y, Ynew)

    mod <- update(mod, Xnew = Xnew, Znew = Ynew, ginit = mod$g * 1.01)
    theta <- mod$theta
    nu.hat <- mod$nu_hat
    noise <- mod$g

    print(sprintf("iteration: %d, Xnew: %f, theta: %f, nu_hat: %f, noise: %f", k, Xnew, theta, nu.hat, noise))
}
print(sprintf("train time: %f", Sys.time() - start))


#################### visualization
par(mfrow = c(2, 1))

# xgrid <- seq(0, 1, length.out = 1000)
# xgrid <- runif(1000, 0, 1)
# xgrid <- rnorm(1000,0.5,0.1)
xgrid <- rbeta(1000, 0.5, 0.5)

xgrid <- sort(xgrid[0 <= xgrid & xgrid <= 1])

pgrid <- c()
fgrid <- c()
for (i in 1:length(xgrid)) {
    fgrid <- c(fgrid, f(xgrid[i]))
    pgrid <- c(pgrid, predict_mean(xgrid[i], XN, Y, cov_type, nu.hat, theta, noise))
}
plot(xgrid, fgrid, col = "black") +
    lines(xgrid, pgrid, col = "blue") +
    abline(v = XN[(N + 1):(N + step)], col = "green", lty = 3) +
    abline(v = XN[1:N], col = "red")

print(rmse(fgrid, pgrid))


##################### train IMSPE
XN <- matrix(XN[1:N, ], nrow = N)
Y <- Y[1:N]

start <- Sys.time()
for (i in 1:step) {
    out <- IMSPE_optim(mod, h = 0)
    Xnew <- matrix(out$par, nrow = 1)
    Ynew <- apply(Xnew, 1, y)
    # print(Xnew)

    XN <- rbind(XN, Xnew)
    Y <- c(Y, Ynew)

    mod <- update(mod, Xnew = Xnew, Znew = Ynew, ginit = mod$g * 1.01)
    print(sprintf("iteration: %d, Xnew: %f, theta: %f, nu_hat: %f, noise: %f", i, Xnew, mod$theta, mod$nu_hat, mod$g))
}
print(sprintf("train time: %f", Sys.time() - start))

pgrid <- c()
for (i in 1:length(xgrid)) {
    pgrid <- c(pgrid, predict_mean(xgrid[i], XN, Y, cov_type, nu.hat, theta, noise))
}
plot(xgrid, fgrid, col = "black") +
    lines(xgrid, pgrid, col = "blue") +
    abline(v = XN[(N + 1):(N + step)], col = "green", lty = 3) +
    abline(v = XN[1:N], col = "red")
print(rmse(fgrid, pgrid))

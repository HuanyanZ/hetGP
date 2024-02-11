library("hetGP")
library("DEoptim")
library("Metrics")

# real function
f <- function(x) {
    return(f1d2(sum(x)))
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

# mse_N_1 = mse_N + diff
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
}

sample_optim <- function(XN, Xnew, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    xgrid <- seq(0, 1, length.out = 2000)
    vgrid <- c()
    for (i in 1:length(xgrid)) {
        vgrid <- c(vgrid, diff(xgrid[i], XN, Xnew, cov_type, nu.hat, theta, noise))
    }
    return(xgrid[which.min(vgrid)])
}

N <- 5
cov_type <- "Gaussian"
nu.hat <- 1
theta <- 0.1
noise <- 0.01

XN <- sort(rnorm(N, mean = 0.5))
XN <- (XN - min(XN)) / (max(XN) - min(XN))
XN <- matrix(XN, nrow = N)
# XN <- rbind(XN, XN)
# print(XN)

Y <- apply(XN, 1, y)


xgrid <- seq(0, 1, length.out = 2000)
fgrid <- c()
mNgrid <- c()
pgrid <- c()
for (i in 1:length(xgrid)) {
    fgrid <- c(fgrid, f(xgrid[i]))
    mNgrid <- c(mNgrid, mse(xgrid[i], XN, cov_type, nu.hat, theta, noise))
    pgrid <- c(pgrid, predict_mean(xgrid[i], XN, Y, cov_type, nu.hat, theta, noise))
}

Xnew <- xgrid[which.max(mNgrid)]
# print(Xnew)


par(mfrow = c(5, 1))
plt1 <- plot(xgrid, fgrid, col = "black", ylim = c(-2, 1)) +
    lines(xgrid, pgrid, col = "blue") +
    abline(v = Xnew, col = "red")
legend(0, 1,
    legend = c("real value", "predict mean", "argmaxMSE(x)"),
    col = c("black", "blue", "red"), lty = 1:1:1, cex = 0.8
)

plt2 <- plot(xgrid, mNgrid, col = "black", ylim = c(0, 0.5)) +
    abline(v = XN, col = "blue") +
    abline(v = Xnew, col = "red")
legend(0, 0.5,
    legend = c("MSE(x)", "current samples", "argmaxMSE(x)"),
    col = c("black", "blue", "red"), lty = 1:1:1, cex = 0.8
)

mN1grid <- c()
dgrid <- c()
XN_new <- rbind(XN, matrix(Xnew))
for (i in 1:length(xgrid)) {
    mN1grid <- c(mN1grid, mse(xgrid[i], XN_new, cov_type, nu.hat, theta, noise))

    d <- diff(xgrid[i], XN, Xnew, cov_type, nu.hat, theta, noise)
    dgrid <- c(dgrid, d)
}

Xnew <- xgrid[which.max(mN1grid)]
# print(Xnew)

plt3 <- plot(xgrid, mN1grid, col = "black", ylim = c(-0.5, 0.2)) +
    lines(xgrid, mNgrid + dgrid, col = "green") +
    lines(xgrid, dgrid, col = "blue") +
    abline(v = Xnew, col = "red")
legend(0, 0.2,
    legend = c("MSE(x)", "diff", "argmaxMSE(x)"),
    col = c("green", "blue", "red"), lty = 1:1:1, cex = 0.8
)

mod <- mleHetGP(X = XN, Z = Y, covtype = cov_type)


# par(mfrow = c(2, 1))

muX <- 0.5
step <- 30
print(Sys.time())
for (i in 1:step) {
    out1 <- DEoptim(
        diff,
        lower = c(0), upper = c(1),
        DEoptim.control(trace = FALSE, itermax = 30),
        XN = XN, Xnew = muX,
    )
    opt1 <- out1$optim$bestmem
    # opt2 <- sample_optim(XN, muX, cov_type, nu.hat, theta, noise)
    # print(c(opt1, opt2))
    XN <- rbind(XN, matrix(opt1))
}
print(Sys.time())
# print(XN)
Y <- apply(XN, 1, y)

pgrid <- c()
for (i in 1:length(xgrid)) {
    pgrid <- c(pgrid, predict_mean(xgrid[i], XN, Y, cov_type, nu.hat, theta, noise))
}
plt4 <- plot(xgrid, fgrid, col = "black", ylim = c(-2, 1)) +
    lines(xgrid, pgrid, col = "blue")
legend(0, 1,
    legend = c("real value", "predict mean"),
    col = c("black", "blue"), lty = 1:1, cex = 0.8
)
print(rmse(fgrid, pgrid))


XN <- XN[1:N, ]
# print(XN)
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
plt5 <- plot(xgrid, fgrid, col = "black", ylim = c(-2, 1)) +
    lines(xgrid, pgrid, col = "blue")
legend(0, 1,
    legend = c("real value", "predict mean"),
    col = c("black", "blue"), lty = 1:1, cex = 0.8
)
print(rmse(fgrid, pgrid))

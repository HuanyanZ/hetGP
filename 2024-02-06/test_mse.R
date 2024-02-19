library(hetGP)
library(gridExtra)
library(ggplot2)

mse <- function(par, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    x0 <- matrix(par, nrow = 1, byrow = TRUE)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type)

    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1]))

    L <- t(chol(K))
    v <- solve(L, t(KX))

    mse <- nu.hat - t(v) %*% v
    return(mse)
}

mse2 <- function(X0, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    K <- c()
    for (i in 1:length(XN)) {
        for (j in 1:length(XN)) {
            K <- c(K, -(XN[i, ] - XN[j, ])^2 / theta)
        }
    }
    K <- matrix(nu.hat * exp(K), nrow = length(XN)) + diag(rep(noise, length(XN)))

    Kn <- c()
    for (i in 1:length(XN)) {
        Kn <- c(Kn, -(X0 - XN[i, ])^2 / theta)
    }

    Kn <- matrix(nu.hat * exp(Kn), ncol = 1)
    mse2 <- nu.hat - t(Kn) %*% solve(K) %*% Kn
    return(mse2)
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

N <- 5
cov_type <- "Gaussian"
nu.hat <- 1
theta <- 0.1
noise <- 0.001

XN <- matrix(runif(N, 0, 1), nrow = N)
# XN <- rbind(XN, XN)
print(XN)

Y <- apply(XN, 1, y)

X0 <- runif(1, 0, 1)
m1 <- mse(X0, XN, cov_type, nu.hat, theta, noise)
m2 <- mse2(X0, XN, cov_type, nu.hat, theta, noise)
print(c(m1, m2))


# grid plot of mse function
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
print(Xnew)


par(mfrow = c(3, 1))
plt1 <- plot(xgrid, fgrid, col = "black", ylim = c(-2, 1)) +
    lines(xgrid, pgrid, col = "blue") +
    abline(v = Xnew, col = "red")
legend(0, 1,
    legend = c("real value", "predict mean", "argmaxMSE(x)"),
    col = c("black", "blue", "red"), lty = 1:1:1, cex = 0.8
)

plt2 <- plot(xgrid, mNgrid, col = "black", ylim = c(0, 0.2)) +
    abline(v = XN, col = "blue") +
    abline(v = Xnew, col = "red")
legend(0, 0.2,
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
print(Xnew)

plt3 <- plot(xgrid, mN1grid, col = "black", ylim = c(-0.05, 0.05)) +
    lines(xgrid, mNgrid + dgrid, col = "green") +
    lines(xgrid, dgrid, col = "blue") +
    abline(v = Xnew, col = "red")

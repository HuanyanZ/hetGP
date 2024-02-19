library("hetGP")
library("DEoptim")

f <- function(x) {
    return(f1d2(sum(x)))
}

y <- function(x) {
    # set.seed(123)
    return(f(x) + rnorm(1, sd = 0.01))
}

mse <- function(par, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    x0 <- matrix(par, nrow = 1, byrow = TRUE)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type)

    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1]))
    # K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type)
    # print(K)

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
}

diff2 <- function(X0, XN, Xnew, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    x0 <- matrix(X0, nrow = 1, byrow = TRUE)
    xN1 <- matrix(Xnew, nrow = 1, byrow = TRUE)

    sigma2 <- mse(Xnew, XN, cov_type, nu.hat, theta, noise)
    sigma2 <- sigma2 + noise

    K1 <- nu.hat * cov_gen(X1 = x0, X2 = xN1, theta = theta, type = cov_type)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type)
    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1]))
    Kn <- nu.hat * cov_gen(X1 = xN1, X2 = XN, theta = theta, type = cov_type)

    return(-1 / sigma2 * (K1 - KX %*% solve(K) %*% t(Kn))^2)
}

N <- 5
cov_type <- "Gaussian"
nu.hat <- 1
theta <- 0.1
noise <- 0.01

XN <- matrix(runif(N, 0, 1), nrow = N)
Y <- apply(XN, 1, y)

X0 <- runif(1, 0, 1)

xgrid <- seq(0, 1, length.out = 2000)
ygrid <- c()
for (i in 1:length(xgrid)) {
    ygrid <- c(ygrid, mse(xgrid[i], XN, cov_type, nu.hat, theta, noise))
}
Xnew <- xgrid[which.max(ygrid)]

d1 <- diff(X0, XN, Xnew, cov_type, nu.hat, theta, noise)
d2 <- diff2(X0, XN, Xnew, cov_type, nu.hat, theta, noise)
XN_new <- rbind(XN, matrix(Xnew))
mN_1 <- mse(X0, XN_new, cov_type, nu.hat, theta, noise)
mN <- mse(X0, XN, cov_type, nu.hat, theta, noise)
d3 <- mN_1 - mN

print(c(d1, d2, d3))
print(c(mN_1, mN, d3))

dgrid <- c()
for (i in 1:length(xgrid)) {
    d <- diff(xgrid[i], XN, Xnew, cov_type, nu.hat, theta, noise)

    dgrid <- c(dgrid, d)
}

muX <- 0.5
out1 <- DEoptim(
    diff,
    lower = c(0), upper = c(1),
    DEoptim.control(trace = FALSE),
    XN = XN, Xnew = muX,
)
opt1 <- out1$optim$bestmem
print(opt1)

out2 <- optim(
    par = c(0.5),
    fn = diff,
    lower = c(0), upper = c(1),
    method = "L-BFGS-B",
    XN = XN, Xnew = muX,
)
opt2 <- out2$par
print(opt2)

plot(xgrid, dgrid) +
    abline(v = opt1, col = "red") +
    abline(v = opt2, col = "blue")

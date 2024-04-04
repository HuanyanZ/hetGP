library(rgl)
library(hetGP)
library(optimx)

f <- function(x) {
    return(sum(cos(8 * pi * x / (x + 0.2))) + sum(x))
}

# real function with constant noise = 0.01
y <- function(x) {
    return(f(x) + rnorm(1, sd = 0.01))
}


h <- function(x) {
    # return(1) # uniform(0, 1)
    return(exp(-x^2 / 2)) # norm(0.5, 1)
    # return((x + 10^(-9))^(-0.5) * (1 - x + 10^(-9))^(-0.5)) # beta(0.5, 0.5)
}

g <- function(XN1, X, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    x <- matrix(X, nrow = 1, byrow = TRUE)
    xN1 <- matrix(XN1, nrow = 1, byrow = TRUE)

    K1 <- nu.hat * cov_gen(X1 = x, X2 = xN1, theta = theta, type = cov_type) # k(x, xn+1)
    KX <- nu.hat * cov_gen(X1 = x, X2 = XN, theta = theta, type = cov_type) # kn(x)^T
    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1])) # Kn
    Kn <- nu.hat * cov_gen(X1 = xN1, X2 = XN, theta = theta, type = cov_type) # kn(xn+1)^T

    L <- t(chol(K)) # K = L L^T
    u <- solve(L, t(Kn)) # L^(-1) kn(xn+1)

    v0 <- solve(L, t(KX)) # L^(-1) kn(x)
    # v1 <- solve(L, t((X[i] - XN[, i]) * KX)) # L^(-1) diag(X[i] - XN[,i]) kn(x)
    # v2 <- solve(L, t((X[i] - XN[, i])^2 * KX)) # L^(-1) diag(X[i] - XN[,i])^2 kn(x)

    g <- K1 - t(u) %*% v0
    return(g)
}

gd1 <- function(i, XN1, X, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    x <- matrix(X, nrow = 1, byrow = TRUE)
    xN1 <- matrix(XN1, nrow = 1, byrow = TRUE)

    K1 <- nu.hat * cov_gen(X1 = x, X2 = xN1, theta = theta, type = cov_type) # k(x, xn+1)
    KX <- nu.hat * cov_gen(X1 = x, X2 = XN, theta = theta, type = cov_type) # kn(x)
    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1])) # Kn
    Kn <- nu.hat * cov_gen(X1 = xN1, X2 = XN, theta = theta, type = cov_type) # kn(xn+1)

    L <- t(chol(K)) # K = L L^T
    u <- solve(L, t(Kn)) # L^(-1) kn(xn+1)

    v0 <- solve(L, t(KX)) # L^(-1) kn(x)
    v1 <- solve(L, t((X[i] - XN[, i]) * KX)) # L^(-1) diag(X[i] - XN[,i]) kn(x)
    v2 <- solve(L, t((X[i] - XN[, i])^2 * KX)) # L^(-1) diag(X[i] - XN[,i])^2 kn(x)

    gdi <- -2 / theta * (X[i] - XN1[i]) * K1 + 2 / theta * t(u) %*% v1
    return(gdi)
}

gd2 <- function(i, j, XN1, X, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    x <- matrix(X, nrow = 1, byrow = TRUE)
    xN1 <- matrix(XN1, nrow = 1, byrow = TRUE)

    K1 <- nu.hat * cov_gen(X1 = x, X2 = xN1, theta = theta, type = cov_type) # k(x, xn+1)
    KX <- nu.hat * cov_gen(X1 = x, X2 = XN, theta = theta, type = cov_type) # kn(x)
    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1])) # Kn
    Kn <- nu.hat * cov_gen(X1 = xN1, X2 = XN, theta = theta, type = cov_type) # kn(xn+1)

    L <- t(chol(K)) # K = L L^T
    u <- solve(L, t(Kn)) # L^(-1) kn(xn+1)

    v0 <- solve(L, t(KX)) # L^(-1) kn(x)
    v1 <- solve(L, t((X[i] - XN[, i]) * KX)) # L^(-1) diag(X[i] - XN[,i]) kn(x)
    v2 <- solve(L, t((X[i] - XN[, i])^2 * KX)) # L^(-1) diag(X[i] - XN[,i])^2 kn(x)
    v3 <- solve(L, t((X[i] - XN[, i]) * (X[j] - XN[, j]) * KX)) # L^(-1) diag((X[i] - XN[,i])(X[j] - XN[, j])) kn(x)

    if (i == j) {
        gdii <- -2 / theta * K1 + 4 / theta^2 * (X[i] - XN1[i])^2 * K1 + 2 / theta * t(u) %*% v0 - 4 / theta^2 * t(u) %*% v2
        return(gdii)
    }
    gdij <- 4 / theta^2 * (X[i] - XN1[i]) * (X[j] - XN1[j]) * K1 - 4 / theta^2 * t(u) %*% v3
    return(gdij)
}

hess_g <- function(XN1, X, XN, cov_type = "Gaussian", nu.hat = 1, theta = 0.1, noise = 0.01) {
    d <- length(X)
    H <- matrix(rep(0, d^2), nrow = d)
    for (i in 1:d) {
        for (j in i:d) {
            H[i, j] <- gd2(i, j, XN1, X, XN, cov_type, nu.hat, theta, noise)
            H[j, i] <- H[i, j]
        }
    }
    return(H)
}

Laplace_method_mv <- function(par, XN, cov_type = "Gaussian", nu.hat = 1, theta = c(0.5, 0.5), noise = 0.01) {
    d <- ncol(XN)

    X <- matrix(par, nrow = 1) # x
    XN1 <- matrix(par, nrow = 1) # xn+1

    K1 <- nu.hat * cov_gen(X1 = X, X2 = XN1, theta = theta, type = cov_type) # k(x, xn+1)
    KX <- nu.hat * cov_gen(X1 = X, X2 = XN, theta = theta, type = cov_type) # kn(x)
    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1])) # Kn
    Kn <- nu.hat * cov_gen(X1 = XN1, X2 = XN, theta = theta, type = cov_type) # kn(xn+1)

    L <- t(chol(K)) # K = L L^T
    u <- solve(L, t(Kn)) # L^(-1) kn(xn+1)

    v0 <- solve(L, t(KX)) # L^(-1) kn(x)

    sigma2 <- nu.hat + noise - t(u) %*% u

    g <- K1 - t(u) %*% v0

    gd <- rep(0, d)
    H <- matrix(rep(0, d^2), nrow = d)
    for (k in 1:d) {
        v1 <- solve(L, t((X[k] - XN[, k]) * KX)) # L^(-1) diag(X[k] - XN[,k]) kn(x)
        gd <- -2 / theta[k] * (X[k] - XN1[k]) * K1 + 2 / theta[k] * t(u) %*% v1
    }

    for (i in 1:d) {
        gdi <- -2 / theta[i] * (X[i] - XN1[i]) * K1 + 2 / theta[i] * t(u) %*% v1

        for (j in i:d) {
            if (i == j) {
                v2 <- solve(L, t((X[i] - XN[, i])^2 * KX)) # L^(-1) diag(X[i] - XN[,i])^2 kn(x)
                gdii <- -2 / theta[i] * K1 + 4 / theta[i]^2 * (X[i] - XN1[i])^2 * K1 + 2 / theta[i] * t(u) %*% v0 - 4 / theta[i]^2 * t(u) %*% v2

                H[i, j] <- -(gdi / g)^2 + gdii / g
            } else {
                gdj <- -2 / theta[j] * (X[j] - XN1[j]) * K1 + 2 / theta[j] * t(u) %*% v1

                v3 <- solve(L, t((X[i] - XN[, i]) * (X[j] - XN[, j]) * KX)) # L^(-1) diag((X[i] - XN[,i])(X[j] - XN[, j])) kn(x)

                gdij <- 4 / theta[j]^2 * (X[i] - XN1[i]) * (X[j] - XN1[j]) * K1 - 4 / theta[j]^2 * t(u) %*% v3

                H[j, i] <- -gdi * gdj / g^2 + gdij / g
                H[j, i] <- H[i, j]
            }
        }
    }

    I <- -1 / sigma2 * pi^(d / 2) * prod(apply(X, 2, h)) * g^2 / sqrt(as.numeric(determinant(H, logarithm = FALSE)$modulus))
    return(I)
}

predict_mean <- function(X0, XN, Y, cov_type = "Gaussian", nu.hat = 1, theta = c(0.5, 0.5), noise = 0.01) {
    x0 <- matrix(X0, nrow = 1)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type)

    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(rep(noise, dim(XN)[1]))

    L <- t(chol(K))
    v <- solve(L, t(KX))
    u <- solve(L, Y)
    return(t(v) %*% u)
}

d <- 2
N <- 20
cov_type <- "Gaussian"

# x1 <- runif(N, 0, 1)
# x2 <- runif(N, 0, 1)

x1 <- rnorm(N, 0.5, 0.1)
x2 <- rnorm(N, 0.5, 0.1)
print(c(length(x1[0 <= x1 & x1 <= 1]), length(x2[0 <= x2 & x2 <= 1])))

# x1 <- rbeta(N, 0.5, 0.5)
# x2 <- rbeta(N, 0.5, 0.5)

# XN <- matrix(x1, nrow = N)
XN <- cbind(x1, x2)
print(XN)

Y <- apply(XN, 1, y)


################# validation
# delta <- 1e-9
# X <- c(runif(1, 0, 1), runif(1, 0, 1))
# XN1 <- c(runif(1, 0, 1), runif(1, 0, 1))
# print(c((g(XN1, X + c(delta, 0), XN) - g(XN1, X, XN)) / delta, gd1(1, XN1, X, XN)))
# print(c((gd1(2, XN1, X + c(delta, 0), XN) - gd1(2, XN1, X, XN)) / delta, gd2(2, 1, XN1, X, XN)))
# print(c((gd1(1, XN1, X + c(delta, 0), XN) - gd1(1, XN1, X, XN)) / delta, gd2(1, 1, XN1, X, XN)))

# H <- hess_g(XN1, X, XN)
# print(H)

# print(Laplace_method_mv(XN1, XN))


# ################ visualization
# x1 <- runif(5000, 0, 1)
# x2 <- runif(5000, 0, 1)
# xgrid <- cbind(x1, x2)
# z <- apply(xgrid, 1, Laplace_method_mv, XN = XN)

# print(sprintf("grid min: (%f, %f)", x1[which.min(z)], x2[which.min(z)]))

# open3d()
# plot3d(x1, x2, z)


# ################# optim
# out <- optim(
#     rep(0.5, d),
#     Laplace_method_mv,
#     method = "L-BFGS-B",
#     lower = rep(0, d), upper = rep(1, d),
#     XN = XN
# )
# print(sprintf("optim min: (%f, %f)", out$par[1], out$par[2]))


# ################# multistart
multi_input <- matrix(rep(0, d * 4^d), nrow = 4^d)
for (c in 1:d) {
    multi_input[, c] <- rep(kronecker(seq(0, 1, length.out = 4), rep(1, 4^(d - c))), 4^(c - 1))
}

# out <- multistart(multi_input,
#     Laplace_method_mv,
#     method = "L-BFGS-B",
#     lower = rep(0, d), upper = rep(1, d),
#     XN = XN
# )

# print(sprintf("multistart min: (%f, %f)", out$p1[which.min(out$value)], out$p2[which.min(out$value)]))


################# train Laplace
mod <- mleHomGP(X = XN, Z = Y, covtype = cov_type)
theta <- mod$theta
nu.hat <- mod$nu_hat
noise <- mod$g
print(sprintf("initial params, theta: (%f, %f), nu.hat: %f, noise: %f", theta[1], theta[2], nu.hat, noise))

step <- 50
start <- Sys.time()
for (k in 1:step) {
    out <- multistart(
        multi_input,
        Laplace_method_mv,
        method = "L-BFGS-B",
        lower = rep(0, d), upper = rep(1, d),
        XN = XN, cov_type = cov_type, nu.hat = nu.hat, theta = theta, noise = noise
    )

    Xnew <- matrix(c(out$p1[which.min(out$value)], out$p2[which.min(out$value)]), nrow = 1)
    Ynew <- apply(Xnew, 1, y)

    XN <- rbind(XN, Xnew)
    Y <- c(Y, Ynew)

    mod <- update(mod, Xnew = Xnew, Znew = Ynew, ginit = mod$g * 1.01)
    theta <- mod$theta
    nu.hat <- mod$nu_hat
    noise <- mod$g

    print(sprintf("iteration: %d, Xnew: (%f, %f), theta: (%f, %f), nu_hat: %f, noise: %f", k, Xnew[1], Xnew[2], theta[1], theta[2], nu.hat, noise))
}
print(sprintf("train time: %f", Sys.time() - start))

par(mfrow = c(2, 1))

plot(XN[, 1], XN[, 2]) + points(XN[1:N, 1], XN[1:N, 2], col = "red")

# x1 <- runif(5000, 0, 1)
# x2 <- runif(5000, 0, 1)

x1 <- rnorm(5000, 0.5, 0.1)
x2 <- rnorm(5000, 0.5, 0.1)

# x1 <- rbeta(5000, 0.5, 0.5)
# x2 <- runif(5000, 0.5, 0.5)

xgrid <- cbind(x1, x2)
z <- apply(xgrid, 1, f)
pm1 <- apply(xgrid, 1, predict_mean, XN = XN, Y = Y, cov_type = cov_type, nu.hat = nu.hat, theta = theta, noise = noise)

open3d()
plot3d(x1, x2, z)
rgl.snapshot("2d-norm-sample.png", fmt = "png")

open3d()
plot3d(x1, x2, pm1)
rgl.snapshot("2d-norm-laplace.png", fmt = "png")

print(sprintf("rmse: %f", mean((pm1 - z)^2)))

################# train IMSPE
XN <- matrix(XN[1:N, ], nrow = N)
Y <- Y[1:N]

start <- Sys.time()
for (k in 1:step) {
    out <- IMSPE_optim(mod, h = 2)

    Xnew <- matrix(out$par, nrow = 1)
    Ynew <- apply(Xnew, 1, y)
    XN <- rbind(XN, Xnew)
    Y <- c(Y, Ynew)

    mod <- update(mod, Xnew = Xnew, Znew = Ynew, ginit = mod$g * 1.01)
    theta <- mod$theta
    nu.hat <- mod$nu_hat
    noise <- mod$g

    print(sprintf("iteration: %d, Xnew: (%f, %f), theta: (%f, %f), nu_hat: %f, noise: %f", k, Xnew[1], Xnew[2], theta[1], theta[2], nu.hat, noise))
}
print(sprintf("train time: %f", Sys.time() - start))

plot(XN[, 1], XN[, 2]) + points(XN[1:N, 1], XN[1:N, 2], col = "red")

pm2 <- apply(xgrid, 1, predict_mean, XN = XN, Y = Y, cov_type = cov_type, nu.hat = nu.hat, theta = theta, noise = noise)
print(sprintf("rmse: %f", mean((pm2 - z)^2)))

open3d()
plot3d(x1, x2, pm2)
rgl.snapshot("2d-norm-imspe.png", fmt = "png")

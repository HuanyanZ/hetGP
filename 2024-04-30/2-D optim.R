library(rgl)
library(hetGP)
library(optimx)
library(mvtnorm)
library(DiceDesign)

############## benchmark functions
sphere_func <- function(x) {
    return(sum(x^2))
}

ackley_func <- function(x) {
    return(-20 * exp(-0.2 * sqrt(0.5 * sum(x^2))) - exp(0.5 * sum(cos(2 * pi * x))) + exp(1) + 20)
}

rastrigin_func <- function(x, A = 10) {
    n <- length(x)
    return(A * n + sum(x^2) - A * sum(cos(2 * pi * x)))
}

gaussian_pdf_func <- function(x) {
    f1 <- dmvnorm(x, mean = c(0.25, 0.75), sigma = diag(rep(0.01, d)))
    f2 <- dmvnorm(x, mean = c(0.75, 0.25), sigma = diag(rep(0.01, d)))

    return(f1 - f2)
}

rastrigin_adj_func <- function(x) {
    return(sum(cos(8 * pi * x / (x + 0.2))) + sum(x^2))
}

# real function with constant noise = 0.01
y <- function(x, f) {
    return(f(x) + rnorm(1, sd = 0.01))
}


# pdf of input covariate X
h <- function(x, dist) {
    if (dist == "unif") {
        return(1) # uniform(0, 1)
    }
    if (dist == "norm") {
        return(exp(-(x - 0.5)^2 / (2 * 0.1^2))) # norm(0.5, 0.1)
    }
    if (dist == "beta") {
        return((x^0.5) * ((1 - x)^3)) # beta(1.5, 4)
    }
    return(1)
}

# Approximation of objective function based on Laplace method
Laplace_method_expl <- function(par, mod, dist = "unif") {
    Xn <- mod$X0
    mult <- mod$mult
    theta <- mod$theta
    cov.type <- mod$covtype
    nu.hat <- mod$nu_hat
    noise <- mod$g

    d <- ncol(Xn)

    X <- matrix(par, nrow = 1) # x
    XN1 <- matrix(par, nrow = 1) # xn1

    K1 <- nu.hat * cov_gen(X1 = X, X2 = XN1, theta = theta, type = cov.type) # k(x, xn+1)
    KX <- nu.hat * cov_gen(X1 = X, X2 = Xn, theta = theta, type = cov.type) # kn(x)
    K <- nu.hat * cov_gen(X1 = Xn, theta = theta, type = cov.type) + diag(rep(noise, nrow(Xn)) / mult) # Kn
    Kn <- nu.hat * cov_gen(X1 = XN1, X2 = Xn, theta = theta, type = cov.type) # kn(xn+1)

    L <- t(chol(K)) # K = L L^T
    u <- solve(L, t(Kn)) # L^(-1) kn(xn+1)

    v0 <- solve(L, t(KX)) # L^(-1) kn(x)

    sigma2 <- nu.hat + noise - t(u) %*% u

    g <- K1 - t(u) %*% v0

    H <- matrix(rep(0, d^2), nrow = d)

    for (i in 1:d) {
        v1 <- solve(L, t((X[i] - Xn[, i]) * KX)) # L^(-1) diag(X[i] - Xn[,i]) kn(x)

        gdi <- -2 / theta * (X[i] - XN1[i]) * K1 + 2 / theta * t(u) %*% v1

        for (j in i:d) {
            if (i == j) {
                v2 <- solve(L, t((X[i] - Xn[, i])^2 * KX)) # L^(-1) diag(X[i] - Xn[,i])^2 kn(x)
                gdii <- -2 / theta * K1 + 4 / theta^2 * (X[i] - XN1[i])^2 * K1 + 2 / theta * t(u) %*% v0 - 4 / theta^2 * t(u) %*% v2

                H[i, j] <- -(gdi / g)^2 + gdii / g
            } else {
                gdj <- -2 / theta * (X[j] - XN1[j]) * K1 + 2 / theta * t(u) %*% v1

                v3 <- solve(L, t((X[i] - Xn[, i]) * (X[j] - Xn[, j]) * KX)) # L^(-1) diag((X[i] - Xn[,i])(X[j] - Xn[, j])) kn(x)

                gdij <- 4 / theta^2 * (X[i] - XN1[i]) * (X[j] - XN1[j]) * K1 - 4 / theta^2 * t(u) %*% v3

                H[j, i] <- -gdi * gdj / g^2 + gdij / g
                H[j, i] <- H[i, j]
            }
        }
    }

    I <- -1 / sigma2 * pi^(d / 2) * prod(apply(X, 2, h, dist = dist)) * g^2 / sqrt(as.numeric(determinant(H, logarithm = FALSE)$modulus))
    return(I)
}

reps_optim <- function(mod, dist = "unif") {
    Xn <- mod$X0
    mult <- mod$mult
    theta <- mod$theta
    cov.type <- mod$covtype
    nu.hat <- mod$nu_hat
    noise <- mod$g

    d <- ncol(Xn)
    n <- nrow(Xn)

    S <- nu.hat * cov_gen(X1 = Xn, theta = theta, type = cov.type) + diag(rep(noise, n) / mult) # Sigma_n

    L <- t(chol(S)) # S = L L^T
    L.inv <- solve(L)

    S.inv <- t(L.inv) %*% L.inv

    deltas <- c()
    for (k in 1:n) {
        xk <- matrix(Xn[k, ], nrow = 1)
        uk <- S.inv[, k]

        hk <- prod(apply(xk, 2, h, dist = dist))

        Knk <- nu.hat * cov_gen(X1 = xk, X2 = Xn, theta = theta, type = cov.type) # kn(xk, Xi)
        g <- uk %*% t(Knk)

        H <- matrix(rep(0, d^2), nrow = d)
        for (i in 1:d) {
            gdi <- -2 / theta * ((xk[i] - Xn[, i]) * Knk) %*% uk
            for (j in i:d) {
                if (i == j) {
                    gdii <- 4 / theta^2 * ((xk[i] - Xn[, i])^2 * Knk) %*% uk

                    H[i, j] <- -(gdi / g)^2 + gdii / g
                } else {
                    gdj <- -2 / theta * ((xk[j] - Xn[, j]) * Knk) %*% uk

                    gdij <- 4 / (theta * theta) * ((xk[i] - Xn[, i]) * (xk[j] - Xn[, j]) * Knk) %*% uk

                    H[j, i] <- -gdi * gdj / g^2 + gdij / g
                    H[j, i] <- H[i, j]
                }
            }
        }


        dk <- -1 / (mult[k] * (mult[k] + 1) / noise - S.inv[k, k]) * pi^(d / 2) * hk * g^2 / sqrt(as.numeric(determinant(H, logarithm = FALSE)$modulus))
        deltas <- c(deltas, dk)
    }

    return(list(value = deltas))
}

delta_func <- function(x, k, mod) {
    Xn <- mod$X0
    mult <- mod$mult
    theta <- mod$theta
    cov.type <- mod$covtype
    nu.hat <- mod$nu_hat
    noise <- mod$g

    Knx <- nu.hat * cov_gen(X1 = matrix(x, nrow = 1), X2 = Xn, theta = theta, type = cov.type) # k(x, Xn)
    S <- nu.hat * cov_gen(X1 = Xn, theta = theta, type = cov.type) + diag(rep(noise, nrow(Xn)) / mult) # Sigma_n

    L <- t(chol(S)) # S = L L^T
    L.inv <- solve(L)

    S.inv <- t(L.inv) %*% L.inv
    uk <- S.inv[, k]

    return((uk %*% t(Knx))^2 / (mult[k] * (mult[k] + 1) / noise - S.inv[k, k]))
}

# mean of GP conditional distribution
predict_mean <- function(X0, XN, Y, cov.type = "Gaussian", nu.hat = 1, theta = c(0.5, 0.5), noise = 0.01) {
    x0 <- matrix(X0, nrow = 1)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov.type)

    K <- nu.hat * cov_gen(X1 = XN, theta = theta, type = cov.type) + diag(rep(noise, dim(XN)[1]))

    L <- t(chol(K))
    v <- solve(L, t(KX))
    u <- solve(L, Y)
    return(t(v) %*% u)
}

train_Laplace <- function(XN, Y, mod, f, step = 50, dist = "unif") {
    start <- Sys.time()
    for (k in 1:step) {
        explore_out <- multistart(
            Xstart,
            Laplace_method_expl,
            method = "L-BFGS-B",
            lower = rep(0, d), upper = rep(1, d),
            mod = mod, dist = dist
        )
        explore_min <- min(explore_out$value)

        reps_out <- reps_optim(mod, dist)
        reps_min <- min(reps_out$value)

        if (explore_min < reps_min) {
            Xnew <- matrix(c(explore_out$p1[which.min(explore_out$value)], explore_out$p2[which.min(explore_out$value)]), nrow = 1)
            is_replication <- FALSE
        } else {
            Xnew <- matrix(mod$X0[which.min(reps_out$value), ], nrow = 1)
            is_replication <- TRUE
        }
        Ynew <- apply(Xnew, 1, y, f = f)

        XN <- rbind(XN, Xnew)
        Y <- c(Y, Ynew)

        mod <- update(mod, Xnew = Xnew, Znew = Ynew, ginit = mod$g * 1.01)
        print(sprintf("iteration: %d, is_replication: %i, Xnew: (%f, %f), theta: %f, nu_hat: %f, noise: %f", k, is_replication, Xnew[1], Xnew[2], mod$theta, mod$nu_hat, mod$g))
    }
    print(sprintf("Laplace train time: %f", Sys.time() - start))

    return(list(X_seq = XN, Y_seq = Y, nu.hat = mod$nu_hat, theta = mod$theta, noise = mod$g))
}

train_IMSPE <- function(XN, Y, mod, f, step = 50) {
    start <- Sys.time()
    for (k in 1:step) {
        out <- IMSPE_optim(mod, h = 2)

        n <- nrow(mod$X0)

        Xnew <- matrix(out$par, nrow = 1)
        Ynew <- apply(Xnew, 1, y, f = f)

        XN <- rbind(XN, Xnew)
        Y <- c(Y, Ynew)

        mod <- update(mod, Xnew = Xnew, Znew = Ynew, ginit = mod$g * 1.01)

        is_replication <- nrow(mod$X0) == n

        print(sprintf("iteration: %d, is_replication: %i, Xnew: (%f, %f), theta: %f, nu_hat: %f, noise: %f", k, is_replication, Xnew[1], Xnew[2], mod$theta, mod$nu_hat, mod$g))
    }
    print(sprintf("IMSE train time: %f", Sys.time() - start))

    return(list(X_seq = XN, Y_seq = Y, nu.hat = nu.hat, theta = theta, noise = noise))
}

run <- function(XN, f, dist, step = 50) {
    Y <- apply(XN, 1, y, f)

    # sample path of real function
    z <- apply(xgrid, 1, f)

    # fit model parameters by MLE
    mod <- mleHomGP(
        X = XN, Z = Y,
        covtype = cov.type,
        lower = c(1e-4, 1e-4), upper = c(1, 0.5),
        known = list(beta0 = 0),
        init = c(list(theta = c(0.1), g = 0.01))
    )
    print(mod)

    res_Laplace <- train_Laplace(XN, Y, mod, f, step, dist = dist)
    res_IMSPE <- train_IMSPE(XN, Y, mod, f, step)

    # predict mean of GP trained by Laplace method
    pm1 <- apply(xgrid, 1, predict_mean, XN = res_Laplace$X_seq, Y = res_Laplace$Y_seq, cov.type = cov.type, nu.hat = res_Laplace$nu.hat, theta = res_Laplace$theta, noise = res_Laplace$noise)
    print(sprintf("rmse1: %f", mean((pm1 - z)^2)))
    # predict mean of GP trained by hetGP.IMSPE
    pm2 <- apply(xgrid, 1, predict_mean, XN = res_IMSPE$X_seq, Y = res_IMSPE$Y_seq, cov.type = cov.type, nu.hat = res_IMSPE$nu.hat, theta = res_IMSPE$theta, noise = res_IMSPE$noise)
    print(sprintf("rmse2: %f", mean((pm2 - z)^2)))

    plot(res_Laplace$X_seq[, 1], res_Laplace$X_seq[, 2]) + points(XN[, 1], XN[, 2], col = "red")
    plot(res_IMSPE$X_seq[, 1], res_IMSPE$X_seq[, 2]) + points(XN[, 1], XN[, 2], col = "red")

    return(list(z = z, pm1 = pm1, pm2 = pm2))
}

d <- 2
multi_start <- 20
N <- 20
cov.type <- "Gaussian"

# initial input X
x1 <- runif(N, 0, 1)
x2 <- runif(N, 0, 1)
X_runif <- cbind(x1, x2)

x1 <- rnorm(N, 0.5, 0.1)
x2 <- rnorm(N, 0.5, 0.1)
X_norm <- cbind(x1, x2)

x1 <- rbeta(N, 1.5, 4)
x2 <- rbeta(N, 1.5, 4)
X_beta <- cbind(x1, x2)

# multistart points
Xstart <- maximinSA_LHS(lhsDesign(multi_start, 2, seed = 2024)$design)$design
# plot(Xstart)
# sample path of real function
x1 <- runif(5000, 0, 1)
x2 <- runif(5000, 0, 1)
xgrid <- cbind(x1, x2)

X <- X_runif
dist <- "unif"
step <- 50
f <- gaussian_pdf_func

par(mfrow = c(2, 1))

res <- run(X, f, dist = dist, step = step)

open3d()
plot3d(x1, x2, res$z, zlab = "real function")
open3d()
plot3d(x1, x2, res$pm1, zlab = "Laplace")
open3d()
plot3d(x1, x2, res$pm2, zlab = "IMSPE")

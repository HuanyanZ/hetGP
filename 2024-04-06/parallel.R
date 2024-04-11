library(parallel)
library(hetGP)


h <- function(x) {
    # return(1) # uniform(0, 1)
    # return(exp(-x^2 / 2)) # norm(0.5, 1)
    return((x + 10^(-9))^(-0.5) * (1 - x + 10^(-9))^(-0.5)) # beta(0.5, 0.5)
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

local_opt_fun <- function(i, Xstart, XN) {
    # print(i)
    out <- try(
        optim(
            Xstart[i, , drop = FALSE],
            Laplace_method_mv,
            XN = XN,
            method = "L-BFGS-B",
            lower = rep(0, 2), upper = rep(1, 2)
        )
    )
    if (is(out, "try-error")) {
        return(NULL)
    }
    return(out)
}

x1 <- rbeta(N, 0.5, 0.5)
x2 <- rbeta(N, 0.5, 0.5)
XN <- cbind(x1, x2)

multi_start <- 20
Xstart <- maximinSA_LHS(lhsDesign(multi_start, 2)$design)$design
print(Xstart)

ncores <- 1

start <- Sys.time()
all_res <- mclapply(1:nrow(Xstart), local_opt_fun, Xstart = Xstart, XN = XN, mc.cores = ncores)
res_min <- which.min(Reduce(c, lapply(all_res, function(x) x$value)))
opt <- all_res[[res_min]]$par
print(sprintf("optim time cost: %f", Sys.time() - start))
print(opt)


out <- multistart(
    Xstart,
    Laplace_method_mv,
    method = "L-BFGS-B",
    lower = rep(0, d), upper = rep(1, d),
    XN = XN
)
print(c(out$p1[which.min(out$value)], out$p2[which.min(out$value)]))


x1 <- runif(5000, 0, 1)
x2 <- runif(5000, 0, 1)
xgrid <- cbind(x1, x2)
# xgrid <- maximinSA_LHS(lhsDesign(1000, 2)$design)$design

z <- apply(xgrid, 1, Laplace_method_mv, XN = XN)
print(xgrid[which.min(z), ])

open3d()
plot3d(xgrid[, 1], xgrid[, 2], z)

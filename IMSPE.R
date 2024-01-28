set.seed(1234)

library("parallel")
library("hetGP")
library("MASS")
library("lhs")

# # negative log-likelihood
# nLL <- function(par, X, Y) {
#     theta <- par[1:ncol(X)]
#     tau2 <- par[ncol(X) + 1]
#     n <- length(Y)
#     K <- cov_gen(X1 = X, theta = theta) + diag(tau2, n)
#     Ki <- solve(K)
#     ldetK <- determinant(K, logarithm = TRUE)$modulus
#     (n / 2) * log(t(Y) %*% Ki %*% Y) + (1 / 2) * ldetK
# }

# # gradient of negative log-likelihood
# gnLL <- function(par, X, Y) {
#     n <- length(Y)
#     theta <- par[1:ncol(X)]
#     tau2 <- par[ncol(X) + 1]
#     K <- cov_gen(X1 = X, theta = theta) + diag(tau2, n)
#     Ki <- solve(K)
#     KiY <- Ki %*% Y
#     dlltheta <- rep(0, length(theta))
#     for (k in 1:length(dlltheta)) {
#         dotK <- K * as.matrix(dist(X[, k]))^2 / (theta[k]^2)
#         dlltheta[k] <- n * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) - sum(diag(Ki %*% dotK))
#     }
#     dlltau2 <- n * t(KiY) %*% KiY / (t(Y) %*% KiY) - sum(diag(Ki))

#     return(-c(dlltheta / 2, dlltau2 / 2))
# }


# X <- 6 * randomLHS(40, 2) - 2
# # print(X)
# X <- rbind(X, X)
# y <- X[, 1] * exp(-X[, 1]^2 - X[, 2]^2) + rnorm(nrow(X), sd = 0.01)

# Lwr <- sqrt(.Machine$double.eps)
# Upr <- 10

# out <- optim(c(rep(0.1, 2), 0.1 * var(y)), nLL, gnLL,
#     method = "L-BFGS-B",
#     lower = Lwr, upper = c(rep(Upr, 2), var(y)), X = X, Y = y
# )

# # estimator of theta and tau2
# print(out$par)

# # another way of estimating parameters
# fit <- mleHomGP(X, y, rep(Lwr, 2), rep(Upr, 2),
#     known = list(beta0 = 0),
#     init = c(list(theta = rep(0.1, 2), g = 0.1 * var(y)))
# )
# print(c(fit$theta, fit$g))

# # estimator of covariance matrix and v
# Ki <- solve(cov_gen(X, theta = out$par[1:2]) + diag(out$par[3], nrow(X)))
# nuhat <- drop(t(y) %*% Ki %*% y / nrow(X))

# # generate test data
# xx <- seq(-2, 4, length = 40)
# XX <- as.matrix(expand.grid(xx, xx))

# # posterior distribution
# KXX <- cov_gen(XX, theta = out$par[1:2]) + diag(out$par[3], nrow(XX))
# KX <- cov_gen(XX, X, theta = out$par[1:2])
# mup <- KX %*% Ki %*% y
# Sigmap <- nuhat * (KXX - KX %*% Ki %*% t(KX))

# # print(dim(XX))
# # print(dim(mup))
# # print(dim(Sigmap))

# # homoskedastic VS heteroskedastic
# hom <- mleHomGP(mcycle$times, mcycle$accel, covtype = "Matern5_2")
# het <- mleHetGP(mcycle$times, mcycle$accel, covtype = "Matern5_2")
# # print(het)
# # print(names(het))
# # print(het$theta) # covariance function parameter
# # print(het$Delta) # log noise GP parameter
# # print(het$nu_hat) # mu hat N
# # print(het$beta0) # estimated constant trend value

# Xgrid <- matrix(seq(0, 60, length = 301), ncol = 1)
# p <- predict(x = Xgrid, object = hom)
# p2 <- predict(x = Xgrid, object = het)

# # plot(Xgrid, p$mean - 2 * sqrt(p$sd2 + p$nugs), col = "red")
# # lines(Xgrid, p$mean + 2 * sqrt(p$sd2 + p$nugs), col = "red")

# # lines(Xgrid, p2$mean - 2 * sqrt(p2$sd2 + p2$nugs), col = "blue")
# # lines(Xgrid, p2$mean + 2 * sqrt(p2$sd2 + p2$nugs), col = "blue")



# Xbar <- randomLHS(200, 2)
# a <- sample(1:100, nrow(Xbar), replace = TRUE)
# X <- matrix(0, ncol = 2, nrow = sum(a))

# nf <- 0
# for (i in 1:nrow(Xbar)) {
#     X[(nf + 1):(nf + a[i]), ] <- matrix(rep(Xbar[i, ], a[i]),
#         ncol = 2,
#         byrow = TRUE
#     )
#     nf <- nf + a[i]
# }
# print(dim(X))
# Y <- apply(X, 1, sirEval)
# print(dim(Y))
# print(length(Y))
# fit <- mleHetGP(X, Y,
#     covtype = "Matern5_2", lower = rep(0.05, 2),
#     upper = rep(10, 2), settings = list(linkThetas = "none"), maxit = 1e4
# )
# # fit <- mleHetGP(X, Y,
# #     covtype = "Matern5_2", lower = rep(0.05, 2),
# #     upper = rep(10, 2), maxit = 1e4
# # )
# print(fit)
# print(fit$time) # 597.56



# data("bfs", package = "hetGP")
# thetas <- matrix(bfs.exp$theta, ncol = 1)
# bfs <- as.matrix(t(bfs.exp[, -1]))

# bfs1 <- mleHetTP(
#     X = list(
#         X0 = log10(thetas), Z0 = colMeans(log(bfs)),
#         mult = rep(nrow(bfs), ncol(bfs))
#     ), Z = log(as.numeric(bfs)),
#     lower = 10^(-4), upper = 5, covtype = "Matern5_2"
# )
# print(bfs1)


# # sequential design
# rn <- c(4.5, 5.5, 6.5, 6, 3.5)
# X0 <- matrix(seq(0.05, 0.95, length.out = length(rn)))

# X1 <- matrix(c(X0, 0.2, 0.4))
# Y1 <- c(rn, 5.2, 6.3)
# r1 <- splinefun(x = X1, y = Y1, method = "natural")

# X2 <- matrix(c(X0, 0.0, 0.3))
# Y2 <- c(rn, 7, 4)
# r2 <- splinefun(x = X2, y = Y2, method = "natural")

# XX <- matrix(seq(0, 1, by = 0.005))

# IMSPE.r <- function(x, X0, theta, r) {
#     x <- matrix(x, nrow = 1)
#     Wijs <- Wij(mu1 = rbind(X0, x), theta = theta, type = "Gaussian")
#     K <- cov_gen(X1 = rbind(X0, x), theta = theta)

#     K <- K + diag(apply(rbind(X0, x), 1, r))
#     return(1 - sum(solve(K) * Wijs))
# }

# imspe1 <- apply(XX, 1, IMSPE.r, X0 = X0, theta = 0.25, r = r1)
# imspe2 <- apply(XX, 1, IMSPE.r, X0 = X0, theta = 0.25, r = r2)
# xstar1 <- which.min(imspe1)
# xstar2 <- which.min(imspe2)

# print(c(xstar1, XX[xstar1]))
# print(c(xstar2, XX[xstar2]))


# rx <- function(x, X0, rn, theta, Ki, kstar, Wijs) {
#     x <- matrix(x, nrow = 1)
#     kn1 <- cov_gen(x, X0, theta = theta)
#     wn <- Wij(mu1 = x, mu2 = X0, theta = theta, type = "Gaussian")
#     a <- kn1 %*% Ki %*% Wijs %*% Ki %*% t(kn1) - 2 * wn %*% Ki %*% t(kn1)
#     a <- a + Wij(mu1 = x, theta = theta, type = "Gaussian")
#     Bk <- tcrossprod(Ki[, kstar], Ki[kstar, ]) / (2 / rn[kstar] - Ki[kstar, kstar])
#     b <- sum(Bk * Wijs)
#     sn <- 1 - kn1 %*% Ki %*% t(kn1)
#     return(a / b - sn)
# }

# bestk <- which.min(apply(X0, 1, IMSPE.r, X0 = X0, theta = 0.25, r = r1))
# Wijs <- Wij(X0, theta = 0.25, type = "Gaussian")
# Ki <- solve(cov_gen(X0, theta = 0.25, type = "Gaussian") + diag(rn))
# rx.thresh <- apply(XX, 1, rx,
#     X0 = X0, rn = rn, theta = 0.25, Ki = Ki,
#     kstar = bestk, Wijs = Wijs
# )
# print(rx.thresh)




# fn <- function(x) {
#     1 / 3 * (exp(sin(2 * pi * x)))
# }
# fr <- function(x) {
#     f1d2(x) + rnorm(length(x), sd = fn(x))
# }

# X <- seq(0, 1, length = 10)
# Y <- fr(X)
# mod <- mleHetGP(X = X, Z = Y, lower = 0.0001, upper = 1)

# opt <- IMSPE_optim(mod, h = 5)

# X <- c(X, opt$par)
# Ynew <- fr(opt$par)
# Y <- c(Y, Ynew)
# mod <- update(mod, Xnew = opt$par, Znew = Ynew, ginit = mod$g * 1.01)
# opt <- IMSPE_optim(mod, h = 5)

# print(Sys.time())
# for (i in 1:489) {
#     opt <- IMSPE_optim(mod, h = 5)

#     X <- c(X, opt$par)
#     Ynew <- fr(opt$par)
#     Y <- c(Y, Ynew)
#     mod <- update(mod, Xnew = opt$par, Znew = Ynew, ginit = mod$g * 1.01)
#     if (i %% 25 == 0) {
#         mod2 <- mleHetGP(
#             X = list(X0 = mod$X0, Z0 = mod$Z0, mult = mod$mult),
#             Z = mod$Z, lower = 0.0001, upper = 1
#         )
#         if (mod2$ll > mod$ll) mod <- mod2
#     }
# }
# print(Sys.time())

# xgrid <- seq(0, 1, length = 1000)
# p <- predict(mod, matrix(xgrid, ncol = 1))
# pvar <- p$sd2 + p$nugs

# plot(xgrid, p$mean, col = "red")
# lines(xgrid, f1d2(xgrid), col = "blue")


# X <- seq(0, 1, length = 10)
# Y <- fr(X)
# mod.a <- mleHetGP(X = X, Z = Y, lower = 0.0001, upper = 1)
# h <- rep(0, 500)

# print(Sys.time())
# for (i in 1:490) {
#     h[i] <- horizon(mod.a)
#     opt <- IMSPE_optim(mod.a, h = h[i])
#     X <- c(X, opt$par)
#     Ynew <- fr(opt$par)
#     Y <- c(Y, Ynew)
#     mod.a <- update(mod.a,
#         Xnew = opt$par, Znew = Ynew,
#         ginit = mod.a$g * 1.01
#     )
#     if (i %% 25 == 0) {
#         mod2 <- mleHetGP(X = list(
#             X0 = mod.a$X0, Z0 = mod.a$Z0,
#             mult = mod.a$mult
#         ), Z = mod.a$Z, lower = 0.0001, upper = 1)
#         if (mod2$ll > mod.a$ll) mod.a <- mod2
#     }
# }
# print(Sys.time())

# xgrid <- seq(0, 1, length = 1000)
# p <- predict(mod.a, matrix(xgrid, ncol = 1))
# pvar <- p$sd2 + p$nugs

# plot(xgrid, p$mean, col = "red")
# lines(xgrid, f1d2(xgrid), col = "blue")




# X <- seq(0, 1, length = 10)
# X <- c(X, X, X)
# Y <- fr(X)
# mod <- mleHetGP(X = X, Z = Y)

# ncores <- detectCores()
# print(ncores)
# print(Sys.time())
# for (i in 1:500) {
#     opt <- crit_optim(mod, crit = "crit_tMSE", h = 5)
#     # opt <- crit_optim(mod, crit = "crit_EI")
#     X <- c(X, opt$par)
#     Ynew <- -fr(opt$par)
#     Y <- c(Y, Ynew)
#     mod <- update(mod, Xnew = opt$par, Znew = Ynew, ginit = mod$g * 1.01)
#     if (i %% 25 == 0) {
#         mod2 <- mleHetGP(
#             X = list(X0 = mod$X0, Z0 = mod$Z0, mult = mod$mult),
#             Z = mod$Z, lower = 0.0001, upper = 1
#         )
#         if (mod2$ll > mod$ll) mod <- mod2
#     }
# }
# print(Sys.time())

# xgrid <- seq(0, 1, length = 1000)
# p <- predict(mod, matrix(xgrid, ncol = 1))
# pvar <- p$sd2 + p$nugs

# plot(xgrid, p$mean, col = "red")
# lines(xgrid, f1d2(xgrid), col = "blue")

# crit_tMSE(x, model, thres = 0, preds = NULL, seps = 0.05)

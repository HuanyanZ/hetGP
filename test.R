# # set.seed(125)

# library("parallel")
# library("hetGP")
# library("MASS")
# library("lhs")
# library("Metrics")
# library("DEoptim")

# cov_type <- "Matern5_2"

# noise.var <- function(x) {
#     return(1 / 3 * (exp(sin(2 * pi * sum(x)))))
# }

# gp.mean <- function(x) {
#     return(f1d2(sum(x)))
# }

# observations <- function(x) {
#     gp.mean(x) + rnorm(1, sd = noise.var(x))
# }

# nLL <- function(par, mod) {
#     XN <- mod$X0
#     nu.hat <- mod$nu_hat
#     theta <- mod$theta
#     Lambda <- mod$Lambda
#     mult <- mod$mult

#     x0 <- matrix(par, nrow = 1, byrow = TRUE)
#     # K0 <- nu.hat * cov_gen(X1 = x0, theta = theta, type = cov_type)
#     KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type)

#     K <- nu.hat * (cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(Lambda) %*% diag(1 / mult))

#     L <- chol(K)
#     v <- solve(L, t(KX))

#     mse <- -t(v) %*% v
#     return(mse)
# }

# MSE_optim <- function(mod, h) {
#     out <- DEoptim(nLL, lower = c(0), upper = c(1), DEoptim.control(trace = FALSE, itermax = 50), mod = mod)
#     return(list(par = out$optim$bestmem))
# }

# random_optim <- function(mod, h) {
#     return(list(par = runif(1, 0, 1)))
# }

# train <- function(X, Y, mod, h, obj) {
#     # print(X)
#     # print(Y)
#     # print(mod)
#     print(Sys.time())
#     for (i in 1:100) {
#         if (obj == "IMSPE") {
#             opt <- IMSPE_optim(mod, h = h)
#         }
#         if (obj == "EI") {
#             opt <- crit_optim(mod, crit = "crit_EI", h = h)
#         }
#         if (obj == "tMSE") {
#             opt <- crit_optim(mod, crit = "crit_tMSE", h = h)
#         }
#         if (obj == "MSE") {
#             opt <- MSE_optim(mod, h = h)
#         }
#         if (obj == "random") {
#             opt <- random_optim(mod, h = h)
#         }
#         # print(opt$par)
#         Xnew <- matrix(opt$par, nrow = 1)
#         X <- c(X, Xnew)
#         Ynew <- apply(Xnew, 1, observations)
#         # print(Xnew)
#         # print(Ynew)
#         Y <- c(Y, Ynew)
#         mod <- update(mod, Xnew = Xnew, Znew = Ynew, ginit = mod$g * 1.01)
#         if (i %% 25 == 0) {
#             mod2 <- mleHetGP(
#                 X = list(X0 = mod$X0, Z0 = mod$Z0, mult = mod$mult),
#                 Z = mod$Z,
#                 covtype = cov_type,
#                 settings = list(checkHom = FALSE),
#             )
#             if (mod2$ll > mod$ll) mod <- mod2
#         }
#     }
#     # print(mod)
#     print(Sys.time())
#     return(mod)
# }

# test <- function(x_test, fit, color) {
#     p <- predict(fit, x_test)
#     # print(p)
#     print(rmse(p$mean, apply(x_test, 1, gp.mean)))
#     lines(x_test, p$mean, col = color)
# }

# N <- 10

# X_seq <- seq(0, 1, length = N)

# X_unif <- sort(runif(N, 0, 1))

# X_norm <- sort(rnorm(N))
# X_norm <- (X_norm - min(X_norm)) / (max(X_norm) - min(X_norm))
# X_norm <- matrix(X_norm, nrow = N)

# X_chisq <- sort(rchisq(N, 10))
# X_chisq <- (X_chisq - min(X_chisq)) / (max(X_chisq) - min(X_chisq))
# X_chisq <- matrix(X_chisq, nrow = N)

# X <- X_norm
# X <- rbind(X, X)

# Y <- apply(X, 1, observations)
# # Y <- observations(X)

# mod <- mleHetGP(X = X, Z = Y, covtype = cov_type, settings = list(checkHom = FALSE))
# # xgrid <- matrix(seq(0, 1, length.out = 1000), ncol = 1)

# xgrid <- sort(rnorm(1000))
# xgrid <- (xgrid - min(xgrid)) / (max(xgrid) - min(xgrid))
# xgrid <- matrix(xgrid, nrow = length(xgrid))

# # xgrid <- sort(rchisq(1000, 10))
# # xgrid <- (xgrid - min(xgrid)) / (max(xgrid) - min(xgrid))
# # xgrid <- matrix(xgrid, nrow = 1000)

# plot(xgrid, apply(xgrid, 1, gp.mean), col = "black")

# fit.tMSE <- train(X, Y, mod, 1, "tMSE")
# print(fit.tMSE)
# test(xgrid, fit.tMSE, "green")

# fit.IMSPE <- train(X, Y, mod, 1, "IMSPE")
# print(fit.IMSPE)
# test(xgrid, fit.IMSPE, "blue")

# # fit.EI <- train(X, Y, mod, 1, "EI")
# # test(xgrid, fit.EI, "yellow")

# # fit.MSE <- train(X, Y, mod, 1, "MSE")
# # test(xgrid, fit.MSE, "red")

# # fit.random <- train(X, Y, mod, 1, "random")
# # test(xgrid, fit.random, "green")

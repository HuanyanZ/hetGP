set.seed(128)
# set.seed(2024)
# set.seed(1999)

library("hetGP")
library("MASS")
library("Metrics")
library("DEoptim")

# mse(x) conditioned on XN
mse <- function(par, mod, cov_type) {
    XN <- mod$X0
    nu.hat <- mod$nu_hat
    theta <- mod$theta
    Lambda <- mod$Lambda
    mult <- mod$mult

    x0 <- matrix(par, nrow = 1, byrow = TRUE)
    KX <- nu.hat * cov_gen(X1 = x0, X2 = XN, theta = theta, type = cov_type)

    K <- nu.hat * (cov_gen(X1 = XN, theta = theta, type = cov_type) + diag(Lambda) %*% diag(1 / mult))

    L <- chol(K)
    v <- solve(L, t(KX))

    mse <- -t(v) %*% v
    return(mse)
}

# find the minimum of mse
MSE_optim <- function(mod, h) {
    out <- DEoptim(mse, lower = c(0), upper = c(1), DEoptim.control(trace = FALSE, itermax = 50), mod = mod, cov_type = "Matern5_2")
    return(list(par = out$optim$bestmem))
}

data <- generate_data("unif")
X_input <- data$X_input
Y <- apply(X_input, 1, observations)
mod <- mleHetGP(X = X_input, Z = Y, covtype = "Matern5_2", settings = list(checkHom = FALSE))

xgrid <- seq(0, 1, length.out = 1000)

y <- rep(0, 1000)
for (i in 1:1000) {
    y[i] <- mse(xgrid[i], mod, cov_type = "Matern5_2")
}
png("MSE_optim.png")
plot(xgrid, y, ylim = c(-4000, 5), col = "black")

opt <- MSE_optim(mod, 1)
print(opt$par)


Xnew <- matrix(opt$par, nrow = 1)
Ynew <- apply(Xnew, 1, observations)
mod <- update(mod, Xnew = Xnew, Znew = Ynew, ginit = mod$g * 1.01)

y <- rep(0, 1000)
for (i in 1:1000) {
    y[i] <- mse(xgrid[i], mod, cov_type = "Matern5_2")
}
lines(xgrid, y, col = "red")

opt <- MSE_optim(mod, 1)
print(opt$par)


Xnew <- matrix(opt$par, nrow = 1)
Ynew <- apply(Xnew, 1, observations)
mod <- update(mod, Xnew = Xnew, Znew = Ynew, ginit = mod$g * 1.01)

y <- rep(0, 1000)
for (i in 1:1000) {
    y[i] <- mse(xgrid[i], mod, cov_type = "Matern5_2")
}
lines(xgrid, y, col = "blue")

opt <- MSE_optim(mod, 1)
print(opt$par)

dev.off()

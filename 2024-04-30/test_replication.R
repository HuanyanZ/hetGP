library(hetGP)


mse <- function(x, Xn, mult, cov.type = "Gaussian", nu.hat = 1, theta = 0.5, noise = 0.01) {
    Knx <- nu.hat * cov_gen(X1 = matrix(x, nrow = 1), X2 = Xn, theta = theta, type = cov.type) # k(x, Xn)
    S <- nu.hat * cov_gen(X1 = Xn, theta = theta, type = cov.type) + diag(rep(noise, nrow(Xn)) / mult) # Sigma_n

    L <- t(chol(S))
    v <- solve(L, t(Knx))

    mse <- nu.hat - t(v) %*% v
    return(mse)
}

delta_func <- function(x, Xn, mult, k, cov.type = "Gaussian", nu.hat = 1, theta = 0.5, noise = 0.01) {
    Knx <- nu.hat * cov_gen(X1 = matrix(x, nrow = 1), X2 = Xn, theta = theta, type = cov.type) # k(x, Xn)
    S <- nu.hat * cov_gen(X1 = Xn, theta = theta, type = cov.type) + diag(rep(noise, nrow(Xn)) / mult) # Sigma_n

    L <- t(chol(S)) # S = L L^T
    L.inv <- solve(L)

    S.inv <- t(L.inv) %*% L.inv
    uk <- S.inv[, k]

    return((uk %*% t(Knx))^2 / (mult[k] * (mult[k] + 1) / noise - S.inv[k, k]))
}


cov.type <- "Gaussian"
nu.hat <- 1
theta <- 0.5
noise <- 0.01

n <- 5
Xn <- matrix(seq(0, 1, length.out = n), nrow = n)

mult <- sample(c(1:5), n, replace = TRUE)
N <- sum(mult)

x <- runif(1)
k <- sample(c(1:n), 1)
print(sprintf("k: %d, xk: %f", k, Xn[k]))

mse.n <- mse(x, Xn, mult, cov.type = cov.type, nu.hat = nu.hat, theta = theta, noise = noise)

print(delta_func(x, Xn, mult, k, cov.type = cov.type, nu.hat = nu.hat, theta = theta, noise = noise))

mult[k] <- mult[k] + 1
mse.n1 <- mse(x, Xn, mult, cov.type = cov.type, nu.hat = nu.hat, theta = theta, noise = noise)

print(mse.n - mse.n1)


############## plot of delta_func
xgrid <- matrix(seq(0, 1, length.out = 1000), nrow = 1000)
ygrid <- apply(xgrid, 1, delta_func, Xn = Xn, mult = mult, k = k, cov.type = cov.type, nu.hat = nu.hat, theta = theta, noise = noise)

plot(xgrid, ygrid) +
    abline(v = Xn[k], col = "blue")

library(hetGP)

set.seed(123)

cov_type <- "Gaussian"
N <- 10
XN <- sort(rbeta(N, 0.5, 0.5))
XN <- matrix(XN, nrow = N)

Y <- apply(XN, 1, function(x) sin(2 * pi * x))
mod <- mleHomGP(X = XN, Z = Y, covtype = cov_type)
theta <- mod$theta
nu.hat <- mod$nu_hat
noise <- mod$g
print(sprintf("theta: %f, nu.hat: %f, noise: %f", theta, nu.hat, noise))

step <- 10
for (i in 1:step) {
    out <- IMSPE_optim(mod, h = 0)
    Xnew <- matrix(out$par, nrow = 1)
    Ynew <- apply(Xnew, 1, function(x) sin(2 * pi * x))
    # print(Xnew)

    XN <- rbind(XN, Xnew)
    Y <- c(Y, Ynew)

    mod <- update(mod, Xnew = Xnew, Znew = Ynew)
    print(sprintf("iteration: %d, Xnew: %f, theta: %f, nu_hat: %f, noise: %f", i, Xnew, mod$theta, mod$nu_hat, mod$g))

    mod <- mleHomGP(X = XN, Z = Y, covtype = cov_type)
    print(sprintf("iteration: %d, Xnew: %f, theta: %f, nu_hat: %f, noise: %f", i, Xnew, mod$theta, mod$nu_hat, mod$g))
}

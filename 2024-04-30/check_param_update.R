library(hetGP)
library("psych")

# set.seed(123)


cov_type <- "Gaussian"
N <- 100
XN <- seq(0, 1, length.out = N)
XN <- matrix(XN, nrow = N)

Y <- apply(XN, 1, function(x) sin(2 * pi * x))
mod <- mleHomGP(
    X = XN, Z = Y, covtype = cov_type,
    lower = c(1e-4, 1e-4), upper = c(1, 0.5),
    known = list(beta0 = 0),
    init = c(list(theta = c(0.1), g = 0.01))
)

step <- 50
update_cost <- 0
mleHomGP_cost <- 0

for (j in 1:step) {
    out <- IMSPE_optim(mod, h = 0)
    Xnew <- matrix(out$par, nrow = 1)
    Ynew <- apply(Xnew, 1, function(x) sin(2 * pi * x))
    # print(Xnew)

    XN <- rbind(XN, Xnew)
    Y <- c(Y, Ynew)

    reps <- find_reps(XN, Y)
    object <- mod

    t1 <- Sys.time()

    # https://github.com/cran/hetGP/blob/master/R/update_hetGP.R
    mod <- update(
        mod,
        Xnew = Xnew, Znew = Ynew,
        lower = c(1e-4, 1e-4), upper = c(1, 0.5),
        known = list(beta0 = 0)
    )
    print(sprintf("iteration: %d, [hetGP update], theta: %f, nu_hat: %f, noise: %f", j, mod$theta, mod$nu_hat, mod$g))

    t2 <- Sys.time()
    update_cost <- update_cost + t2 - t1

    t3 <- Sys.time()
    # https://github.com/cran/hetGP/blob/master/R/hetGP.R
    mod1 <- mleHomGP(
        X = XN, Z = Y, covtype = "Gaussian",
        lower = c(1e-4, 1e-4), upper = c(1, 0.5),
        known = list(beta0 = 0),
        # init = c(list(theta = c(0.1), g = 0.01))
        init = c(list(theta = c(object$theta), g = object$g))
    )
    print(sprintf("iteration: %d, [hetGP mleHomGP], theta: %f, nu_hat: %f, noise: %f", j, mod1$theta, mod1$nu_hat, mod1$g))

    t4 <- Sys.time()
    mleHomGP_cost <- mleHomGP_cost + t4 - t3
}

print(sprintf("update_cost: %f, mleHomGP_cost: %f", update_cost, mleHomGP_cost))

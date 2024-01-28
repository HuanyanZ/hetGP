set.seed(128)
# set.seed(2024)
# set.seed(1999)

library("hetGP")
library("MASS")
library("Metrics")
library("DEoptim")


generate_data <- function(type) {
    if (type == "seq") {
        X_input <- seq(0, 1, length = N)
        X_test <- seq(0, 1, length = 1000)
    }
    if (type == "unif") {
        X_input <- sort(runif(N, 0, 1))
        X_test <- sort(runif(1000, 0, 1))
    }
    if (type == "norm") {
        X_input <- sort(rnorm(N))
        X_input <- (X_input - min(X_input)) / (max(X_input) - min(X_input))
        X_test <- sort(rnorm(1000))
        X_test <- (X_test - min(X_test)) / (max(X_test) - min(X_test))
    }
    if (type == "chisq") {
        X_input <- sort(rchisq(N, 100))
        X_input <- (X_input - min(X_input)) / (max(X_input) - min(X_input))
        X_test <- sort(rchisq(1000, 100))
        X_test <- (X_test - min(X_test)) / (max(X_test) - min(X_test))
    }

    X_input <- matrix(X_input, nrow = N)
    # X_input <- rbind(X_input, X_input)

    X_test <- matrix(X_test, nrow = 1000)

    return(list(X_input = X_input, X_test = X_test))
}


# real noise variance function
noise.var <- function(x) {
    return(1 / 3 * (exp(sin(2 * pi * sum(x)))))
}

generating_func <- function(x) {
    return(f1d2(sum(x)))
}

# observations(x) = y(x) + noise
observations <- function(x) {
    generating_func(x) + rnorm(1, sd = noise.var(x))
}

train <- function(X, step, h, obj) {
    Y <- apply(X, 1, observations)
    mod <- mleHetGP(X = X, Z = Y, covtype = cov_type, settings = list(checkHom = FALSE))

    for (i in 1:step) {
        if (obj == "IMSPE") {
            opt <- IMSPE_optim(mod, h = h)
        }
        if (obj == "EI") {
            opt <- crit_optim(mod, crit = "crit_EI", h = h)
        }
        if (obj == "tMSE") {
            opt <- crit_optim(mod, crit = "crit_tMSE", h = h)
        }
        if (obj == "MSE") {
            opt <- MSE_optim(mod, h = h)
        }
        Xnew <- matrix(opt$par, nrow = 1)
        X <- c(X, Xnew)
        Ynew <- apply(Xnew, 1, observations)
        Y <- c(Y, Ynew)
        mod <- update(mod, Xnew = Xnew, Znew = Ynew, ginit = mod$g * 1.01)
        if (i %% 25 == 0) {
            mod2 <- mleHetGP(
                X = list(X0 = mod$X0, Z0 = mod$Z0, mult = mod$mult),
                Z = mod$Z,
                covtype = cov_type,
                settings = list(checkHom = FALSE),
            )
            if (mod2$ll > mod$ll) mod <- mod2
        }
    }
    return(mod)
}

test <- function(x_test, fit, color) {
    p <- predict(fit, x_test)
    evaluator <- rmse(p$mean, apply(x_test, 1, generating_func))

    png(paste(paste(obj, dist, cov_type, "no_repli", sep = "_"), ".png", sep = ""))
    plot(
        x_test, apply(x_test, 1, generating_func),
        xlab = paste("rmse", as.character(round(evaluator, 6)), sep = "="), ylab = "Y",
        col = "black"
    )
    lines(x_test, p$mean, col = color)
    dev.off()
}

# initial input size
N <- 10
step <- 300

dist <- "norm"
obj <- "IMSPE"

# type of covariance function, Matern5_2, Matern3_2, and square exponential
cov_type <- "Matern5_2"

data <- generate_data(dist)
X_input <- data$X_input
X_test <- data$X_test

fit <- train(X_input, step, 1, obj)
# print(fit)

test(X_test, fit, "blue")


# fit.tMSE <- train(X, step, 1, "tMSE")
# print(fit.tMSE)
# test(x_test, fit.tMSE, "green")

# fit.MSE <- train(X, step, 1, "MSE")
# print(fit.MSE)
# test(x_test, fit.MSE, "red")

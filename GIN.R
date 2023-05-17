library(BAS)
N <- 1e4
a <- 1
d <- 4
y <- rnorm(N, mean = a*d, sd = d)

asd <- function(x, y) {-mean(dnorm(y, mean = x[1]*x[2], sd = x[2], log = TRUE))}
mles <- optim(c(0, 1), fn = asd, y = y)
mles$par

# Generalized inverse normal density ------------------------------------------
Dgin1 <- function(z, v = 2, mu = 0) {
    if (z == 0 | v <= 1) {
        res <- 0
    } else {
        if (mu == 0) {
            hyp <- 1
        } else {
            hyp <- hypergeometric1F1(v12, .5, .5*(mu/tau)^2, log = FALSE)
        }
        v12 <- .5*(v-1)
        kern <- abs(z)^(-v)*exp(-.5/tau^2*(1/z^2 - 2*mu/z))
        cons <- (2*tau^2)^v12*gamma(v12)*hyp
        res <- kern/cons
    }
    return(res)
}
Dgin <- function(z, v = 2, mu = 0, tau = 1) {
    if (z == 0 | v <= 1 | tau <= 0) {
        res <- 0
    } else {
        if (mu == 0) {
            hyp <- 1
        } else {
            hyp <- hypergeometric1F1(v12, .5, .5*(mu/tau)^2, log = FALSE)
        }
        v12 <- .5*(v-1)
        kern <- abs(z)^(-v)*exp(-.5/tau^2*(1/z^2 - 2*mu/z))
        cons <- (2*tau^2)^v12*gamma(v12)*hyp
        res <- kern/cons
    }
    return(res)
}
dgin <- Vectorize(Dgin)

x <- seq(-5, 5, length = 200)
densx <- matrix(NA, 200, 6)
for (j in 1:200) {
    densx[j, 1] <- dgin(x[j], 1.5, 0, 1)
    densx[j, 2] <- dgin(x[j], 2, 0, 1)
    densx[j, 3] <- dgin(x[j], 4, 0, 1)
    
    densx[j, 4] <- dgin(x[j], 2, 0, 1)
    densx[j, 5] <- dgin(x[j], 2, 1, 1)
    densx[j, 6] <- dgin(x[j], 2, 2, 1)
}

par(mfrow = c(1, 2))
plot(x, densx[, 1], type = 'l', ylim = range(densx[, 1:3]))
lines(x, densx[, 2], type = 'l', lty = 2)
lines(x, densx[, 3], type = 'l', lty = 3)

plot(x, densx[, 4], type = 'l', ylim = range(densx[, 4:6]))
lines(x, densx[, 5], type = 'l', lty = 2)
lines(x, densx[, 6], type = 'l', lty = 3)


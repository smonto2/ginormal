# Choose the configuration of values you want to plot
library(ginormal)
alpha <- 5
mu <- 0
tau <- 1
n_draws <- 1000
n_grid <- 200

# Draw from the standard and truncated distributions
z_unc <- rgin(n_draws, alpha, mu, tau)
z_p <- rtgin(n_draws, alpha, mu, tau, sign = TRUE)
z_n <- rtgin(n_draws, alpha, mu, tau, sign = FALSE)

# Compare histogram for draws with density
z_vals <- seq(-5, 5, length.out = n_grid)
fz_unc <- sapply(z_vals, function(z) dgin(z, alpha, mu, tau, FALSE))
fz_p <- sapply(z_vals[z_vals > 0], function(z) dtgin(z, alpha, mu, tau, TRUE, FALSE))
fz_p <- c(rep(0, n_grid - sum(z_vals > 0)), fz_p)
fz_n <- sapply(z_vals[z_vals < 0], function(z) dtgin(z, alpha, mu, tau, FALSE, FALSE))
fz_n <- c(fz_n, rep(0, n_grid - sum(z_vals < 0)))

# Standard density on the full real line
temp <- hist(z_unc, breaks = 100, plot = FALSE)
plot(temp, freq = FALSE, xlim = c(-5, 5), ylim = range(c(fz_unc, temp$density)),
     main = '', xlab = 'Values', ylab = 'Density', col = 'blue')
lines(z_vals, fz_unc, col = 'red', lwd = 2)

# Density truncated to positive reals
temp <- hist(z_p, breaks = 100, plot = FALSE)
plot(temp, freq = FALSE, xlim = c(-5, 5), ylim = range(c(fz_p, temp$density)),
     main = '', xlab = 'Values', ylab = 'Density', col = 'blue')
lines(z_vals, fz_p, col = 'red', lwd = 2)

# Density truncated to negative reals
temp <- hist(z_n, breaks = 100, plot = FALSE)
plot(temp, freq = FALSE, xlim = c(-5, 5), ylim = range(c(fz_n, temp$density)),
     main = '', xlab = 'Values', ylab = 'Density', col = 'blue')
lines(z_vals, fz_n, col = 'red', lwd = 2)

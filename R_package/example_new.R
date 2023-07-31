# Try it! Let's plot the true pdf and the sampled histogram for many draws
# Set-up with suitable parameter values
n_draws <- 1000
a1 <- 2.5
a2 <- 5
m1 <- 0
m2 <- 1
m3 <- -1
t1 <- 1
t2 <- 2

# Choose the configuration of values you want to plot
n_draws <- 200
a <- a2
m <- m1
t <- t1
z_vals <- seq(-5, 5, length.out = n_draws)
fz_unc <- sapply(z_vals, function(z) dgin(z, a, m, t, FALSE))
fz_p <- sapply(z_vals[z_vals > 0], function(z) dtgin(z, a, m, t, TRUE, FALSE))
fz_p <- c(rep(0, n_draws - sum(z_vals > 0)), fz_p)
fz_n <- sapply(z_vals[z_vals < 0], function(z) dtgin(z, a, m, t, FALSE, FALSE))
fz_n <- c(fz_n, rep(0, n_draws - sum(z_vals < 0)))
z_unc <- rgin(n_draws, a, m, t)
z_p <- rtgin(n_draws, a, m, t, TRUE)
z_n <- rtgin(n_draws, a, m, t, FALSE)

# Full density on real line
hist(z_unc, breaks = 50, freq = FALSE, xlim = range(z_vals),
     col = "blue", main = "", xlab = "Values", ylab = "Density") #FALSE means density not frequency
lines(z_vals, fz_unc, col = "red") #we can add type='l'(line) or type='p'(point)

# Try it! Let's plot the true pdf and the sampled histogram for many draws
import ginormal as gin
import numpy as np
from matplotlib import pyplot as plt

# Choose the configuration of values you want to plot
alpha = 5
mu = 0
tau = 1
n_draws = 1000
n_grid = 200

# Draw from the standard and truncated distributions
z_unc = gin.rgin(n_draws, alpha, mu, tau)
z_p = gin.rtgin(n_draws, alpha, mu, tau, True)
z_n = gin.rtgin(n_draws, alpha, mu, tau, False)

# Compare histogram for draws with density
z_vals = np.linspace(-5, 5, n_grid)
fz_unc = np.array([gin.dgin(z, alpha, mu, tau, False) for z in z_vals])
fz_p = np.array([gin.dtgin(z, alpha, mu, tau, True, False) for z in z_vals[z_vals > 0]])
fz_p = np.append(np.zeros(n_grid - len(fz_p)), fz_p)
fz_n = np.array([gin.dtgin(z, alpha, mu, tau, False, False) for z in z_vals[z_vals < 0]])
fz_n = np.append(fz_n, np.zeros(n_grid - len(fz_n)))

# Standard density on the full real line
plt.figure()
plt.hist(x = z_unc, bins = 100, density = True)
plt.plot(z_vals, fz_unc, 'r')

# Density truncated to positive reals
plt.figure()
plt.hist(x = z_p, bins = 100, density = True)
plt.plot(z_vals, fz_p, 'r')

# Density truncated to negative reals
plt.figure()
plt.hist(x = z_n, bins = 100, density = True)
plt.plot(z_vals, fz_n, 'r')

plt.show()

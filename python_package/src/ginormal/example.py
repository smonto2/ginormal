# Try it! Let's plot the true pdf and the sampled histogram for many draws
import ginormal as gin
import numpy as np
from matplotlib import pyplot as plt

# Set-up with suitable parameter values
n_draws = 1000
a1 = 2.5
a2 = 5
m1 = 0
m2 = 1
m3 = -1
t1 = 1
t2 = 2

# Choose the configuration of values you want to plot
n_grid = 200
a = a2
m = m1
t = t1
z_vals = np.linspace(-5, 5, n_grid)
fz_unc = np.array([gin.dgin(z, a, m, t, False) for z in z_vals])
fz_p = np.array([gin.dtgin(z, a, m, t, True, False) for z in z_vals[z_vals > 0]])
fz_p = np.append(np.zeros(n_grid - len(fz_p)), fz_p)
fz_n = np.array([gin.dtgin(z, a, m, t, False, False) for z in z_vals[z_vals < 0]])
fz_n = np.append(fz_n, np.zeros(n_grid - len(fz_n)))
z_unc = gin.rgin(n_draws, a, m, t, 'hormann')
z_p = gin.rtgin(n_draws, a, m, t, True, 'hormann')
z_n = gin.rtgin(n_draws, a, m, t, False, 'hormann')

# Full density on real line
plt.figure()
plt.hist(x = z_unc, bins = 50, density = True)
plt.plot(z_vals, fz_unc, 'r')
plt.show()

# Truncated to positive reals
# plt.figure()
# plt.hist(x = z_p, bins = 50, density = True)
# plt.plot(z_vals, fz_p, 'r')
# plt.show()

# Truncated to negative reals
# plt.figure()
# plt.hist(x = z_n, bins = 50, density = True)
# plt.plot(z_vals, fz_n, 'r')
# plt.show()

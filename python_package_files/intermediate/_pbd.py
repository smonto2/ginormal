from scipy import special

# Truncated generalized inverse normal log-densities
def pbdam(a, m):
    am = -a + 1
    return special.pbdv(am, m)[0]
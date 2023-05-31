# Generalized Inverse Normal (GIN) distribution
The GIN package provides the density function and random variable generation from the generalized inverse normal (GIN) distribution introduced by [Robert (1991)](#1). The GIN distribution is a way to generalize the distribution of the reciprocal of a normal random variable. That is, the distribution generalizes the distribution of the random variable $1/X$ where $X \sim \text{Normal}(\mu, \sigma^2)$. This distribution is *different* from the generalized inverse Gaussian (GIG) distribution [(Jørgensen, 2012)](#2) despite the similarities in naming.

This package is the first to provide an efficient sampling algorithm for drawing from the GIN distribution. Further details of the distribution, theoretical guarantees and pseudo-code for the sampling algorithms, as well as an application to Bayesian estimation of network formation models can be found in [Ding, Estrada and Montoya-Blandón (2023)](#3).

## Density functions

The GIN distribution is supported on the entire real line $(-\infty, \infty)$ and takes three parameters:
- $\alpha > 1$, a degrees-of-freedom parameter,
- $\mu \in (-\infty, \infty)$, similar to a location parameter, it shifts the density of the distribution left and right,
- $\tau > 0$, similar to a scale parameter, it spreads the density of the distribution.

Let $Z \sim \text{GIN}(\alpha, \mu, \tau)$. The density function is given by $$f_Z(z) = \frac{1}{C(\alpha, \mu, \tau)} |z|^{-\alpha}\exp\left[-\frac{1}{2\tau^2} \left( \frac{1}{z} - \mu \right)^2 \right]$$ where the proportionality constant can be written in closed form as $$C(\nu, \gamma, \tau) = (\sqrt{2} \tau)^{\alpha-1} \exp\left(- \frac{\mu^2}{2\tau^2} \right) \Gamma\left(\frac{\alpha-1}{2}\right){}_1F_1\left(\frac{\alpha-1}{2}; \frac{1}{2}; \frac{\mu^2}{\tau^2}\right)$$ and ${}_1F_1(a, b; x)$ is the confluent hypergeometric function.

TRUNCATED DISTRIBUTION TOO

## Random variable generation


## References
1. <a id="1"> [Robert, C. (1991). Generalized inverse normal distributions. Statistics & Probability Letters, 11(1), 37-41.](https://doi.org/10.1016/0167-7152%2891%2990174-P) </a>
2. <a id="2"> [Jørgensen, B. (2012). Statistical properties of the generalized inverse Gaussian distribution (Vol. 9). Springer Science & Business Media.](https://link.springer.com/book/10.1007/978-1-4612-5698-4) </a>
3. <a id="3"> [Ding, C., Estrada, J., and Montoya-Blandón, S. (2023). Bayesian Inference of Network Formation Models with Payoff Externalities. Working Paper.](https://www.smontoyablandon.com/publication/networks/network_externalities.pdf) </a>

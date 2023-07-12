# Generalized Inverse Normal (GIN) distribution
The `ginormal` package provides the density function and random variable generation from the generalized inverse normal (GIN) distribution introduced by [Robert (1991)](#2). The GIN distribution is a way to generalize the distribution of the reciprocal of a normal random variable. That is, the distribution generalizes the distribution of the random variable $1/X$ where $X \sim \text{Normal}(\mu, \sigma^2)$. This distribution is *different* from the generalized inverse Gaussian (GIG) distribution [(Jørgensen, 2012)](#3) despite the similarities in naming (see [below](#digression)).

This package is the first to provide an efficient sampling algorithm for drawing from the GIN distribution. We provide similar routines for the GIN distribution truncated to the positive or negative reals. Further details of the distribution, theoretical guarantees and pseudo-code for the sampling algorithms, as well as an application to Bayesian estimation of network formation models can be found in [Ding, Estrada and Montoya-Blandón (2023)](#1).

## Routines

The GIN distribution is supported on the entire real line $(-\infty, \infty)$ and takes three parameters:
- $\alpha > 1$, a degrees-of-freedom parameter,
- $\mu \in (-\infty, \infty)$, similar to a location parameter, it shifts the density of the distribution left and right,
- $\tau > 0$, similar to a scale parameter, it spreads the density of the distribution.

Provided with the package are four main routines:
1. `dgin(z, alpha, mu, tau, log = TRUE, quasi = FALSE)`
2. `dtgin(z, alpha, mu, tau, sign, log = TRUE, quasi = FALSE)`
3. `rgin(size, alpha, mu, tau, algo)`
4. `rtgin(size, alpha, mu, tau, sign, algo)`

The first two compute the densities and the last two are used for random number generation. Density routines take in the quantile `z`, parameters, and two optional logical arguments:
- `log`, should the logarithm of the density be returned? Defaults to `TRUE`.
- `quasi`, should the value of the kernel (or quasi-density) be returned? Defaults to `FALSE`.

Generation routines take the same parameters but require a `size` argument determining the amount of random variates to generate. These routines only admit a parameter `alpha` larger than 2. They take an additional argument `algo`, which can be either `"hormann"` or `"leydold"`, and defaults to `"hormann"` as our prefered method. See [below for details](#rvgeneration) on both points.

Those routines including "`t`" in their name work for the truncated variants. They take an additional logical argument `sign`, where `sign = TRUE` implies truncation to positive numbers $(z > 0)$ and `sign = FALSE` to negative numbers $(z < 0)$.

## Density functions

Let $Z \sim \text{GIN}(\alpha, \mu, \tau)$. The GIN density function is given by
$$f_Z(z) = \frac{1}{C(\alpha, \mu, \tau)} |z|^{-\alpha}\exp\left[-\frac{1}{2\tau^2} \left( \frac{1}{z} - \mu \right)^2 \right] \equiv \frac{g(z; \alpha, \mu, \tau)}{C(\alpha, \mu, \tau)}$$
where $g(z; \alpha, \mu, \tau)$ is the kernel or quasi-density and the proportionality constant can be written in closed form as
```math
C(\alpha, \mu, \tau) = (\sqrt{2} \tau)^{\alpha-1} \exp\left(- \frac{\mu^2}{2\tau^2} \right) \Gamma\left(\frac{\alpha-1}{2}\right) {}_1F_1\left(\frac{\alpha-1}{2}; \frac{1}{2}; \frac{\mu^2}{\tau^2}\right)
```
where $`{}_1F_1(a, b; x)`$ is the [confluent hypergeometric function](https://mathworld.wolfram.com/ConfluentHypergeometricFunctionoftheFirstKind.html). In addition to the density and generation routines for the GIN distribution, we provide similar routines for the GIN distribution truncated to positive or negative numbers. These are denoted by $\text{GIN}^{+}$ when truncated to $(0, \infty)$ and by $\text{GIN}^{-}$ when truncated to $(-\infty, 0)$. Let $Z^{+} \sim \text{GIN}^{+}(\alpha, \mu, \tau)$ and $Z^{-} \sim \text{GIN}^{-}(\alpha, \mu, \tau)$. Their densities are given by
$$f_{Z^{+}}(z) = \frac{g(z; \alpha, \mu, \tau)}{C^{+}(\alpha, \mu, \tau)} \mathbb{I}(z > 0)$$
$$f_{Z^{-}}(z) = \frac{g(z; \alpha, \mu, \tau)}{C^{-}(\alpha, \mu, \tau)} \mathbb{I}(z < 0)$$
with proportionality constants
$$C^{-}(\alpha, \mu) = e^{-\frac{\mu^2}{4}} \Gamma(\alpha - 1) D_{-(\alpha-1)}(-\mu)$$
$$C^{+}(\alpha, \mu) = e^{-\frac{\mu^2}{4}} \Gamma(\alpha - 1) D_{-(\alpha-1)}(\mu)$$
where $\mathbb{I}(\cdot)$ is the indicator function that is 1 when its argument is true and 0 otherwise, and $D_\nu(x)$ is the [parabolic cylinder function](https://mathworld.wolfram.com/ParabolicCylinderFunction.html).

## Random variable generation

<a id="rvgeneration"> </a> [Ding, Estrada and Montoya-Blandón (2023)](#1) provide an efficient sampling algorithm for the GIN distribution and its truncated variants for the case of $\alpha > 2$. This restriction is not of concern if the goal is the perform Bayesian estimation using this distribution (see [below for more details](#digression) and Remark 2 in the paper). Generation is done using the ratio-of-uniforms method with mode shift ([Kinderman and Monahan, 1977](#4)), which requires the computation of the minimal bounding rectangle. We implement two alternatives found in the literature:
1. [Leydold (2001)](#5) that requires information on the proportionality constants.
2. [Hörmann and Leydold (2014)](#6) that requires solving a cubic equation. This is our prefered method and the default in the package.

## Digression: Difference between GIN and GIG distributions

<a id="digression"> </a> While the kernels &mdash; and therefore the sampling techniques &mdash; for the GIN and GIG distribution are similar, these two distribution share some important differences. The main is their conceptualization, as they both attempt to generalize the idea of an inverse normal distribution in different ways. The GIG distribution does so by choosing cumulants that are inverses to those of the normal distribution. The GIN distribution does so by directly using the density of the reciprocal after a change of variables. Another important difference comes from their use as conjugate priors in Bayesian analysis:
- $\theta \sim \text{GIN}(\alpha, \mu, \tau)$ is the conjugate prior if observations are random samples from $Y \sim \text{Normal}(\theta \delta, \theta^2)$
- $\theta \sim \text{GIG}(\alpha, \mu, \tau)$ is the conjugate prior if observations are random samples from $Y \sim \text{Normal}(\theta \delta, \theta)$

These are similar mixture models but carry different interpretations and thus require different posterior sampling algorithms. This interpretation also shows why the restriction of $\alpha \geq 2$ is not binding if the goal is to perform Bayesian analysis. A prior $\gamma \sim \text{GIN}(\alpha_0, \mu_0, \tau_0)$ with $\alpha_0 = 1 + \varepsilon$ is non-informative when $\varepsilon > 0$ is arbitrarily small. However, the posterior distribution will have degrees-of-freedom parameter $\alpha_N = N + 1 + \varepsilon$ where $N$ is the sample size. As $N \geq 1$ implies $\alpha_N > 2$, for a Bayesian analysis we always need to draw from the GIN distribution with $\alpha > 2$.

## References
1. <a id="1"> [Ding, C., Estrada, J., and Montoya-Blandón, S. (2023). Bayesian Inference of Network Formation Models with Payoff Externalities. Working Paper.](https://www.smontoyablandon.com/publication/networks/network_externalities.pdf) </a>
2. <a id="2"> [Robert, C. (1991). Generalized inverse normal distributions. Statistics & Probability Letters, 11(1), 37-41.](https://doi.org/10.1016/0167-7152%2891%2990174-P) </a>
3. <a id="3"> [Jørgensen, B. (2012). Statistical properties of the generalized inverse Gaussian distribution (Vol. 9). Springer Science & Business Media.](https://link.springer.com/book/10.1007/978-1-4612-5698-4) </a>
4. <a id="4"> [Kinderman, A. J., & Monahan, J. F. (1977). Computer generation of random variables using the ratio of uniform deviates. ACM Transactions on Mathematical Software (TOMS), 3(3), 257-260.](https://doi.org/10.1145/355744.355750)
5. <a id="5"> [Leydold, J. (2001). A simple universal generator for continuous and discrete univariate T-concave distributions. ACM Transactions on Mathematical Software (TOMS), 27(1), 66-82.](https://doi.org/10.1145/382043.382322) </a>
6. <a id="6"> [Hörmann, W., & Leydold, J. (2014). Generating generalized inverse Gaussian random variates. Statistics and Computing, 24, 547-557.](https://doi.org/10.1007/s11222-013-9387-3) </a>

# Time model for TEHS (dynamic number of LULC covariates)
# Expected data:
# - nobs: number of steps
# - ncov: number of LULC covariates
# - delta.time: observed travel time
# - dist: step distance
# - Miss: missing-fix indicator (0/1)
# - xmat: matrix [nobs x ncov] of LULC covariates

model{
  for (i in 1:nobs) {
    # Mean of gamma distribution
    log_mu[i] <- b0 + inprod(xmat[i, 1:ncov], b[1:ncov])
    mu[i] <- dist[i] * exp(log_mu[i])

    # Dispersion terms
    rate[i] <- exp(g0 + g1 * Miss[i])
    shape[i] <- mu[i] * rate[i]

    # Likelihood
    delta.time[i] ~ dgamma(shape[i], rate[i])
  }

  # Priors
  b0 ~ dnorm(0, 0.01)
  for (k in 1:ncov) {
    b[k] ~ dnorm(0, 1)
  }
  g0 ~ dnorm(0, 0.01)
  g1 ~ dnorm(0, 0.01)
}

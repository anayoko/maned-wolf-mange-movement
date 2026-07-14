# TEHS selection model
# Expected data:
# - nobs: number of observed steps
# - ncov: number of LULC covariates
# - xmat0..xmat4: covariate matrices for realized + 4 alternatives
# - pmov0..pmov4: movement probabilities from time model
# - y: ones vector (for one-trick)

model{
  eps <- 1.0E-12

  for (i in 1:nobs) {
    # Relative probability of each alternative
    p0[i] <- pmov0[i] * exp(inprod(xmat0[i, 1:ncov], betas[1:ncov]))
    p1[i] <- pmov1[i] * exp(inprod(xmat1[i, 1:ncov], betas[1:ncov]))
    p2[i] <- pmov2[i] * exp(inprod(xmat2[i, 1:ncov], betas[1:ncov]))
    p3[i] <- pmov3[i] * exp(inprod(xmat3[i, 1:ncov], betas[1:ncov]))
    p4[i] <- pmov4[i] * exp(inprod(xmat4[i, 1:ncov], betas[1:ncov]))

    denom[i] <- p0[i] + p1[i] + p2[i] + p3[i] + p4[i] + eps
    pi[i] <- p0[i] / denom[i]

    # One-trick likelihood
    y[i] ~ dbern(pi[i])
  }

  # Priors
  for (k in 1:ncov) {
    betas[k] ~ dnorm(0, 1)
  }
}

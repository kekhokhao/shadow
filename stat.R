# Comparing two proportions
pro = prop.test(x=c(312,408), n=c(1000,1200), alternative = "less")
# One-propotion z-test:
pro = prop.test(x = 273, n = 300, p = 0.95, alternative = "less")


# Two sample t-test
(93.47-96.15)/sqrt(10.24^2/68+14.61^2/80)
# One sample t-test
(449.4 - 450) * sqrt(120) / 3.8

## Interval estimation of population mean:
xbar = 122
SE = 4/sqrt(100)
E = qt(.975, df= 100-1)*SE
# SE is standard error
xbar + c(E,-E)

## Estimation of population proportion:
pbar = 102/400
SE = sqrt(pbar*(1-pbar)/400)
E = qnorm(.99)*SE
pbar + c(-E, E) 

## Distributions
# Poisson
# expectation = lambda, var = lambda
# Calculating density
dpois(10, lambda = 12)
# Calculating distribution:
ppois(n, lambda = ...)

# Binomial
# E = np, Var = np(1-p)
dbinom(x = , size = , prob = )
pbinom(x = , size =,  prob = , lower.tail = TRUE)

# Normal
# E = mean, Var = sd^2
dnorm(x, mean = , sd =)
pnorm(q, mean =, sd = , lower.tail = TRUE)
qnorm(p, mean = 0, sd = 1)

# Exponential
# E = 1/rate, V = 1/rate^2
dexp(x, rate = 1)
pexp(1, rate = 4, lower.tail = T)

# Covariance
Cov = E(XY) - E(X)*E(Y)
# Calclating covariance with joint probability
x <- 2:4
y <- 1:3
joint_pmf <- matrix(c(1/12, 1/6, 1/12, 1/6, 0, 1/6, 0, 1/3, 0), ncol = 3, byrow = T)
# E(X), E(Y)
mu_x <- rowSums(joint_pmf) %*% x
mu_y <- colSums(joint_pmf) %*% y
cov_xy <- x %*% joint_pmf %*% y - mu_x - mu_y
v_x <- rowSums(joint_pmf) %*% x^2 - mu_x^2
v_y <- colSums(joint_pmf) %*% y^2 - mu_y^2
# Calculating correlation
p <- cov_xy/(sqrt(v_x)*sqrt(v_y))

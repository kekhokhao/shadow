# Comparing two proportions k1/n1 and k2/n2:
# Alternative hypothesis is k1/n1 < k2/n2
pro = prop.test(x=c(k1,k2), n=c(n1,n2), alternative = "less")
# One-propotion z-test:
# Alternative hypothesis is k/n < p
pro = prop.test(x = k, n = ..., p = 0.95, alternative = "less")


# Two sample t-test
(mu1-mu2)/sqrt(s1^2/n1+s^2/n2)
# One sample t-test
(x - mu) * sqrt(n) / s

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

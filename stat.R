# Comparing two proportions
prop.test(x=c(490,400), n=c(500,500))
# One-propotion z-test:
prop.test(x, n, p)

# Generate data for t-test:
rnorm(n, mean, sd)

# Two sample t-test
(mean1-mean2)/sqrt(var1/m+var2/n)
# One sample t-test
(mean - e) * sqrt(n) / s

## Interval estimation of population mean:
xbar = mean(data)
SE = s/sqrt(n)
E = qt(.975, df= n-1)*SE
# SE is standard error
xbar + c(E,-E)

## Estimation of population proportion:
pbar = k/n
SE = sqrt(pbar*(1-pbar)/n)
E = qnorm(.975)*SE
pbar + c(-E, E) 

## Distributions
# Poisson
# expectation = lambda, var = lambda
# Calculating density
dpois(x, lambda)
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
pexp(7000, rate = 1/3300, lower.tail = FALSE)

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

################################################ Bears with genotype frequency of 0.2 (N=100) ####################################################

#################### Case 1: alpha and beta are equal
# Creating the Sequence
gfg = seq(0,1, by=0.01)

plot(gfg, dbeta(gfg, 1,1), xlab="X",
     ylab = "Beta Density", type = "l",
     col = "Red") 

#################### Case 2: probability of observing 40 alleles out of 200 observations

# Creating the Sequence
gfg = seq(0,1, by=0.01)

plot(gfg, dbeta(gfg, 41,161), xlab="X",
     ylab = "Beta Density", type = "l",
     col = "Red")

# Calculating 95 % interval
d1 <- qbeta(0.025, 41,161, ncp = 0, lower.tail = TRUE, log.p = FALSE)
d1 <- qbeta(0.975, 41,161, ncp = 0, lower.tail = TRUE, log.p = FALSE)

# Calculating notable quantiles
gfg = seq(0.025,0.975, by=0.01)
qbeta(gfg, 41,161, ncp = 0, lower.tail = TRUE, log.p = FALSE)

gfg = seq(0,1, by=0.01)

alpha <- 41
beta <- 161
MaximumAlleleProbability <- (alpha - 1) / (alpha + beta - 2)
print(paste("In a beta distribution with alpha ", alpha, "and beta ", beta, " we have a credible interval between", d1, " and ", d2, " with highest value at ", MaximumAlleleProbability))

#################### Case 3: probability of observing 2 alleles out of 20 observations

plot(gfg, dbeta(gfg, 3,19), xlab="X",
     ylab = "Beta Density", type = "l",
     col = "Red")

# Calculating 95 % interval
d1 <- qbeta(0.025, 3,19, ncp = 0, lower.tail = TRUE, log.p = FALSE)
d2 <- qbeta(0.975, 3,19, ncp = 0, lower.tail = TRUE, log.p = FALSE)

# Calculating notable quantiles
gfg = seq(0.025,0.975, by=0.01)
qbeta(gfg, 3,19, ncp = 0, lower.tail = TRUE, log.p = FALSE)

alpha <- 3
beta <- 19
MaximumAlleleProbability <- (alpha - 1) / (alpha + beta - 2)
print(paste("In a beta distribution with alpha ", alpha, "and beta ", beta, " we have a credible interval between", d1, " and ", d2, " with highest value at ", MaximumAlleleProbability))


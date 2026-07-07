# pkGP_recodeDiagonal.R
# Andrew Brown
# July 7, 2026
#
# This file is to clean up the diagonal example in the pkGP paper. Section 5.2
################################################################################

rm(list = ls())

library(orthopolynom)
library(statmod)
library(lhs)
library(plot3D)
library(viridis)
#sessionInfo()

set.seed(420)


# Define the test function
f <- function(x) {
  result <- x[2] *
    sqrt(1 + x[1]) *
    cos(pi * x[2]) *
    sin(0.5 * pi * (x[1] - x[2]) + 1) *
    exp(0.5 * (x[1] + x[2])^2)

  return(result)
}


K2 <- function(X, Y) {
  params[1] * Matern(rdist(X, Y), range = sqrt(params[2]) / sqrt(3), nu = 3 / 2)
}


################################################################################

J <- 15 # Max order of polynomial approximation for Rayleigh-Ritz method

# This will hold the sup norm relative errors and RMSE's associated with
# Legendre polyn. approximations of varying orders, given by J
leg.max <- rep(0, J - 4)
leg.rmse <- leg.max
#leg.rmse=leg.max=rep(0,J-4);

params <- c(1, 1, 1) # (scale param., length-scale in x, length-scale in y)


# Will bury this in a loop later, looping over polyn. order m
m <- 15


# Get nodes/weights for Gaussian quadrature
quad <- gauss.quad(4 * m, kind = "legendre")
# quad <- gauss.quad(4*m, kind= "chebyshev1")
n.nodes <- length(quad[[1]])


#list of legendre polynomials
leg <- lapply(legendre.polynomials(m, normalized = TRUE), as.function)
basis <- leg


# Going across the diagonal is like the boundary example, but the subspace is a
# single smooth "edge" (the main diagonal) instead of the 4 edges of the square.
# Because there are no corners, we need neither the composite (per-edge)
# quadrature nor the piecewise basis of the boundary case: one Gauss-Legendre
# rule on [-1, 1] and the plain normalized Legendre basis suffice.

# Parameterize the diagonal by s in [-1, 1]: point(s) = (s, s). Single
# n.nodes x 2 matrix (no length-one lists).
xg <- cbind(quad$nodes, quad$nodes)

# Kernel matrix over all diagonal quadrature points (n.nodes x n.nodes). Only
# one block now, so no K.bas[["ij"]] cross-edge bookkeeping.
K.bas <- K2(xg, xg)

# Evaluate the basis functions at the diagonal nodes. Normalized Legendre polys
# are already orthonormal w.r.t. ds on [-1, 1], so no rescaling is needed
# (unlike the boundary case's per-edge subintervals, which needed a factor of 2).
L <- length(basis)
phi <- matrix(0, nrow = n.nodes, ncol = L)
for (i in 1:L) {
  phi[, i] <- sapply(quad$nodes, basis[[i]])
}

# Discretized integral operator: phi^T W K.bas W phi  (L x L). K.int = A matrix
# of Oya et al. Building it from the same full kernel used for Cov below is what
# makes the resulting covariance diagonal.
K.int <- t(phi) %*% diag(quad$weights) %*% K.bas %*% diag(quad$weights) %*% phi

# Solves the (ordinary) e-val problem in eq (5) of Oya et al., where
# K.int = A and C = identity since the basis functions are orthonormal.
# symmetric = TRUE avoids spurious complex eigenvectors from tiny asymmetries.
E <- eigen(K.int, symmetric = TRUE) # columns are eigenvectors in R^L

# V[j, ] = j-th discretized eigenfunction evaluated at all diagonal nodes
V <- t(E$vectors) %*% t(phi) # L x n.nodes

# Covariance in the eigenbasis. By orthogonality of the eigenfunctions this is
# diagonal to machine precision; enforce it explicitly.
Cov <- V %*% diag(quad$weights) %*% K.bas %*% diag(quad$weights) %*% t(V)
Cov <- diag(diag(Cov))

# Diagonal => no nugget and no eigendecomposition needed to invert.
E.cov <- list(values = diag(Cov), vectors = diag(L))
Cov.inv <- E.cov$vectors %*% diag(1 / abs(E.cov$values)) %*% t(E.cov$vectors)


################################################################################

# Use Cov.inv to define the updated GP prior model, K.post and mean.post.
#
# Because the whole boundary is integrated at once (single set of nodes xg with
# weights quad$weights), there is no need to loop over / block by edge as in
# pkGP-section51.R. Each boundary "feature map" is a single n x L matrix:
#   K2(X, xg) %*% diag(quad$weights) %*% t(V)
# and Cov (hence Cov.inv) is a single L x L (diagonal) matrix.

# Approximate updated (projected) kernel with the "nuggetized" RKHS inner
# product. K.x %*% Cov.inv %*% t(K.y) = <k_x, k_y>_{H(T_0)}, so
# Cov.inv = K_{T_0}^{-1}, where K_{T_0} is the projection of K onto H_0^\perp.
K.post <- function(X, Y) {
  K.x <- K2(X, xg) %*% diag(quad$weights) %*% t(V)
  K.y <- K2(Y, xg) %*% diag(quad$weights) %*% t(V)
  K2(X, Y) - K.x %*% Cov.inv %*% t(K.y)
} # End fn K.post


# Corresponds to eq (13) in the RKHS text, assuming prior mean = 0.
mean.post <- function(X, g) {
  g.x <- matrix(apply(xg, 1, g), nrow = 1) # g evaluated at boundary nodes
  K.e <- K2(X, xg) %*% diag(quad$weights) %*% t(V)
  K.g <- g.x %*% diag(quad$weights) %*% t(V)
  K.e %*% Cov.inv %*% t(K.g)
} # End fn mean.post


N <- 20 # Number of interior observations = 10*d
D <- 2 * randomLHS(n = N, k = 2) - 1 # LHS in 2D, converting [0, 1] -> [-1, 1]


reg.bayes <- seq(-1, 1, length.out = length(basis)) # Sequence of points
reg.vals <- cbind(reg.bayes, reg.bayes) # Convert to 2D boundary points
D.reg <- rbind(D, reg.vals) # Combine boundary points with interior points.
# These are the training points for the usual kriging predictor - analagous
# the 4*L nodes around the boundary that were used to compute the spectral
# projection

t <- seq(-0.9, 0.9, 0.1)
#test.b <- t(sapply(t,l)) # Again, comes back around to repeat [-1, 1]
test.b <- cbind(t, t) # Convert to 2D boundary points
n.test <- dim(test.b)[1]
# Pushing away from the boundary
test <- rbind(cbind(t, t - 0.1), cbind(t, t + 0.1))


y <- apply(D, 1, f) #+ rnorm(N, mean=0, sd=0.05)  # Just the interior points
y.reg <- apply(D.reg, 1, f) #+ rnorm((N + n.reg), mean=0, sd=0.05)  # Interior and
# boundary observations

########## pkGP ##############

mu.post <- mean.post(test, f) # Computing eq. (13) for t = test points. This,
# and eq (15) determine the prior GP: GP(\mu_0, k_0)

# mu.post and K.post determine the (projected) GP "PRIOR" distribution. Which
# is updated to obtain the posterior predictive distribuiton. Below is the
# mean of this distribution == the linear predictor using pGP
# post.mean <- mu.post +
#     K.post(test, D)%*%solve(K.post(D, D) + 0.05*diag(dim(D)[1]))%*%
#     (y-mean.post(D,f))

post.mean <- mu.post +
  K.post(test, D) %*% solve(K.post(D, D)) %*% (y - mean.post(D, f)) # We only need the k_0(.,.) at the interior points
# because the boundary has already been absorbed
# into K.post (the basis functions and the i.p.
# on H(T_0)

########### ordinary kriging ############
# Ordinary kriging fails completely so no need to try it

# The typical predictor by conditioning on all points, including those
# at the boundary
#E.reg=eigen(k(D.reg,D.reg))#+1e-10*diag(N+n.reg));
#E.reg <- eigen(K2(D.reg, D.reg) + 0.05*diag(dim(D.reg)[1]))
# E.reg <- eigen(K2(D.reg, D.reg))
# K.reg.inv <- E.reg$vectors%*%diag(1/abs(E.reg$values))%*%t(E.reg$vectors)
# #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
# post.mean.reg <- K2(test,D.reg)%*%K.reg.inv%*%y.reg  # The usual predictor

# Pseudo-kriging
E.reg.nug <- eigen(K2(D.reg, D.reg) + 1e-6 * diag(dim(D.reg)[1]))
#E.reg <- eigen(K2(D.reg, D.reg))
K.reg.inv.nug <- E.reg.nug$vectors %*%
  diag(1 / abs(E.reg.nug$values)) %*%
  t(E.reg.nug$vectors)
#post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
post.mean.reg.nug <- K2(test, D.reg) %*% K.reg.inv.nug %*% y.reg # The usual pred.


f.test <- apply(test, 1, f) # True values


# Metrics
# cor(f.test, post.mean)
# cor(f.test, post.mean.reg.nug)

# # Relative errors
# max(abs(f.test-post.mean))/max(f.test)
# max(abs(f.test-post.mean.reg.nug))/max(f.test)

sqrt(sum((f.test - post.mean)^2) / dim(test)[1])
sqrt(sum((f.test - post.mean.reg.nug)^2) / dim(test)[1])


# Relative L2 error
# sqrt(sum(((f.test-post.mean)^2)))/sqrt(sum(f.test^2))
# sqrt(sum(((f.test-post.mean.reg.nug)^2)))/sqrt(sum(f.test^2))

###############################################################################
############################# Visuals           ###############################

#### Plot Figure 7
x.vals <- seq(-1, 1, length = 50)
y.vals <- seq(-1, 1, length = 50)
grid <- expand.grid(x.vals, y.vals)

z <- apply(grid, 1, f)

x.mat <- matrix(
  rep(x.vals, times = length(y.vals)),
  ncol = length(x.vals),
  byrow = TRUE
)
y.mat <- matrix(
  rep(y.vals, times = length(x.vals)),
  ncol = length(y.vals),
  byrow = TRUE
)
z.mat <- matrix(z, ncol = length(x.vals), nrow = length(y.vals), byrow = TRUE)

#pdf("diagonalExPerspPlot.pdf", width = 25, height = 25)
x11(width = 25, height = 25)
persp3D(
  x.vals,
  y.vals,
  z.mat,
  colvar = z.mat,
  col = viridis(100),
  alpha = 0.7,
  theta = 40,
  phi = 30,
  expand = 0.7,
  xlab = "x",
  ylab = "y",
  zlab = "f(x, y)",
  main = "",
  contour = TRUE,
  colkey = FALSE,
  lighting = TRUE,
  ticktype = "detailed",
  axes = TRUE,
  cex.lab = 3, # Increases size of x, y, z labels
  cex.axis = 2.6, # Increases size of axis tick labels
  cex.main = 4
)


points3D(
  test[, 2],
  test[, 1],
  post.mean.reg.nug,
  add = TRUE,
  pch = 25,
  cex = 2.0,
  col = "blue"
)


points3D(
  test[, 2],
  test[, 1],
  post.mean,
  add = TRUE,
  pch = 8,
  cex = 2.0,
  col = "red"
)


points3D(
  D.reg[, 2],
  D.reg[, 1],
  y.reg,
  add = TRUE,
  pch = 20,
  cex = 2.0,
  col = "black"
)


legend(
  "topleft",
  legend = c("pkGP", "Pseudo-Kriging"),
  col = c("red", "blue"),
  pch = c(8, 25),
  cex = 4,
  bty = "n"
)
dev.off()


# Figure 8
#pdf("diagonalExPredsVsTrue.pdf", width= 30, height= 15)
x11(width = 30, height = 15)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2) + 0.1)
plot(
  post.mean,
  f.test,
  xlab = "Prediction",
  ylab = "True f",
  pch = 20,
  cex = 3,
  main = "pkGP",
  cex.lab = 3, # Increases size of x, y, z labels
  cex.axis = 2.6, # Increases size of axis tick labels
  cex.main = 4
)
abline(a = 0, b = 1, lwd = 2)
legend(
  "topleft",
  legend = sprintf(
    "RMSE = %.4f",
    sqrt(sum((f.test - post.mean)^2) / dim(test)[1])
  ),
  cex = 4,
  bty = "n"
)

plot(
  post.mean.reg.nug,
  f.test,
  pch = 20,
  xlab = "Prediction",
  cex = 3,
  ylab = "",
  main = "Pseudo-Kriging",
  cex.lab = 3, # Increases size of x, y, z labels
  cex.axis = 2.6, # Increases size of axis tick labels
  cex.main = 4
)
abline(a = 0, b = 1, lwd = 2)
legend(
  "topleft",
  legend = sprintf(
    "RMSE = %.4f",
    sqrt(sum((f.test - post.mean.reg.nug)^2) / dim(test)[1])
  ),
  cex = 4,
  bty = "n"
)
dev.off()


################################################################################
################################################################################

# Do ordinary kriging with an increasing number of pseudo-observations and plot
# results vs. pkGP

# No. of pseudo-observations is in reg.vals

r.errs <- rep(0, 10)
r.errs.full <- rep(0, 10)
n.ps <- seq(6, 55, length = 10) # Number of pseudo-observations


x.vals <- seq(-1, 1, length = 30)
y.vals <- seq(-1, 1, length = 30)
test.grid <- as.matrix(expand.grid(x.vals, y.vals))
f.test.grid <- apply(test.grid, 1, f) # True values

for (i in 1:10) {
  reg.bayes.it <- seq(-1, 1, length.out = n.ps[i]) # Sequence of points
  reg.vals.it <- cbind(reg.bayes.it, reg.bayes.it) # Convert to 2D boundary points
  D.reg.it <- rbind(D, reg.vals.it)

  y.reg.it <- apply(D.reg.it, 1, f)

  # E.reg <- eigen(K2(D.reg, D.reg))
  # K.reg.inv <- E.reg$vectors%*%diag(1/abs(E.reg$values))%*%t(E.reg$vectors)
  # #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
  # post.mean.reg <- K2(test,D.reg)%*%K.reg.inv%*%y.reg  # The usual predictor

  # Try with and without nugget
  E.reg.nug <- eigen(K2(D.reg.it, D.reg.it) + 1e-6 * diag(dim(D.reg.it)[1]))
  #E.reg <- eigen(K2(D.reg, D.reg))
  K.reg.inv.nug <- E.reg.nug$vectors %*%
    diag(1 / abs(E.reg.nug$values)) %*%
    t(E.reg.nug$vectors)
  #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
  post.mean.reg.nug <- K2(test, D.reg.it) %*% K.reg.inv.nug %*% y.reg.it

  post.mean.reg.full <- K2(test.grid, D.reg.it) %*% K.reg.inv.nug %*% y.reg.it

  #r.errs[i] <- max(abs(f.test-post.mean.reg.nug))/max(f.test)
  r.errs[i] <- sqrt(sum(((f.test - post.mean.reg.nug)^2))) / sqrt(sum(f.test^2))
  #r.errs.full[i] <- max(abs(f.test.grid-post.mean.reg.full))/max(f.test.grid)
  r.errs.full[i] <- sqrt(sum(((f.test.grid - post.mean.reg.full)^2))) /
    sqrt(sum(f.test.grid^2))
}


mu.post.grid <- mean.post(test.grid, f) # Computing eq. (13) for t = test points. This,
# and eq (15) determine the prior GP: GP(\mu_0, k_0)

# mu.post and K.post determine the (projected) GP "PRIOR" distribution. Which
# is updated to obtain the posterior predictive distribuiton. Below is the
# mean of this distribution == the linear predictor using pGP
# post.mean <- mu.post +
#     K.post(test, D)%*%solve(K.post(D, D) + 0.05*diag(dim(D)[1]))%*%
#     (y-mean.post(D,f))

post.mean.grid <- mu.post.grid +
  K.post(test.grid, D) %*% solve(K.post(D, D)) %*% (y - mean.post(D, f)) # We only need the k_0(.,.) at the interior points
# because the boundary has already been absorbed
# into K.post (the basis functions and the i.p.
# on H(T_0)

# Figure 9
#pdf("diagonalExPseudos.pdf", width = 25, height = 25)
x11(width = 25, height = 25)
#par(mfrow= c(1, 2), mar = c(5, 5, 4, 2) + 0.1)
par(mar = c(5, 5, 4, 2) + 0.1)


plot(
  n.ps,
  r.errs.full,
  type = "b",
  pch = 20,
  cex = 3,
  ylim = c(0, 0.8),
  lwd = 3,
  xlab = "Number of Pseudo-observations",
  ylab = "Relative L2 Error",
  main = "",
  cex.lab = 3, # Increases size of x, y, z labels
  cex.axis = 2.6, # Increases size of axis tick labels
  cex.main = 4
)
abline(
  h = sqrt(sum(((f.test.grid - post.mean.grid)^2))) / sqrt(sum(f.test.grid^2)),
  lty = 2,
  lwd = 3
)
legend(
  "topright",
  legend = c("pkGP", "Pseudo-Kriging"),
  lty = c(2, 1),
  lwd = 3,
  bty = "n",
  cex = 4
)
dev.off()

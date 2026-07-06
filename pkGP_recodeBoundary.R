# pkGP_recodeBoundary.R
# Andrew Brown
# 18 June 2026

# This file is to replicate the pkGP-section51 boundary example, but with a
# cleaner parameterization of the boundary to simplify inner product computation
#
#
# minor changes to check git tracking.
################################################################################

rm(list = ls())

library(orthopolynom)
library(statmod)
library(lhs)
library(plot3D)
library(viridis)
library(fields)

#sessionInfo()

set.seed(420) # Same as in *-section51.R to verify replication

################################################################################

# Define functions

# Lim et al. nonpolynomial function. (Using code directly from
# https://www.sfu.ca/~ssurjano/Code/limetal02nonr.html)
limetal02non <- function(xx) {
  ##########################################################################
  #
  # LIM ET AL. (2002) FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it
  # and/or modify it under the terms of the GNU General Public License as
  # published by the Free Software Foundation; version 2.0 of the License.
  # Accordingly, this program is distributed in the hope that it will be
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2)
  #
  # The function as defined assumes x \in [0, 1]^2, but the input is
  # x \in [-1, 1]^2, so needs to be rescaled.
  #
  ##########################################################################

  # Map [-1, 1] -> [0, 1]
  xs <- xx
  xs[1] <- (xx[1] + 1) / 2
  xs[2] <- (xx[2] + 1) / 2

  x1 <- xs[1]
  x2 <- xs[2]

  fact1 <- 30 + 5 * x1 * sin(5 * x1)
  fact2 <- 4 + exp(-5 * x2)

  y <- (fact1 * fact2 - 100) / 6
  return(y)
} # End fn limetal02non


# Covariance function
# K2 <- function(X, Y) {
#   # X \in R^2, Y \in R^2

#   # out= params[1]*kern.x(X[, 1], Y[, 1])*kern.y(X[, 2], Y[, 2])

#   # Matern kernel w/ nu = 3/2, correlation length = sqrt(params[2])
#   dist.mat = as.matrix(dist(rbind(X, Y)))
#   pairwise_dist = dist.mat[1:nrow(X), (nrow(X) + 1):(nrow(X) + nrow(Y))]
#   len.par = sqrt(params[2])

#   out = params[1] *
#     (1 + (pairwise_dist * sqrt(3) / len.par)) *
#     exp(-1 * pairwise_dist * sqrt(3) / len.par)

#   return(out)
# } # End fn K2
K2 <- function(X, Y) {
  params[1] * Matern(rdist(X, Y), range = sqrt(params[2]), nu = 3 / 2)
}


# Boundary parameterization: s in [-1, 1] traces the perimeter of [-1, 1]^2
# counter-clockwise, with corners at s = -1, -0.5, 0, 0.5 (back to -1 at s=1).
#
#   s in [-1,  -0.5)  ->  bottom edge: (4s+3,  -1),  x: -1 ->  1,  y = -1
#   s in [-0.5,  0)   ->  right edge:  (1,  4s+1),   x =  1,       y: -1 ->  1
#   s in [0,    0.5)  ->  top edge:    (1-4s,   1),  x:  1 -> -1,  y =  1
#   s in [0.5,   1]   ->  left edge:   (-1, 3-4s),   x = -1,       y:  1 -> -1
l <- function(s) {
  if (s >= -1 & s < -0.5) {
    return(c(4 * s + 3, -1))
  } else if (s >= -0.5 & s < 0) {
    return(c(1, 4 * s + 1))
  } else if (s >= 0 & s < 0.5) {
    return(c(1 - 4 * s, 1))
  } else {
    return(c(-1, 3 - 4 * s))
  }
} # End fn l


# Inverse: map a 2D boundary point back to s in [-1, 1]
l.inv <- function(x) {
  if (x[2] == -1) {
    return((x[1] - 3) / 4) # Bottom edge
  } else if (x[1] == 1) {
    return((x[2] - 1) / 4) # Right edge
  } else if (x[2] == 1) {
    return((1 - x[1]) / 4) # Top edge
  } else {
    return((3 - x[2]) / 4) # Left edge
  }
} # End fn l.inv


################################################################################

J <- 15 # Max order of polynomial approximation for RR method

# In the original file, we used m = 15 quadrature nodes along each of 4 edges.
# This is equivalent to using 4*15 = 60 nodes around the entire square.

params <- c(1, 1, 1) # (scale, length-scale x, length-scale y)

m <- 15 # Polynomial order (basis has m+1 functions)

# Gauss-Legendre nodes and weights on [-1, 1] covering the entire boundary.
# Using 4*m nodes total; increase for finer resolution along the boundary.
quad <- gauss.quad(4 * m, kind = "legendre")
n.nodes <- length(quad$nodes)

# Map quadrature nodes to 2D boundary points: single n.nodes x 2 matrix
xg <- t(sapply(quad$nodes, l))

# Legendre polynomial basis, evaluated at the boundary quadrature nodes
basis <- lapply(legendre.polynomials(m, normalized = TRUE), as.function)
L <- length(basis) # = m + 1

phi <- matrix(0, nrow = n.nodes, ncol = L)
for (i in 1:L) {
  phi[, i] <- sapply(quad$nodes, basis[[i]])
}

# Unified kernel matrix over all boundary quadrature points (n.nodes x n.nodes)
K.bas <- K2(xg, xg)

# Discretized boundary integral operator: phi^T W K.bas W phi  (L x L)
# Eigenvectors of this are proper eigenfunctions of the full boundary kernel
# operator
K.int <- t(phi) %*% diag(quad$weights) %*% K.bas %*% diag(quad$weights) %*% phi # K.int = A matrix in Oya
# = {<K\varphi_i, \varphi_j>}

# Solves the (ordinary) e-val problem in eq (5) of
# Oya et al., where K.int = A and C = identity since basis fns are
# orthonormal. E$vectors contains the e-vectors as COLUMN vectors.
E <- eigen(K.int) # E$vectors: columns are eigenvectors in R^L


# V[j, ] = j-th discretized eigenfunction evaluated at all n.nodes boundary points
V <- t(E$vectors) %*% t(phi) # L x n.nodes

# Boundary covariance matrix: now a single L x L block (vs. 4L x 4L in original)
Cov <- V %*% diag(quad$weights) %*% K.bas %*% diag(quad$weights) %*% t(V)

# By orthogonal eigenfunctions, we know Cov is diagonal. Code it as such.
Cov <- diag(diag(Cov))


# I don't think we need the nugget if the matrix is diagonal.
#E.cov <- eigen(Cov + 1e-6 * diag(L), symmetric = TRUE)
# E.cov <- eigen(Cov, symmetric = TRUE) don't need to do it at all!
E.cov <- list(values = diag(Cov), vectors = diag(L))
E.cov$values <- diag(Cov)
E.cov$vectors <- diag(L)

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


################################################################################
################# Putting it to work ###########################################
################################################################################

N <- 20 # Number of interior observations = 10*d
D <- 2 * randomLHS(n = N, k = 2) - 1 # LHS in 2D, converting [0, 1] -> [-1, 1]

# We do not use [0, 8] anymore because of how the boundary is coded now. Just
# use [-1, 1]
# reg.bayes <- seq(0, 8, length.out = 4 * length(basis)) # Sequence of points
# to map to the edges of [-1, 1]

reg.bayes <- seq(-1, 1, length.out = 4 * length(basis)) # Sequence of points
n.reg <- length(reg.bayes)
reg.vals <- t(sapply(reg.bayes, l)) # Convert to 2D boundary points

D.reg <- rbind(D, reg.vals) # Combine boundary points with interior points.
# These are the training points for the usual kriging predictor - analagous
# the 4*L nodes around the boundary that were used to compute the spectral
# projection

# Here also using sequence along [-1, 1] instead of [0, 8]
#t <- seq(0, 8, 0.1)
t <- seq(-1, 1, 0.1)
test.b <- t(sapply(t, l)) # Again, comes back around to repeat [-1, 1]
n.test <- dim(test.b)[1]
# Pushing away from the boundary
test <- rbind(0.9 * test.b, 0.5 * test.b)


y <- apply(D, 1, limetal02non) #+ rnorm(N, mean=0, sd=0.05)  # Just the interior points
y.reg <- apply(D.reg, 1, limetal02non) #+ rnorm((N + n.reg), mean=0, sd=0.05)  # Interior
# and boundary observations

############################ pkGP ###################################

mu.post <- mean.post(test, limetal02non) # Computing eq. (13) for
# t = test points. This, and eq (15) determine the prior GP: GP(\mu_0, k_0)

# mu.post and K.post determine the (projected) GP "PRIOR" distribution. Which
# is updated to obtain the posterior predictive distribuiton. Below is the
# mean of this distribution == the linear predictor using pkGP
post.mean <- mu.post +
  K.post(test, D) %*% solve(K.post(D, D)) %*% (y - mean.post(D, limetal02non))
# We only need the k_0(.,.) at the interior points because the boundary has
# already been absorbed into K.post (the basis functions and the i.p. on
# H(T_0)

######################### ordinary kriging ############################

# The typical predictor (no pseduo observations)
# E.reg=eigen(k(D.reg,D.reg))#+1e-10*diag(N+n.reg));
E.reg <- eigen(K2(D, D) + 1e-6 * diag(dim(D)[1]))
# E.reg <- eigen(K2(D.reg, D.reg))
K.reg.inv <- E.reg$vectors %*% diag(1 / abs(E.reg$values)) %*% t(E.reg$vectors)
#post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
post.mean.reg <- K2(test, D) %*% K.reg.inv %*% y # The usual predictor


# Pseudo-kriging with additional training data
E.reg.nug <- eigen(K2(D.reg, D.reg) + 1e-6 * diag(dim(D.reg)[1]))
#E.reg <- eigen(K2(D.reg, D.reg))
K.reg.inv.nug <- E.reg.nug$vectors %*%
  diag(1 / abs(E.reg.nug$values)) %*%
  t(E.reg.nug$vectors)
#post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
post.mean.reg.nug <- K2(test, D.reg) %*% K.reg.inv.nug %*% y.reg # The usual pred.


################################################################################

# Comparing performance

f.test <- apply(test, 1, limetal02non) # True values


# # Metrics
# cor(f.test, post.mean)
# #cor(f.test, bdry.mean)
# #cor(f.test, post.mean.reg)
# cor(f.test, post.mean.reg.nug)
# # summary(lm(f.test~post.mean))
# # summary(lm(f.test~bdry.mean))
# # summary(lm(f.test~post.mean.reg))
#
# max(abs(f.test-post.mean))
# #max(abs(f.test-bdry.mean))
# #max(abs(f.test-post.mean.reg))
# max(abs(f.test-post.mean.reg.nug))

sqrt(sum((f.test - post.mean)^2) / dim(test)[1])
sqrt(sum((f.test - post.mean.reg)^2) / dim(test)[1])
sqrt(sum((f.test - post.mean.reg.nug)^2) / dim(test)[1])


# Analagous to Figure 3 in the pkGP paper
x.vals <- seq(-1, 1, length = 50)
y.vals <- seq(-1, 1, length = 50)
grid <- expand.grid(x.vals, y.vals)

z <- apply(grid, 1, limetal02non)

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


#pdf("boundaryExPerspPlot.pdf", width = 25, height = 25)
x11()
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
  post.mean.reg,
  add = TRUE,
  pch = 14,
  cex = 1.5,
  col = "blue"
)


points3D(
  test[, 2],
  test[, 1],
  post.mean.reg.nug,
  add = TRUE,
  pch = 25,
  cex = 1.5,
  col = "blue"
)


points3D(
  test[, 2],
  test[, 1],
  post.mean,
  add = TRUE,
  pch = 8,
  cex = 1.5,
  col = "red"
)


points3D(
  D.reg[, 2],
  D.reg[, 1],
  y.reg,
  add = TRUE,
  pch = 20,
  cex = 1.5,
  col = "black"
)


legend(
  "topleft",
  legend = c("pkGP", "Kriging", "Pseudo-Kriging"),
  col = c(
    "red",
    #"green",
    "blue",
    "blue"
  ),
  pch = c(
    8,
    #13,
    14,
    25
  ),
  cex = 4,
  bty = "n"
)
dev.off()

# pkGP_recodeBoundary.R
# Andrew Brown
# 18 June 2026

# This file is to replicate the pkGP-section51 boundary example, but with a
# cleaner parameterization of the boundary to simplify inner product computation
#
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
  ##########################################################################

  x1 <- xx[1]
  x2 <- xx[2]

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
# Eigenvectors of this are proper eigenfunctions of the full boundary kernel operator
K.int <- t(phi) %*% diag(quad$weights) %*% K.bas %*% diag(quad$weights) %*% phi

E <- eigen(K.int) # E$vectors: columns are eigenvectors in R^L

# V[j, ] = j-th discretized eigenfunction evaluated at all n.nodes boundary points
V <- t(E$vectors) %*% t(phi) # L x n.nodes

# Boundary covariance matrix: now a single L x L block (vs. 4L x 4L in original)
Cov <- V %*% diag(quad$weights) %*% K.bas %*% diag(quad$weights) %*% t(V)

E.cov <- eigen(Cov + 1e-6 * diag(L), symmetric = TRUE)
Cov.inv <- E.cov$vectors %*% diag(1 / abs(E.cov$values)) %*% t(E.cov$vectors)

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
K2 <- function(X, Y) {
  # X \in R^2, Y \in R^2

  # out= params[1]*kern.x(X[, 1], Y[, 1])*kern.y(X[, 2], Y[, 2])

  # Matern kernel w/ nu = 3/2, correlation length = sqrt(params[2])
  dist.mat = as.matrix(dist(rbind(X, Y)))
  pairwise_dist = dist.mat[1:nrow(X), (nrow(X) + 1):(nrow(X) + nrow(Y))]
  len.par = sqrt(params[2])

  out = params[1] *
    (1 + (pairwise_dist * sqrt(3) / len.par)) *
    exp(-1 * pairwise_dist * sqrt(3) / len.par)

  return(out)
} # End fn K2


################################################################################

J <- 15 # Max order of polynomial approximation for RR method

# In the original file, we used m = 15 quadrature nodes along each of 4 edges.
# This is equivalent to using 4*15 = 60 nodes around the entire square.

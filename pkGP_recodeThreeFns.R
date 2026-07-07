# pkGP-section51-ThreeFns.R
# Andrew Brown and John Nicholson
# Sept. 19, 2025.
#
# This file is to replicate pkGP-section51-ThreeFns.R, but using the cleaner
# boundary representation to facilitate a diagonal cov matrix.
################################################################################

rm(list = ls())


library(orthopolynom)
library(statmod)
library(lhs)
library(plot3D)
library(viridis)
library(fields)
#sessionInfo()

set.seed(420)

################################################################################

# Define the functions
f.corn <- function(x) {
  # The 'corner peak' function with d = 2

  xs <- x

  # Map [-1, 1] -> [0, 1]
  xs[1] <- (x[1] + 1) / 2
  xs[2] <- (x[2] + 1) / 2

  result <- (1 + (1 / 2) * (xs[1] + xs[2]))^(-2 - 1)

  return(result)
} # End f.corn


f.prod <- function(x) {
  # The 'product peak' function with d = 2

  xs <- x

  # Map [-1, 1] -> [0, 1]
  xs[1] <- (x[1] + 1) / 2
  xs[2] <- (x[2] + 1) / 2

  result <- ((1 + 10 * (xs[1] - 0.25)^2)^(-1)) *
    ((1 + 10 * (xs[2] - 0.25)^2)^(-1))

  return(result)
} # End f.prod


f.rosen <- function(x) {
  # The Rosenbrock function with d = 2 dimensions

  xs <- x

  # Map [-1, 1] -> [0, 1]
  xs[1] <- (x[1] + 1) / 2
  xs[2] <- (x[2] + 1) / 2

  result <- 100 * (xs[2] - xs[1]^2)^2 + (1 - xs[1])^2

  return(result)
} # End f.rosen


# Covariance function (Note that params are defined below.)
K2 <- function(X, Y) {
  params[1] * Matern(rdist(X, Y), range = sqrt(params[2]) / sqrt(3), nu = 3 / 2)
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


################################################################################

################################################################################

params <- c(1, 1, 1) # (scale, length-scale x, length-scale y)

m <- 15 # Polynomial order PER EDGE (basis has m+1 functions per edge,
# L = 4*(m+1) total)

n.per <- 60 # Quadrature nodes per edge (matches pkGP-section51.R, which used
# a 60-node rule on each edge)

# --- Per-edge composite Gauss-Legendre quadrature -----------------------------
#
# The boundary functions s |-> K(x, l(s)) are smooth along each edge but have
# kinks at the corners (s = -0.5, 0, 0.5). A single Gauss rule spanning the
# corners converges poorly, so instead we use four Gauss panels, one per edge,
# with panel boundaries at the corners. Each edge occupies an s-subinterval of
# half-width 0.25 centered at:
centers <- c(-0.75, -0.25, 0.25, 0.75)

gq <- gauss.quad(n.per, kind = "legendre") # Reference rule on [-1, 1]

# Composite rule over s in [-1, 1]: affinely map the reference rule into each
# edge's subinterval (Jacobian = 0.25)
quad <- list(
  nodes = as.vector(outer(0.25 * gq$nodes, centers, "+")),
  weights = rep(0.25 * gq$weights, 4)
)
n.nodes <- length(quad$nodes) # = 4 * n.per

# Map quadrature nodes to 2D boundary points: single n.nodes x 2 matrix
xg <- t(sapply(quad$nodes, l))

# --- Piecewise (per-edge) Legendre basis --------------------------------------
#
# 16 normalized Legendre polynomials on each edge, zero off the edge, giving
# L = 64 basis functions total. Piecewise polynomials capture the corner kinks
# that global polynomials on [-1, 1] cannot, and keep the polynomial degree low
# (high-degree Legendre evaluation via orthopolynom is numerically unstable).
basis <- lapply(legendre.polynomials(m, normalized = TRUE), as.function)
L <- 4 * length(basis) # = 4 * (m + 1)

# On an edge subinterval of half-width 0.25, the normalized basis function is
# 2 * P_k(u), where u = 4*(s - center) maps the subinterval to [-1, 1] and the
# factor 2 = 1/sqrt(0.25) restores unit L2(ds) norm. Since the quadrature
# nodes are the mapped reference nodes, u = gq$nodes exactly, so each edge's
# block of phi is the same (n.per x (m+1)) matrix.
phi.edge <- 2 * sapply(basis, function(b) sapply(gq$nodes, b))
phi <- matrix(0, nrow = n.nodes, ncol = L)
for (j in 1:4) {
  phi[
    ((j - 1) * n.per + 1):(j * n.per),
    ((j - 1) * (m + 1) + 1):(j * (m + 1))
  ] <- phi.edge
}

# Unified kernel matrix over all boundary quadrature points (n.nodes x n.nodes)
K.bas <- K2(xg, xg)

# Discretized boundary integral operator: phi^T W K.bas W phi  (L x L)
# Eigenvectors of this are proper eigenfunctions of the full boundary kernel
# operator
#
# K.int = A matrix in Oya = {<K\varphi_i, \varphi_j>}
K.int <- t(phi) %*% diag(quad$weights) %*% K.bas %*% diag(quad$weights) %*% phi


# Solves the (ordinary) e-val problem in eq (5) of
# Oya et al., where K.int = A and C = identity since basis fns are
# orthonormal. E$vectors contains the e-vectors as COLUMN vectors.
# symmetric = TRUE is important: the default general solver can return
# spurious complex eigenvectors from tiny numerical asymmetries.
E <- eigen(K.int, symmetric = TRUE) # E$vectors: columns are eigenvectors in R^L


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
#E.cov$values <- diag(Cov)
#E.cov$vectors <- diag(L)

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

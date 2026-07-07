# pkGP-section51-ThreeFns.R
# Andrew Brown and John Nicholson
# July 7, 2026
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


################# Test with training data ######################################
################################################################################
################################################################################

# We are going to loop over training data size, but keeping the number of
# psuedo-observations the same for ordinary kriging

# corner peak fn

#Ns <- c(20, 40, 60, 80, 100)  # Number of interior observations
Ns <- seq(10, 200, by = 10) # Number of interior observations


pk.err <- rep(0, length(Ns))
ps.krig <- rep(0, length(Ns))
ps.k.reg <- rep(0, length(Ns))

for (sz.ind in 1:length(Ns)) {
  N <- Ns[sz.ind] # Number of interior observations = 10*d
  D <- 2 * randomLHS(n = N, k = 2) - 1 # LHS in 2D, converting [0, 1] -> [-1, 1]

  reg.bayes <- seq(-1, 1, length.out = 4 * length(basis)) # Sequence of points
  # to map to the edges of [-1, 1]
  n.reg <- length(reg.bayes)
  reg.vals <- t(sapply(reg.bayes, l)) # Convert to 2D boundary points
  D.reg <- rbind(D, reg.vals) # Combine boundary points with interior points.
  # These are the training points for the usual kriging predictor - analagous
  # the 4*L nodes around the boundary that were used to compute the spectral
  # projection

  # t <- seq(0, 8, 0.1)
  # test.b <- t(sapply(t,l)) # Again, comes back around to repeat [-1, 1]
  # n.test <- dim(test.b)[1]
  # # Pushing away from the boundary
  # test <- rbind(0.9*test.b, 0.5*test.b)

  # For this example, we want to approximate the function over the entire grid
  x.vals <- seq(-1, 1, length = 30)
  y.vals <- seq(-1, 1, length = 30)
  test.grid <- as.matrix(expand.grid(x.vals, y.vals))

  z <- apply(test.grid, 1, f.corn) # True values

  y <- apply(D, 1, f.corn) #+ rnorm(N, mean=0, sd=0.05)  # Just the interior points
  y.reg <- apply(D.reg, 1, f.corn) #+ rnorm((N + n.reg), mean=0, sd=0.05)  # Interior and
  # boundary observations

  ########## pkGP ##############

  mu.post <- mean.post(test.grid, f.corn) # Computing eq. (13) for t = test points.
  # This, and eq (15) determine the prior GP: GP(\mu_0, k_0)

  # mu.post and K.post determine the (projected) GP "PRIOR" distribution. Which
  # is updated to obtain the posterior predictive distribuiton. Below is the
  # mean of this distribution == the linear predictor using pGP
  # post.mean <- mu.post +
  #     K.post(test, D)%*%solve(K.post(D, D) + 0.05*diag(dim(D)[1]))%*%
  #     (y-mean.post(D,f))
  post.mean <- mu.post +
    K.post(test.grid, D) %*% solve(K.post(D, D)) %*% (y - mean.post(D, f.corn)) # We only need the k_0(.,.) at the interior points
  # because the boundary has already been absorbed
  # into K.post (the basis functions and the i.p.
  # on H(T_0)

  ########### ordinary kriging ############

  # The typical predictor (no pseduo observations)
  # E.reg=eigen(k(D.reg,D.reg))#+1e-10*diag(N+n.reg));
  E.reg <- eigen(K2(D, D) + 1e-6 * diag(dim(D)[1]))
  # E.reg <- eigen(K2(D.reg, D.reg))
  K.reg.inv <- E.reg$vectors %*%
    diag(1 / abs(E.reg$values)) %*%
    t(E.reg$vectors)
  #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
  post.mean.reg <- K2(test.grid, D) %*% K.reg.inv %*% y # The usual predictor

  # Pseudo with additional training data
  E.reg.nug <- eigen(K2(D.reg, D.reg) + 1e-6 * diag(dim(D.reg)[1]))
  #E.reg <- eigen(K2(D.reg, D.reg))
  K.reg.inv.nug <- E.reg.nug$vectors %*%
    diag(1 / abs(E.reg.nug$values)) %*%
    t(E.reg.nug$vectors)
  #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
  post.mean.reg.nug <- K2(test.grid, D.reg) %*% K.reg.inv.nug %*% y.reg # The usual pred.

  f.test <- apply(test.grid, 1, f.corn) # True values

  # Metrics
  # cor(f.test, post.mean)
  # #cor(f.test, post.mean.reg)
  # cor(f.test, post.mean.reg.nug)
  # # summary(lm(f.test~post.mean))
  # # summary(lm(f.test~bdry.mean))
  # # summary(lm(f.test~post.mean.reg))
  #
  # max(abs(f.test-post.mean))
  # #max(abs(f.test-post.mean.reg))
  # max(abs(f.test-post.mean.reg.nug))
  #
  #
  # sqrt(sum((f.test-post.mean)^2)/dim(test)[1])
  # #sqrt(sum((f.test-post.mean.reg)^2)/dim(test)[1])
  # sqrt(sum((f.test-post.mean.reg.nug)^2)/dim(test)[1])

  # Compute the relative sup errors
  pk.err[sz.ind] <- max(abs(f.test - post.mean)) / max(abs(f.test))
  ps.krig[sz.ind] <- max(abs(f.test - post.mean.reg.nug)) / max(abs(f.test))
  ps.k.reg[sz.ind] <- max(abs(f.test - post.mean.reg)) / max(abs(f.test))

  # Compute relative L1 errors
  # pk.err[sz.ind] <- norm(f.test - post.mean, type= "o")/
  #   norm(as.matrix(f.test), type= "o")
  # ps.krig[sz.ind] <- norm(f.test - post.mean.reg.nug, type= "o")/
  #   norm(as.matrix(f.test), type= "o")
} # End loop over training data size


################################################################################
################################################################################
################################################################################

# Repeat the experiment above with the product peak function

# product peak fn

#Ns <- c(20, 40, 60, 80, 100)  # Number of interior observations
Ns <- seq(10, 200, by = 10) # Number of interior observations


pk.err.prod <- rep(0, length(Ns))
ps.krig.prod <- rep(0, length(Ns))
ps.k.reg.prod <- rep(0, length(Ns))

for (sz.ind in 1:length(Ns)) {
  N <- Ns[sz.ind] # Number of interior observations = 10*d
  D <- 2 * randomLHS(n = N, k = 2) - 1 # LHS in 2D, converting [0, 1] -> [-1, 1]

  reg.bayes <- seq(-1, 1, length.out = 4 * length(basis)) # Sequence of points
  # to map to the edges of [-1, 1]
  n.reg <- length(reg.bayes)
  reg.vals <- t(sapply(reg.bayes, l)) # Convert to 2D boundary points
  D.reg <- rbind(D, reg.vals) # Combine boundary points with interior points.
  # These are the training points for the usual kriging predictor - analagous
  # the 4*L nodes around the boundary that were used to compute the spectral
  # projection

  # t <- seq(0, 8, 0.1)
  # test.b <- t(sapply(t,l)) # Again, comes back around to repeat [-1, 1]
  # n.test <- dim(test.b)[1]
  # # Pushing away from the boundary
  # test <- rbind(0.9*test.b, 0.5*test.b)

  # For this example, we want to approximate the function over the entire grid
  x.vals <- seq(-1, 1, length = 30)
  y.vals <- seq(-1, 1, length = 30)
  test.grid <- as.matrix(expand.grid(x.vals, y.vals))

  #z <- apply(test.grid, 1, f.corn)  # True values

  y <- apply(D, 1, f.prod) #+ rnorm(N, mean=0, sd=0.05)  # Just the interior points
  y.reg <- apply(D.reg, 1, f.prod) #+ rnorm((N + n.reg), mean=0, sd=0.05)  # Interior and
  # boundary observations

  ########## pkGP ##############

  mu.post <- mean.post(test.grid, f.prod) # Computing eq. (13) for t = test points.
  # This, and eq (15) determine the prior GP: GP(\mu_0, k_0)

  # mu.post and K.post determine the (projected) GP "PRIOR" distribution. Which
  # is updated to obtain the posterior predictive distribuiton. Below is the
  # mean of this distribution == the linear predictor using pGP
  # post.mean <- mu.post +
  #     K.post(test, D)%*%solve(K.post(D, D) + 0.05*diag(dim(D)[1]))%*%
  #     (y-mean.post(D,f))
  post.mean <- mu.post +
    K.post(test.grid, D) %*% solve(K.post(D, D)) %*% (y - mean.post(D, f.prod)) # We only need the k_0(.,.) at the interior points
  # because the boundary has already been absorbed
  # into K.post (the basis functions and the i.p.
  # on H(T_0)

  ########### ordinary kriging ############

  # The typical predictor (no pseduo observations)
  # E.reg=eigen(k(D.reg,D.reg))#+1e-10*diag(N+n.reg));
  E.reg <- eigen(K2(D, D) + 1e-6 * diag(dim(D)[1]))
  # E.reg <- eigen(K2(D.reg, D.reg))
  K.reg.inv <- E.reg$vectors %*%
    diag(1 / abs(E.reg$values)) %*%
    t(E.reg$vectors)
  #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
  post.mean.reg <- K2(test.grid, D) %*% K.reg.inv %*% y # The usual predictor

  # Pseudo with additional training data
  E.reg.nug <- eigen(K2(D.reg, D.reg) + 1e-6 * diag(dim(D.reg)[1]))
  #E.reg <- eigen(K2(D.reg, D.reg))
  K.reg.inv.nug <- E.reg.nug$vectors %*%
    diag(1 / abs(E.reg.nug$values)) %*%
    t(E.reg.nug$vectors)
  #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
  post.mean.reg.nug <- K2(test.grid, D.reg) %*% K.reg.inv.nug %*% y.reg # The usual pred.

  f.test <- apply(test.grid, 1, f.prod) # True values

  # Metrics
  # cor(f.test, post.mean)
  # #cor(f.test, post.mean.reg)
  # cor(f.test, post.mean.reg.nug)
  # # summary(lm(f.test~post.mean))
  # # summary(lm(f.test~post.mean.reg))
  #
  # max(abs(f.test-post.mean))
  # #max(abs(f.test-post.mean.reg))
  # max(abs(f.test-post.mean.reg.nug))
  #
  #
  # sqrt(sum((f.test-post.mean)^2)/dim(test)[1])
  # #sqrt(sum((f.test-post.mean.reg)^2)/dim(test)[1])
  # sqrt(sum((f.test-post.mean.reg.nug)^2)/dim(test)[1])

  # Compute the relative sup errors
  pk.err.prod[sz.ind] <- max(abs(f.test - post.mean)) / max(abs(f.test))
  ps.krig.prod[sz.ind] <- max(abs(f.test - post.mean.reg.nug)) /
    max(abs(f.test))
  ps.k.reg.prod[sz.ind] <- max(abs(f.test - post.mean.reg)) / max(abs(f.test))
} # End loop over training data size


################################################################################
################################################################################
################################################################################

# Repeat the experiment above with the Rosenbrock function

#Ns <- c(20, 40, 60, 80, 100)  # Number of interior observations
Ns <- seq(10, 200, by = 10) # Number of interior observations


pk.err.rosen <- rep(0, length(Ns))
ps.krig.rosen <- rep(0, length(Ns))
ps.k.reg.rosen <- rep(0, length(Ns))

for (sz.ind in 1:length(Ns)) {
  N <- Ns[sz.ind] # Number of interior observations = 10*d
  D <- 2 * randomLHS(n = N, k = 2) - 1 # LHS in 2D, converting [0, 1] -> [-1, 1]

  reg.bayes <- seq(-1, 1, length.out = 4 * length(basis)) # Sequence of points
  # to map to the edges of [-1, 1]
  n.reg <- length(reg.bayes)
  reg.vals <- t(sapply(reg.bayes, l)) # Convert to 2D boundary points
  D.reg <- rbind(D, reg.vals) # Combine boundary points with interior points.
  # These are the training points for the usual kriging predictor - analagous
  # the 4*L nodes around the boundary that were used to compute the spectral
  # projection

  # t <- seq(0, 8, 0.1)
  # test.b <- t(sapply(t,l)) # Again, comes back around to repeat [-1, 1]
  # n.test <- dim(test.b)[1]
  # # Pushing away from the boundary
  # test <- rbind(0.9*test.b, 0.5*test.b)

  # For this example, we want to approximate the function over the entire grid
  x.vals <- seq(-1, 1, length = 30)
  y.vals <- seq(-1, 1, length = 30)
  test.grid <- as.matrix(expand.grid(x.vals, y.vals))

  #z <- apply(test.grid, 1, f.corn)  # True values

  y <- apply(D, 1, f.rosen) #+ rnorm(N, mean=0, sd=0.05)  # Just the interior points
  y.reg <- apply(D.reg, 1, f.rosen) #+ rnorm((N + n.reg), mean=0, sd=0.05)  # Interior and
  # boundary observations

  ########## pkGP ##############

  mu.post <- mean.post(test.grid, f.rosen) # Computing eq. (13) for t = test points.
  # This, and eq (15) determine the prior GP: GP(\mu_0, k_0)

  # mu.post and K.post determine the (projected) GP "PRIOR" distribution. Which
  # is updated to obtain the posterior predictive distribuiton. Below is the
  # mean of this distribution == the linear predictor using pGP
  # post.mean <- mu.post +
  #     K.post(test, D)%*%solve(K.post(D, D) + 0.05*diag(dim(D)[1]))%*%
  #     (y-mean.post(D,f))
  post.mean <- mu.post +
    K.post(test.grid, D) %*% solve(K.post(D, D)) %*% (y - mean.post(D, f.rosen)) # We only need the k_0(.,.) at the interior points
  # because the boundary has already been absorbed
  # into K.post (the basis functions and the i.p.
  # on H(T_0)

  ########### ordinary kriging ############

  # The typical predictor (no pseduo observations)
  # E.reg=eigen(k(D.reg,D.reg))#+1e-10*diag(N+n.reg));
  E.reg <- eigen(K2(D, D) + 1e-6 * diag(dim(D)[1]))
  # E.reg <- eigen(K2(D.reg, D.reg))
  K.reg.inv <- E.reg$vectors %*%
    diag(1 / abs(E.reg$values)) %*%
    t(E.reg$vectors)
  #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
  post.mean.reg <- K2(test.grid, D) %*% K.reg.inv %*% y # The usual predictor

  # Pseudo with additional training data
  E.reg.nug <- eigen(K2(D.reg, D.reg) + 1e-6 * diag(dim(D.reg)[1]))
  #E.reg <- eigen(K2(D.reg, D.reg))
  K.reg.inv.nug <- E.reg.nug$vectors %*%
    diag(1 / abs(E.reg.nug$values)) %*%
    t(E.reg.nug$vectors)
  #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
  post.mean.reg.nug <- K2(test.grid, D.reg) %*% K.reg.inv.nug %*% y.reg # The usual pred.

  f.test <- apply(test.grid, 1, f.rosen) # True values

  # Metrics
  # cor(f.test, post.mean)
  # #cor(f.test, post.mean.reg)
  # cor(f.test, post.mean.reg.nug)
  # # summary(lm(f.test~post.mean))
  # # summary(lm(f.test~bdry.mean))
  # # summary(lm(f.test~post.mean.reg))
  #
  # max(abs(f.test-post.mean))
  # #max(abs(f.test-post.mean.reg))
  # max(abs(f.test-post.mean.reg.nug))
  #
  #
  # sqrt(sum((f.test-post.mean)^2)/dim(test)[1])
  # #sqrt(sum((f.test-post.mean.reg)^2)/dim(test)[1])
  # sqrt(sum((f.test-post.mean.reg.nug)^2)/dim(test)[1])

  # Compute the relative sup errors
  pk.err.rosen[sz.ind] <- max(abs(f.test - post.mean)) / max(abs(f.test))
  ps.krig.rosen[sz.ind] <- max(abs(f.test - post.mean.reg.nug)) /
    max(abs(f.test))
  ps.k.reg.rosen[sz.ind] <- max(abs(f.test - post.mean.reg)) / max(abs(f.test))

  # Compute relative L1 errors
  # pk.err.rosen[sz.ind] <- norm(f.test - post.mean, type= "o")/
  #   norm(as.matrix(f.test), type= "o")
  # ps.krig.rosen[sz.ind] <- norm(f.test - post.mean.reg.nug, type= "o")/
  #   norm(as.matrix(f.test), type= "o")
} # End loop over training data size


################################################################################
################################################################################
################################################################################

#### Plot the three functions, then plot relative errors vs. training data size
x.vals <- seq(-1, 1, length = 50)
y.vals <- seq(-1, 1, length = 50)
grid <- expand.grid(x.vals, y.vals)

z <- apply(grid, 1, f.corn) # True values

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


# Figure 5
#pdf("threeFunctionsPlot.pdf", width= 60, height= 15)
x11(width = 60, height = 15)
par(mfrow = c(1, 3))
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
  xlab = "",
  ylab = "",
  zlab = "",
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


# points3D(test[, 2], test[, 1], bdry.mean, add= TRUE, pch= 13, cex= 1.5,
#          col= "green")
#
#
# # points3D(test[, 2], test[, 1], post.mean.reg, add= TRUE, pch= 14, cex= 1.5,
# #          col= "blue")
#
#
# points3D(test[, 2], test[, 1], post.mean.reg.nug, add= TRUE, pch= 25, cex= 1.5,
#          col= "blue")
#
#
# points3D(test[, 2], test[, 1], post.mean, add= TRUE, pch= 8, cex= 1.5,
#          col= "red")
#
#
# points3D(D.reg[, 2], D.reg[, 1], y.reg, add= TRUE, pch= 20, cex= 1.5,
#          col= "black")
#
#
# legend("topleft", legend= c("pkGP", "bdryGP", "Kriging"),
#        col= c("red", "green", "blue"),
#        pch= c(8, 13, 25), cex= 4, bty= "n")

z <- apply(grid, 1, f.prod) # True values

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


#x11(width= 45, height= 45)
#pdf("boundaryExPerspPlot.pdf", width= 25, height= 25)
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
  xlab = "",
  ylab = "",
  zlab = "",
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


z <- apply(grid, 1, f.rosen) # True values

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


#x11(width= 45, height= 45)
#pdf("boundaryExPerspPlot.pdf", width= 25, height= 25)
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
  xlab = "",
  ylab = "",
  zlab = "",
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
dev.off()


#################################### Plot the relative errors ##################

# Figure 6
#pdf("threeSupErrsPlot.pdf", width= 60, height= 15)
x11(width = 60, height = 15)
# corner peak
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2) + 0.1)
plot(
  Ns,
  log(ps.krig),
  type = "b",
  pch = 8,
  cex = 2,
  col = "cyan",
  xlab = "Training data size",
  ylab = "Relative sup error (log scale)",
  main = "Corner Peak",
  cex.main = 4,
  ylim = c(-10, 1.5),
  cex.lab = 2.5,
  cex.axis = 2.5,
  lwd = 3
)


lines(Ns, log(ps.k.reg), type = "b", pch = 14, cex = 2, col = "blue", lwd = 3)

lines(Ns, log(pk.err), type = "b", pch = 25, cex = 2, col = "red", lwd = 3)


legend(
  "topright",
  legend = c(
    "pkGP",
    #"bdryGP",
    "Kriging",
    "Pseudo-Kriging"
  ),
  col = c(
    "red",
    #"green",
    "blue",
    "cyan"
  ),
  pch = c(
    8,
    #13,
    14,
    25
  ),
  cex = 2.5,
  bty = "n",
  lwd = 3
)


# product peak
plot(
  Ns,
  log(ps.krig.prod),
  type = "b",
  pch = 8,
  cex = 2,
  col = "cyan",
  xlab = "Training data size",
  ylab = "",
  main = "Product Peak",
  cex.main = 4,
  ylim = c(-10, 1.5),
  cex.lab = 2.5,
  cex.axis = 2.5,
  lwd = 3
)


lines(
  Ns,
  log(ps.k.reg.prod),
  type = "b",
  pch = 14,
  cex = 2,
  col = "blue",
  lwd = 3
)

lines(Ns, log(pk.err.prod), type = "b", pch = 25, cex = 2, col = "red", lwd = 3)

legend(
  "topright",
  legend = c(
    "pkGP",
    #"bdryGP",
    "Kriging",
    "Pseudo-Kriging"
  ),
  col = c(
    "red",
    #"green",
    "blue",
    "cyan"
  ),
  pch = c(
    8,
    #13,
    14,
    25
  ),
  cex = 2.5,
  bty = "n",
  lwd = 3
)


# Rosenbrock
plot(
  Ns,
  log(pk.err.rosen),
  type = "b",
  pch = 8,
  cex = 2,
  col = "red",
  xlab = "Training data size",
  ylab = "",
  main = "Rosenbrock",
  cex.main = 4,
  ylim = c(-10, 1.5),
  cex.lab = 2.5,
  cex.axis = 2.5,
  lwd = 3
)


lines(
  Ns,
  log(ps.k.reg.rosen),
  type = "b",
  pch = 14,
  cex = 2,
  col = "blue",
  lwd = 3
)

lines(
  Ns,
  log(ps.krig.rosen),
  type = "b",
  pch = 25,
  cex = 2,
  col = "cyan",
  lwd = 3
)

legend(
  "topright",
  legend = c(
    "pkGP",
    #"bdryGP",
    "Kriging",
    "Pseudo-Kriging"
  ),
  col = c(
    "red",
    #"green",
    "blue",
    "cyan"
  ),
  pch = c(
    8,
    #13,
    14,
    25
  ),
  cex = 2.5,
  bty = "n",
  lwd = 3
)

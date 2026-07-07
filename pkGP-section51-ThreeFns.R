# pkGP-section51-ThreeFns.R
# Andrew Brown and John Nicholson
# Sept. 19, 2025.
#
# This reproduces the second example in Subsection 5.1 in which we interpolate
# the corner peak, product peak, and Rosenbrock functions.


rm(list= ls())


library(orthopolynom)
library(statmod)
library(lhs)
library(plot3D)
library(viridis)
sessionInfo()

set.seed(420)



################################################################################


# Define the functions
f.corn <- function(x){
  # The 'corner peak' function with d = 2
  
  xs= x
  
  # Map [-1, 1] -> [0, 1]
  xs[1]= (x[1] + 1)/2
  xs[2]= (x[2] + 1)/2
  
  result= (1 + (1/2)*(xs[1] + xs[2]))^(-2-1)
  
  return(result)
} # End f.corn



f.prod <- function(x){
  # The 'product peak' function with d = 2
  
  xs= x
  
  # Map [-1, 1] -> [0, 1]
  xs[1]= (x[1] + 1)/2
  xs[2]= (x[2] + 1)/2
  
  result= ((1 + 10*(xs[1] - 0.25)^2)^(-1))*
    ((1 + 10*(xs[2] - 0.25)^2)^(-1))
  
  return(result)
  
}  # End f.prod



f.rosen <- function(x){
  # The Rosenbrock function with d = 2 dimensions
  
  xs= x
  
  # Map [-1, 1] -> [0, 1]
  xs[1]= (x[1] + 1)/2
  xs[2]= (x[2] + 1)/2
  
  result= 100*(xs[2] - xs[1]^2)^2 + (1 - xs[1])^2
  
  return(result)
  
}  # End f.rosen



# Boundary parameterization
l <- function(t) {
  
  if (t <2 & t>=0) {
    
    return(c(t-1, -1))
    
  } else if(t <4 & t>=2) {
    
    return(c(1, -3+t))
    
  } else if (t<6 & t>=4) {
    
    return(c(5-t, 1))
    
  } else {
    
    return(c(-1, 7-t))
    
  }
  
}  # End l




################################################################################


J <- 15  # Max order of polynomial approximation for RR method

# This will hold the sup norm relative errors and RMSE's associated with 
# Legendre polyn. approximations of varying orders, given by J
leg.max <- rep(0, J-4)  # Legendre?
leg.rmse <- leg.max
#leg.rmse=leg.max=rep(0,J-4);


params <- c(1, 1, 1)  # (scale param., length-scale in x, length-scale in y)


# Will bury this in a loop later, looping over polyn. order m
m <- 15



K2 <- function(X, Y){
  # X \in R^2, Y \in R^2
  
  # out= params[1]*kern.x(X[, 1], Y[, 1])*kern.y(X[, 2], Y[, 2])
  
  # Matern kernel w/ nu = 3/2, correlation length = sqrt(params[2])
  dist.mat= as.matrix(dist(rbind(X, Y)))
  pairwise_dist= dist.mat[1:nrow(X), (nrow(X)+1):(nrow(X) + nrow(Y))]
  len.par= sqrt(params[2])
  
  out= params[1]*(1 + (pairwise_dist*sqrt(3)/len.par))*
    exp(-1*pairwise_dist*sqrt(3)/len.par)
  
  return(out)
  
}  # End fn K2


# John was using 3*m, for one edge? It's not clear why he was doing that. Use
# these nodes for one edge, but repeated around each edge ([-1, 1])
quad <- gauss.quad(4*m, kind= "legendre")
# quad <- gauss.quad(4*m, kind= "chebyshev1")  
n.nodes <- length(quad[[1]])


#list of legendre polynomials
leg <- lapply(legendre.polynomials(m, normalized=TRUE), as.function)
basis <-leg
# cheb.base=lapply(chebyshev.c.polynomials(m,normalized=T),as.function);
# var.change=function(f,scale,f.arg) {
#   function(x) scale*f(f.arg(x));
# }
# 
# cheb <- lapply(cheb.base,function(x) var.change(x,2^(3/2),function(y) 2*y))
# basis <- cheb


# # xg = 'x grid'
xg <- list()
xg[["1"]] <- cbind(quad[[1]], -1)
xg[["2"]] <- cbind(1, quad[[1]]) #goes left to right and down to up
xg[["3"]] <- cbind(quad[[1]], 1)
xg[["4"]] <- cbind(-1, quad[[1]])


# I think this loop is finding the pairwise cov. matrices for the different
# combinations of boundary points (side 1, side 1), (side 1, side 2), ...
# and naming them accordingly.
K.bas=list() 
for (i in 1:4) {
  
  for (j in 1:4) {
    
    name <- paste(i,j,sep="")  # Naming each combination: "11", "12", "13",...
    
    #K.bas[[name]] <- k(xg[[paste(i)]], xg[[paste(j)]])
    K.bas[[name]] <- K2(xg[[paste(i)]], xg[[paste(j)]])
    
    # Add a nugget if pairing the same side since they are so close together
    # => near rank deficiency
    if (i==j) {
      
      K.bas[[name]] <- K.bas[[name]] #+ 1e-6*diag(dim(K.bas[[name]])[1])
      
    }
  }
}




# Now evaluate the basis functions (l.c.'s of which are used to approximate
# the e-vects of K.bas) at the nodes. Only want the nodes along one edge.
L <- length(basis)
phi <- matrix(0,ncol=L,nrow=length(quad[[1]]))
for (i in 1:L) {
  
  phi[,i] <- sapply(quad[[1]],basis[[i]])  # Evaluate the basis functions at
  # the nodes
  
}

K.phi <- K.bas[['11']]%*%diag(quad[[2]])%*%phi
K.int <- t(phi)%*%diag(quad[[2]])%*%K.phi  # \approx A matrix of Oya et al.,
# but only along one edge of the boundary



E=eigen(K.int)  # Solves the (ordinary) e-val problem in eq (5) of 
# Oya et al., where K.int = A and C = identity since basis fns are
# orthonormal. E$vectors contains the e-vectors as COLUMN vectors.

###contains eigenvectors of integral operator evaluated at quad[[1]]
V=t(E$vectors)%*%t(phi)  # jth row of V = discretized appx e-fn j of the
# integral operator, evaluated at x_1, ..., x_n. These guys are orthogonal,
# but not normalized.


Cov <- matrix(0, nrow=4*length(basis), ncol=4*length(basis));
sparse <- 0
for (i in 1:4) {
  
  for (j in 1:4) {
    
    name=paste(i,j,sep="")
    # print(name)
    
    if( abs(i-j)==2 & sparse==1) {
      # "sparse" option zeros out A_{ij} for i,j more than 1 apart. But the
      # labeling is arbitrary, so makes no sense. (Bottom-right pairs [12] 
      # incl., but bottom-left [14] pairs not)
      
      print('independence')
      
    } else {
      # Cov[((i-1)*length(basis)+1):(i*length(basis)),
      #     ((j-1)*length(basis)+1):(j*length(basis))]=K.1=V%*%diag(quad[[2]])%*%K.bas[[name]]%*%diag(quad[[2]])%*%t(V);
      
      
      
      Cov[((i-1)*length(basis)+1):(i*length(basis)),
          ((j-1)*length(basis)+1):(j*length(basis))] <- V%*%diag(quad[[2]])%*%
        K.bas[[name]]%*%diag(quad[[2]])%*%t(V)
      # The e-vects (E) are computed from K.bas[["11"]] above
      
      
      
    }  # End ifelse
  }  # End j loop
}  # End i loop




E.cov <- eigen(Cov + 1e-6*diag(4*length(basis)), symmetric=TRUE)  # Then
# using SVD to invert

Cov.inv <- E.cov$vectors%*%diag(1/abs(E.cov$values))%*%t(E.cov$vectors)


# Use cov.inv to define the updated GP prior model, k.post and mean.post

# This function is the (rough) approximation of, and corresponds to, eq. (15)
# in the RKHS text - the updated kernel
K.post <- function(X,Y) {
  
  n.x= dim(X)[1]
  n.y= dim(Y)[1]
  K.x= matrix(0, nrow=n.x, ncol=4*length(basis))
  K.y= matrix(0, nrow=n.y, ncol=4*length(basis))
  
  for (i in 1:4) {
    
    #k.x= k(X,xg[[paste(i)]])
    k.x= K2(X, xg[[paste(i)]])
    #k.y= k(Y,xg[[paste(i)]])
    k.y= K2(Y, xg[[paste(i)]])
    K.x[,((i-1)*length(basis)+1):(i*length(basis))]= 
      k.x%*%diag(quad[[2]])%*%t(V)
    K.y[,((i-1)*length(basis)+1):(i*length(basis))]= 
      k.y%*%diag(quad[[2]])%*%t(V)
    
  }
  
  #return(k(X,Y) - K.x%*%Cov.inv%*%t(K.y))
  return(K2(X,Y) - K.x%*%Cov.inv%*%t(K.y))  # K.x%*%Cov.inv%*%K.y =
  # < k_x, k_y >_{H(T_0)}. Meaning that Cov.inv = K_{T_0}_^{-1}. K_{T_0}
  # is obtained via projection of K onto H_0^\perp
  
}  # End k.post


# I bet this corresponds to eq (13) in the RKHS text, assuming prior mean = 0
mean.post=function(X, g) {
  
  # X defines the NEW locations (away from the known set)
  # In computing the inner product, we are only evaluating g on the known
  # subset
  
  n.x= dim(X)
  K.g= matrix(0,nrow=1, ncol=4*length(basis))
  K.e= matrix(0,nrow=n.x, ncol=4*length(basis))
  
  
  for (i in 1:4) {
    
    #k.x=k(X,xg[[paste(i)]])
    k.x= K2(X, xg[[paste(i)]])
    g.x= matrix(apply(xg[[paste(i)]], 1, g), nrow=1)  # This is where g is
    # evaluated on the known subset
    K.e[,((i-1)*length(basis)+1):(i*length(basis))]= 
      k.x%*%diag(quad[[2]])%*%t(V)
    K.g[1,((i-1)*length(basis)+1):(i*length(basis))]=
      g.x%*%diag(quad[[2]])%*%t(V);
    
  }
  
  return(K.e%*%Cov.inv%*%t(K.g))
  
}  # End mean.post





################# Test with training data ######################################
################################################################################
################################################################################

# We are going to loop over training data size, but keeping the number of
# psuedo-observations the same for ordinary kriging

# corner peak fn

#Ns <- c(20, 40, 60, 80, 100)  # Number of interior observations
Ns <- seq(10, 200, by= 10)  # Number of interior observations


pk.err <- rep(0, length(Ns))
ps.krig <- rep(0, length(Ns))
ps.k.reg <- rep(0, length(Ns))

for(sz.ind in 1:length(Ns)){
  
  N <- Ns[sz.ind]  # Number of interior observations = 10*d
  D <- 2*randomLHS(n=N, k=2)-1  # LHS in 2D, converting [0, 1] -> [-1, 1]
  
  reg.bayes <- seq(0, 8, length.out= 4*length(basis))  # Sequence of points
  # to map to the edges of [-1, 1]
  n.reg <- length(reg.bayes)
  reg.vals <- t(sapply(reg.bayes,l))  # Convert to 2D boundary points
  D.reg <- rbind(D, reg.vals)  # Combine boundary points with interior points. 
  # These are the training points for the usual kriging predictor - analagous
  # the 4*L nodes around the boundary that were used to compute the spectral
  # projection
  
  
  
  # t <- seq(0, 8, 0.1)
  # test.b <- t(sapply(t,l)) # Again, comes back around to repeat [-1, 1]
  # n.test <- dim(test.b)[1]
  # # Pushing away from the boundary
  # test <- rbind(0.9*test.b, 0.5*test.b)
  
  # For this example, we want to approximate the function over the entire grid
  x.vals <- seq(-1, 1, length= 30)
  y.vals <- seq(-1, 1, length= 30)
  test.grid <- as.matrix(expand.grid(x.vals, y.vals))
  
  z <- apply(test.grid, 1, f.corn)  # True values
  
  
  y <- apply(D, 1, f.corn) #+ rnorm(N, mean=0, sd=0.05)  # Just the interior points
  y.reg <- apply(D.reg, 1, f.corn) #+ rnorm((N + n.reg), mean=0, sd=0.05)  # Interior and
  # boundary observations
  
  
  
  
  
  ########## pkGP ##############
  
  mu.post <- mean.post(test.grid, f.corn)  # Computing eq. (13) for t = test points. 
  # This, and eq (15) determine the prior GP: GP(\mu_0, k_0)
  
  # mu.post and K.post determine the (projected) GP "PRIOR" distribution. Which
  # is updated to obtain the posterior predictive distribuiton. Below is the 
  # mean of this distribution == the linear predictor using pGP
  # post.mean <- mu.post + 
  #     K.post(test, D)%*%solve(K.post(D, D) + 0.05*diag(dim(D)[1]))%*%
  #     (y-mean.post(D,f))
  post.mean <- mu.post + K.post(test.grid, D)%*%solve(K.post(D, D))%*%
    (y-mean.post(D, f.corn)) # We only need the k_0(.,.) at the interior points
  # because the boundary has already been absorbed
  # into K.post (the basis functions and the i.p.
  # on H(T_0)
  
  
  
  
  ########### ordinary kriging ############
  
  # The typical predictor (no pseduo observations)
  # E.reg=eigen(k(D.reg,D.reg))#+1e-10*diag(N+n.reg));
  E.reg <- eigen(K2(D, D) + 1e-6*diag(dim(D)[1]))
  # E.reg <- eigen(K2(D.reg, D.reg))
  K.reg.inv <- E.reg$vectors%*%diag(1/abs(E.reg$values))%*%t(E.reg$vectors)
  #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
  post.mean.reg <- K2(test.grid, D)%*%K.reg.inv%*%y  # The usual predictor
  
  
  # Pseudo with additional training data
  E.reg.nug <- eigen(K2(D.reg, D.reg) + 1e-6*diag(dim(D.reg)[1]))
  #E.reg <- eigen(K2(D.reg, D.reg))
  K.reg.inv.nug <- E.reg.nug$vectors%*%diag(1/abs(E.reg.nug$values))%*%
    t(E.reg.nug$vectors)
  #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
  post.mean.reg.nug <- K2(test.grid, D.reg)%*%K.reg.inv.nug%*%y.reg  # The usual pred.
  
  
  
  
  
  
  f.test <- apply(test.grid, 1, f.corn)  # True values
  
  
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
  pk.err[sz.ind] <- max(abs(f.test-post.mean))/max(abs(f.test))
  ps.krig[sz.ind] <- max(abs(f.test-post.mean.reg.nug))/max(abs(f.test))
  ps.k.reg[sz.ind] <- max(abs(f.test-post.mean.reg))/max(abs(f.test))
  
  # Compute relative L1 errors
  # pk.err[sz.ind] <- norm(f.test - post.mean, type= "o")/
  #   norm(as.matrix(f.test), type= "o")
  # ps.krig[sz.ind] <- norm(f.test - post.mean.reg.nug, type= "o")/
  #   norm(as.matrix(f.test), type= "o")
  
  
}  # End loop over training data size







################################################################################
################################################################################
################################################################################

# Repeat the experiment above with the product peak function


# product peak fn

#Ns <- c(20, 40, 60, 80, 100)  # Number of interior observations
Ns <- seq(10, 200, by= 10)  # Number of interior observations


pk.err.prod <- rep(0, length(Ns))
ps.krig.prod <- rep(0, length(Ns))
ps.k.reg.prod <- rep(0, length(Ns))

for(sz.ind in 1:length(Ns)){
  
  N <- Ns[sz.ind]  # Number of interior observations = 10*d
  D <- 2*randomLHS(n=N, k=2)-1  # LHS in 2D, converting [0, 1] -> [-1, 1]
  
  reg.bayes <- seq(0, 8, length.out= 4*length(basis))  # Sequence of points
  # to map to the edges of [-1, 1]
  n.reg <- length(reg.bayes)
  reg.vals <- t(sapply(reg.bayes,l))  # Convert to 2D boundary points
  D.reg <- rbind(D, reg.vals)  # Combine boundary points with interior points. 
  # These are the training points for the usual kriging predictor - analagous
  # the 4*L nodes around the boundary that were used to compute the spectral
  # projection
  
  
  
  # t <- seq(0, 8, 0.1)
  # test.b <- t(sapply(t,l)) # Again, comes back around to repeat [-1, 1]
  # n.test <- dim(test.b)[1]
  # # Pushing away from the boundary
  # test <- rbind(0.9*test.b, 0.5*test.b)
  
  # For this example, we want to approximate the function over the entire grid
  x.vals <- seq(-1, 1, length= 30)
  y.vals <- seq(-1, 1, length= 30)
  test.grid <- as.matrix(expand.grid(x.vals, y.vals))
  
  #z <- apply(test.grid, 1, f.corn)  # True values
  
  
  y <- apply(D, 1, f.prod) #+ rnorm(N, mean=0, sd=0.05)  # Just the interior points
  y.reg <- apply(D.reg, 1, f.prod) #+ rnorm((N + n.reg), mean=0, sd=0.05)  # Interior and
  # boundary observations
  
  
  
  
  
  ########## pkGP ##############
  
  mu.post <- mean.post(test.grid, f.prod)  # Computing eq. (13) for t = test points. 
  # This, and eq (15) determine the prior GP: GP(\mu_0, k_0)
  
  # mu.post and K.post determine the (projected) GP "PRIOR" distribution. Which
  # is updated to obtain the posterior predictive distribuiton. Below is the 
  # mean of this distribution == the linear predictor using pGP
  # post.mean <- mu.post + 
  #     K.post(test, D)%*%solve(K.post(D, D) + 0.05*diag(dim(D)[1]))%*%
  #     (y-mean.post(D,f))
  post.mean <- mu.post + K.post(test.grid, D)%*%solve(K.post(D, D))%*%
    (y-mean.post(D, f.prod)) # We only need the k_0(.,.) at the interior points
  # because the boundary has already been absorbed
  # into K.post (the basis functions and the i.p.
  # on H(T_0)
  
  
  
  
  ########### ordinary kriging ############
  
  # The typical predictor (no pseduo observations)
  # E.reg=eigen(k(D.reg,D.reg))#+1e-10*diag(N+n.reg));
  E.reg <- eigen(K2(D, D) + 1e-6*diag(dim(D)[1]))
  # E.reg <- eigen(K2(D.reg, D.reg))
  K.reg.inv <- E.reg$vectors%*%diag(1/abs(E.reg$values))%*%t(E.reg$vectors)
  #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
  post.mean.reg <- K2(test.grid, D)%*%K.reg.inv%*%y  # The usual predictor
  
  
  # Pseudo with additional training data
  E.reg.nug <- eigen(K2(D.reg, D.reg) + 1e-6*diag(dim(D.reg)[1]))
  #E.reg <- eigen(K2(D.reg, D.reg))
  K.reg.inv.nug <- E.reg.nug$vectors%*%diag(1/abs(E.reg.nug$values))%*%
    t(E.reg.nug$vectors)
  #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
  post.mean.reg.nug <- K2(test.grid, D.reg)%*%K.reg.inv.nug%*%y.reg  # The usual pred.
  
  
  
  
  
  
  
  f.test <- apply(test.grid, 1, f.prod)  # True values
  
  
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
  pk.err.prod[sz.ind] <- max(abs(f.test-post.mean))/max(abs(f.test))
  ps.krig.prod[sz.ind] <- max(abs(f.test-post.mean.reg.nug))/max(abs(f.test))
  ps.k.reg.prod[sz.ind] <- max(abs(f.test-post.mean.reg))/max(abs(f.test))
  
 
  
  
}  # End loop over training data size





################################################################################
################################################################################
################################################################################

# Repeat the experiment above with the product peak function


#Ns <- c(20, 40, 60, 80, 100)  # Number of interior observations
Ns <- seq(10, 200, by= 10)  # Number of interior observations


pk.err.rosen <- rep(0, length(Ns))
ps.krig.rosen <- rep(0, length(Ns))
ps.k.reg.rosen <- rep(0, length(Ns))

for(sz.ind in 1:length(Ns)){
  
  N <- Ns[sz.ind]  # Number of interior observations = 10*d
  D <- 2*randomLHS(n=N, k=2)-1  # LHS in 2D, converting [0, 1] -> [-1, 1]
  
  reg.bayes <- seq(0, 8, length.out= 4*length(basis))  # Sequence of points
  # to map to the edges of [-1, 1]
  n.reg <- length(reg.bayes)
  reg.vals <- t(sapply(reg.bayes,l))  # Convert to 2D boundary points
  D.reg <- rbind(D, reg.vals)  # Combine boundary points with interior points. 
  # These are the training points for the usual kriging predictor - analagous
  # the 4*L nodes around the boundary that were used to compute the spectral
  # projection
  
  
  
  # t <- seq(0, 8, 0.1)
  # test.b <- t(sapply(t,l)) # Again, comes back around to repeat [-1, 1]
  # n.test <- dim(test.b)[1]
  # # Pushing away from the boundary
  # test <- rbind(0.9*test.b, 0.5*test.b)
  
  # For this example, we want to approximate the function over the entire grid
  x.vals <- seq(-1, 1, length= 30)
  y.vals <- seq(-1, 1, length= 30)
  test.grid <- as.matrix(expand.grid(x.vals, y.vals))
  
  #z <- apply(test.grid, 1, f.corn)  # True values
  
  
  y <- apply(D, 1, f.rosen) #+ rnorm(N, mean=0, sd=0.05)  # Just the interior points
  y.reg <- apply(D.reg, 1, f.rosen) #+ rnorm((N + n.reg), mean=0, sd=0.05)  # Interior and
  # boundary observations
  
  
  
  
  
  ########## pkGP ##############
  
  mu.post <- mean.post(test.grid, f.rosen)  # Computing eq. (13) for t = test points. 
  # This, and eq (15) determine the prior GP: GP(\mu_0, k_0)
  
  # mu.post and K.post determine the (projected) GP "PRIOR" distribution. Which
  # is updated to obtain the posterior predictive distribuiton. Below is the 
  # mean of this distribution == the linear predictor using pGP
  # post.mean <- mu.post + 
  #     K.post(test, D)%*%solve(K.post(D, D) + 0.05*diag(dim(D)[1]))%*%
  #     (y-mean.post(D,f))
  post.mean <- mu.post + K.post(test.grid, D)%*%solve(K.post(D, D))%*%
    (y-mean.post(D, f.rosen)) # We only need the k_0(.,.) at the interior points
  # because the boundary has already been absorbed
  # into K.post (the basis functions and the i.p.
  # on H(T_0)
  
  
  
  
  ########### ordinary kriging ############
  
  # The typical predictor (no pseduo observations)
  # E.reg=eigen(k(D.reg,D.reg))#+1e-10*diag(N+n.reg));
  E.reg <- eigen(K2(D, D) + 1e-6*diag(dim(D)[1]))
  # E.reg <- eigen(K2(D.reg, D.reg))
  K.reg.inv <- E.reg$vectors%*%diag(1/abs(E.reg$values))%*%t(E.reg$vectors)
  #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
  post.mean.reg <- K2(test.grid, D)%*%K.reg.inv%*%y  # The usual predictor
  
  
  # Pseudo with additional training data
  E.reg.nug <- eigen(K2(D.reg, D.reg) + 1e-6*diag(dim(D.reg)[1]))
  #E.reg <- eigen(K2(D.reg, D.reg))
  K.reg.inv.nug <- E.reg.nug$vectors%*%diag(1/abs(E.reg.nug$values))%*%
    t(E.reg.nug$vectors)
  #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
  post.mean.reg.nug <- K2(test.grid, D.reg)%*%K.reg.inv.nug%*%y.reg  # The usual pred.
  

  
  
  f.test <- apply(test.grid, 1, f.rosen)  # True values
  
  
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
  pk.err.rosen[sz.ind] <- max(abs(f.test-post.mean))/max(abs(f.test))
  ps.krig.rosen[sz.ind] <- max(abs(f.test-post.mean.reg.nug))/max(abs(f.test))
  ps.k.reg.rosen[sz.ind] <- max(abs(f.test-post.mean.reg))/max(abs(f.test))
  
  # Compute relative L1 errors
  # pk.err.rosen[sz.ind] <- norm(f.test - post.mean, type= "o")/
  #   norm(as.matrix(f.test), type= "o")
  # ps.krig.rosen[sz.ind] <- norm(f.test - post.mean.reg.nug, type= "o")/
  #   norm(as.matrix(f.test), type= "o")
  
  
}  # End loop over training data size








################################################################################
################################################################################
################################################################################

#### Plot the three functions, then plot relative errors vs. training data size
x.vals <- seq(-1, 1, length= 50)
y.vals <- seq(-1, 1, length= 50)
grid <- expand.grid(x.vals, y.vals)

z <- apply(grid, 1, f.corn)  # True values

x.mat <- matrix(rep(x.vals, times= length(y.vals)), ncol= length(x.vals), byrow= TRUE)
y.mat <- matrix(rep(y.vals, times= length(x.vals)), ncol= length(y.vals), byrow= TRUE)
z.mat <- matrix(z, ncol= length(x.vals), nrow= length(y.vals), byrow= TRUE)


# Figure 5
pdf("threeFunctionsPlot.pdf", width= 60, height= 15)
par(mfrow= c(1,3))
persp3D(x.vals, y.vals, z.mat, 
        colvar = z.mat,
        col = viridis(100),
        alpha = 0.7,
        theta = 40, phi = 30, expand = 0.7,
        xlab = "", ylab = "", zlab = "", main = "",
        contour = TRUE, colkey = FALSE,
        lighting = TRUE, ticktype = "detailed", axes = TRUE,
        cex.lab = 3,    # Increases size of x, y, z labels
        cex.axis = 2.6,   # Increases size of axis tick labels
        cex.main = 4)



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



z <- apply(grid, 1, f.prod)  # True values

x.mat <- matrix(rep(x.vals, times= length(y.vals)), ncol= length(x.vals), byrow= TRUE)
y.mat <- matrix(rep(y.vals, times= length(x.vals)), ncol= length(y.vals), byrow= TRUE)
z.mat <- matrix(z, ncol= length(x.vals), nrow= length(y.vals), byrow= TRUE)


#x11(width= 45, height= 45)
#pdf("boundaryExPerspPlot.pdf", width= 25, height= 25)
persp3D(x.vals, y.vals, z.mat, 
        colvar = z.mat,
        col = viridis(100),
        alpha = 0.7,
        theta = 40, phi = 30, expand = 0.7,
        xlab = "", ylab = "", zlab = "", main = "",
        contour = TRUE, colkey = FALSE,
        lighting = TRUE, ticktype = "detailed", axes = TRUE,
        cex.lab = 3,    # Increases size of x, y, z labels
        cex.axis = 2.6,   # Increases size of axis tick labels
        cex.main = 4)



z <- apply(grid, 1, f.rosen)  # True values

x.mat <- matrix(rep(x.vals, times= length(y.vals)), ncol= length(x.vals), byrow= TRUE)
y.mat <- matrix(rep(y.vals, times= length(x.vals)), ncol= length(y.vals), byrow= TRUE)
z.mat <- matrix(z, ncol= length(x.vals), nrow= length(y.vals), byrow= TRUE)


#x11(width= 45, height= 45)
#pdf("boundaryExPerspPlot.pdf", width= 25, height= 25)
persp3D(x.vals, y.vals, z.mat, 
        colvar = z.mat,
        col = viridis(100),
        alpha = 0.7,
        theta = 40, phi = 30, expand = 0.7,
        xlab = "", ylab = "", zlab = "", main = "",
        contour = TRUE, colkey = FALSE,
        lighting = TRUE, ticktype = "detailed", axes = TRUE,
        cex.lab = 3,    # Increases size of x, y, z labels
        cex.axis = 2.6,   # Increases size of axis tick labels
        cex.main = 4)
dev.off()


#################################### Plot the relative errors

# Figure 6
pdf("threeSupErrsPlot.pdf", width= 60, height= 15)
# corner peak
par(mfrow= c(1,3), mar = c(5, 5, 4, 2) + 0.1)
plot(Ns, log(ps.krig), type= "b", pch= 8, cex= 2, col= "cyan",
     xlab= "Training data size", ylab= "Relative sup error (log scale)",
     main= "Corner Peak", cex.main= 4,
     ylim= c(-10, 1.5), cex.lab= 2.5, cex.axis= 2.5,
     lwd= 3)



lines(Ns, log(ps.k.reg), type= "b", pch= 14, cex= 2, col= "blue",
      lwd= 3)

lines(Ns, log(pk.err), type= "b", pch= 25, cex= 2, col= "red",
      lwd= 3)


legend("topright", legend= c("pkGP", 
                             #"bdryGP", 
                             "Kriging", 
                             "Pseudo-Kriging"),
       col= c("red", 
              #"green", 
              "blue", 
              "cyan"), 
       pch= c(8, 
              #13, 
              14, 
              25), cex= 2.5, bty= "n",
       lwd= 3)



# product peak
plot(Ns, log(ps.krig.prod), type= "b", pch= 8, cex= 2, col= "cyan",
     xlab= "Training data size", ylab= "",
     main= "Product Peak", cex.main= 4,
     ylim= c(-10, 1.5), cex.lab= 2.5, cex.axis= 2.5,
     lwd= 3)



lines(Ns, log(ps.k.reg.prod), type= "b", pch= 14, cex= 2, col= "blue",
      lwd= 3)

lines(Ns, log(pk.err.prod), type= "b", pch= 25, cex= 2, col= "red",
      lwd= 3)

legend("topright", legend= c("pkGP", 
                             #"bdryGP", 
                             "Kriging", 
                             "Pseudo-Kriging"),
       col= c("red", 
              #"green", 
              "blue", 
              "cyan"), 
       pch= c(8, 
              #13, 
              14, 
              25), cex= 2.5, bty= "n",
       lwd= 3)



# Rosenbrock
plot(Ns, log(pk.err.rosen), type= "b", pch= 8, cex= 2, col= "red",
     xlab= "Training data size", ylab= "",
     main= "Rosenbrock", cex.main= 4,
     ylim= c(-10, 1.5), cex.lab= 2.5, cex.axis= 2.5,
     lwd= 3)



lines(Ns, log(ps.k.reg.rosen), type= "b", pch= 14, cex= 2, col= "blue",
      lwd= 3)

lines(Ns, log(ps.krig.rosen), type= "b", pch= 25, cex= 2, col= "cyan",
      lwd= 3)

legend("topright", legend= c("pkGP", 
                             #"bdryGP", 
                             "Kriging", 
                             "Pseudo-Kriging"),
       col= c("red", 
              #"green", 
              "blue", 
              "cyan"), 
       pch= c(8, 
              #13, 
              14, 
              25), cex= 2.5, bty= "n",
       lwd= 3)
dev.off()

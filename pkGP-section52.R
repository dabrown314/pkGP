# pkGP-section52.R
# Andrew Brown and John Nicholson
# Sept. 19, 2025
#
# This script reproduces the simulation study with diagonal conditions, 
# presented in subection 5.2 of the pkGP JUQ paper.


rm(list= ls())

library(orthopolynom)
library(statmod)
library(lhs)
library(plot3D)
library(viridis)
sessionInfo()

set.seed(420)



################################################################################

# Define the test function
f <- function(x) {
  
  result= x[2]*sqrt(1 + x[1])*cos(pi*x[2])*sin(0.5*pi*(x[1] - x[2]) + 1)*
    exp(0.5*(x[1] + x[2])^2)
  
  return(result)
  
}




################################################################################

J <- 15  # Max order of polynomial approximation for Rayleigh-Ritz method

# This will hold the sup norm relative errors and RMSE's associated with 
# Legendre polyn. approximations of varying orders, given by J
leg.max <- rep(0, J-4)  
leg.rmse <- leg.max
#leg.rmse=leg.max=rep(0,J-4);


params <- c(1, 1, 1)  # (scale param., length-scale in x, length-scale in y)


# Will bury this in a loop later, looping over polyn. order m
m <- 15


# Matern kernel w/ nu = 3/2, correlation length = sqrt(params[2])
K2 <- function(X, Y){
  # X \in R^2, Y \in R^2
  
  dist.mat= as.matrix(dist(rbind(X, Y)))
  pairwise_dist= dist.mat[1:nrow(X), (nrow(X)+1):(nrow(X) + nrow(Y))]
  len.par= sqrt(params[2])
  
  out= params[1]*(1 + (pairwise_dist*sqrt(3)/len.par))*
    exp(-1*pairwise_dist*sqrt(3)/len.par)
  
  return(out)
  
}  # End fn K2


# Get nodes/weights for Gaussian quadrature
quad <- gauss.quad(4*m, kind= "legendre")
# quad <- gauss.quad(4*m, kind= "chebyshev1")  
n.nodes <- length(quad[[1]])


#list of legendre polynomials
leg <- lapply(legendre.polynomials(m, normalized=TRUE), as.function)
basis <- leg


# To go across the diagonal, it is like the boundary example, but we
# only need to go along one "edge" instead of around 4 edges

# # xg = 'x grid'
xg <- list()
xg[["1"]] <- cbind(quad[[1]], quad[[1]])  # Diagonal is of the form (x_i, x_i),
# x_i \in [-1, 1]


# Use lists of length one to simplify things right now
K.bas=list() 
K.bas[["11"]] <- K2(xg[["1"]], xg[["1"]])  # K(x_i, x_i) for i= 1, ..., n.nodes



# Now evaluate the basis functions (linear combinations of which are used to 
# approximate the e-vects of K.bas) at the nodes. Only want the nodes along the'
# diagonal
L <- length(basis)
phi <- matrix(0, ncol= L, nrow= length(quad[[1]]))
for (i in 1:L) {
  
  phi[, i] <- sapply(quad[[1]], basis[[i]])  # Evaluate the basis functions at
  # the nodes
  
}


K.phi <- K.bas[['11']]%*%diag(quad[[2]])%*%phi
K.int <- t(phi)%*%diag(quad[[2]])%*%K.phi  # \approx A matrix of Oya et al.,
# along the diagonal.



E <- eigen(K.int)  # Solves the (ordinary) e-val problem in eq (5) of 
# Oya et al., where K.int = A and C = identity since basis fns are
# orthonormal. E$vectors contains the e-vectors as COLUMN vectors.

###contains eigenvectors of integral operator evaluated at quad[[1]]
V <- t(E$vectors)%*%t(phi)  # jth row of V = discretized appx e-fn j of the
# integral operator, evaluated at x_1, ..., x_n. These guys are orthogonal,
# but not normalized.


# To (hopefully) carry over from the boundary example, structure cov.mat the
# same way, but with only 1 edge rather than four
Cov <- matrix(0, nrow= length(basis), ncol= length(basis));
sparse <- 0
for (i in 1:1) {
  
  for (j in 1:1) {
    
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
      
      
      # This is the VD_wKD_wV' expression.
      #
      # Here, V is not changing w/ i,j, but K.bas is??
      Cov[((i-1)*length(basis)+1):(i*length(basis)),
          ((j-1)*length(basis)+1):(j*length(basis))] <- V%*%diag(quad[[2]])%*%
        K.bas[[name]]%*%diag(quad[[2]])%*%t(V)
      # The e-vects (E) are computed from K.bas[["11"]] above
      
      
      
    }  # End ifelse
  }  # End j loop
}  # End i loop



E.cov <- eigen(Cov + 1e-6*diag(length(basis)), symmetric=TRUE)  # Then
# using SVD to invert

Cov.inv <- E.cov$vectors%*%diag(1/abs(E.cov$values))%*%t(E.cov$vectors)




# Use cov.inv to define the updated GP prior model, k.post and mean.post

# This function is the (rough) approximation of, and corresponds to, eq. (15)
# in the RKHS text - the updated kernel. Again adjusting the arguments to use
# 1 "edge" instead of 4
K.post <- function(X,Y) {
  
  n.x= dim(X)[1]
  n.y= dim(Y)[1]
  K.x= matrix(0, nrow=n.x, ncol= length(basis))
  K.y= matrix(0, nrow=n.y, ncol= length(basis))
  
  for (i in 1:1) {
    
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
  
  n.x= dim(X)
  K.g= matrix(0,nrow=1, ncol= length(basis))
  K.e= matrix(0,nrow=n.x, ncol= length(basis))
  
  
  for (i in 1:1) {
    
    #k.x=k(X,xg[[paste(i)]])
    k.x= K2(X, xg[[paste(i)]])
    g.x= matrix(apply(xg[[paste(i)]], 1, g), nrow=1)
    K.e[,((i-1)*length(basis)+1):(i*length(basis))]= 
      k.x%*%diag(quad[[2]])%*%t(V)
    K.g[1,((i-1)*length(basis)+1):(i*length(basis))]=
      g.x%*%diag(quad[[2]])%*%t(V);
    
  }
  
  return(K.e%*%Cov.inv%*%t(K.g))
  
}  # End mean.post



################# Putting it to work ###########################################

N <- 20  # Number of interior observations = 10*d
D <- 2*randomLHS(n=N, k=2)-1  # LHS in 2D, converting [0, 1] -> [-1, 1]



reg.bayes <- seq(-1, 1, length.out= length(basis))# Sequence of points
reg.vals <- cbind(reg.bayes, reg.bayes)  # Convert to 2D boundary points
D.reg <- rbind(D, reg.vals)  # Combine boundary points with interior points. 
# These are the training points for the usual kriging predictor - analagous
# the 4*L nodes around the boundary that were used to compute the spectral
# projection



t <- seq(-0.9, 0.9, 0.1)
#test.b <- t(sapply(t,l)) # Again, comes back around to repeat [-1, 1]
test.b <- cbind(t, t)  # Convert to 2D boundary points
n.test <- dim(test.b)[1] 
# Pushing away from the boundary
test <- rbind(cbind(t, t-0.1), cbind(t, t+0.1))



y <- apply(D, 1, f) #+ rnorm(N, mean=0, sd=0.05)  # Just the interior points
y.reg <- apply(D.reg, 1, f) #+ rnorm((N + n.reg), mean=0, sd=0.05)  # Interior and
# boundary observations



########## pkGP ##############

mu.post <- mean.post(test, f)  # Computing eq. (13) for t = test points. This, 
# and eq (15) determine the prior GP: GP(\mu_0, k_0)

# mu.post and K.post determine the (projected) GP "PRIOR" distribution. Which
# is updated to obtain the posterior predictive distribuiton. Below is the 
# mean of this distribution == the linear predictor using pGP
# post.mean <- mu.post + 
#     K.post(test, D)%*%solve(K.post(D, D) + 0.05*diag(dim(D)[1]))%*%
#     (y-mean.post(D,f))

post.mean <- mu.post + K.post(test, D)%*%solve(K.post(D, D))%*%
  (y - mean.post(D,f)) # We only need the k_0(.,.) at the interior points
# because the boundary has already been absorbed
# into K.post (the basis functions and the i.p.
# on H(T_0)




########### ordinary kriging ############

# The typical predictor by conditioning on all points, including those
# at the boundary
#E.reg=eigen(k(D.reg,D.reg))#+1e-10*diag(N+n.reg));
#E.reg <- eigen(K2(D.reg, D.reg) + 0.05*diag(dim(D.reg)[1]))
# E.reg <- eigen(K2(D.reg, D.reg))
# K.reg.inv <- E.reg$vectors%*%diag(1/abs(E.reg$values))%*%t(E.reg$vectors)
# #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
# post.mean.reg <- K2(test,D.reg)%*%K.reg.inv%*%y.reg  # The usual predictor


# Pseudo-kriging 
E.reg.nug <- eigen(K2(D.reg, D.reg) + 1e-6*diag(dim(D.reg)[1]))
#E.reg <- eigen(K2(D.reg, D.reg))
K.reg.inv.nug <- E.reg.nug$vectors%*%diag(1/abs(E.reg.nug$values))%*%
  t(E.reg.nug$vectors)
#post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
post.mean.reg.nug <- K2(test, D.reg)%*%K.reg.inv.nug%*%y.reg  # The usual pred.




f.test <- apply(test, 1, f)  # True values



# Metrics
cor(f.test, post.mean)
cor(f.test, post.mean.reg.nug)


# Relative errors
max(abs(f.test-post.mean))/max(f.test)
max(abs(f.test-post.mean.reg.nug))/max(f.test)


sqrt(sum((f.test-post.mean)^2)/dim(test)[1])
sqrt(sum((f.test-post.mean.reg.nug)^2)/dim(test)[1])


# Relative L2 error
sqrt(sum(((f.test-post.mean)^2)))/sqrt(sum(f.test^2))
sqrt(sum(((f.test-post.mean.reg.nug)^2)))/sqrt(sum(f.test^2))


###############################################################################
############################# Visuals           ###############################


# Figure 8
pdf("diagonalExPredsVsTrue.pdf", width= 30, height= 15)
par(mfrow= c(1,2), mar = c(5, 5, 4, 2) + 0.1)
plot(post.mean, f.test, xlab= "Prediction", ylab= "True f", 
     pch= 20, cex= 3, main= "pkGP",
     cex.lab = 3,    # Increases size of x, y, z labels
     cex.axis = 2.6,   # Increases size of axis tick labels
     cex.main = 4)
abline(a= 0, b= 1, lwd= 2)
legend("topleft", 
       legend= "RMSE = 0.1413", cex= 4, bty= "n")

plot(post.mean.reg.nug, f.test, pch= 20, xlab= "Prediction", cex= 3,
     ylab= "", main= "Pseudo-Kriging",
     cex.lab = 3,    # Increases size of x, y, z labels
     cex.axis = 2.6,   # Increases size of axis tick labels
     cex.main = 4)
abline(a= 0, b= 1, lwd= 2)
legend("topleft", 
       legend= "RMSE = 0.2240", cex= 4, bty= "n")
dev.off()





#### Plot
x.vals <- seq(-1, 1, length= 50)
y.vals <- seq(-1, 1, length= 50)
grid <- expand.grid(x.vals, y.vals)

z <- apply(grid, 1, f)

x.mat <- matrix(rep(x.vals, times= length(y.vals)), ncol= length(x.vals), 
                byrow= TRUE)
y.mat <- matrix(rep(y.vals, times= length(x.vals)), ncol= length(y.vals), 
                byrow= TRUE)
z.mat <- matrix(z, ncol= length(x.vals), nrow= length(y.vals), byrow= TRUE)



# Figure 7
pdf("diagonalExPerspPlot.pdf", width= 25, height= 25)
persp3D(x.vals, y.vals, z.mat, 
        colvar = z.mat,
        col = viridis(100),
        alpha = 0.7,
        theta = 40, phi = 30, expand = 0.7,
        xlab = "x", ylab = "y", zlab = "f(x, y)", main = "",
        contour = TRUE, colkey = FALSE,
        lighting = TRUE, ticktype = "detailed", axes = TRUE,
        cex.lab = 3,    # Increases size of x, y, z labels
        cex.axis = 2.6,   # Increases size of axis tick labels
        cex.main = 4)



points3D(test[, 2], test[, 1], post.mean.reg.nug, add= TRUE, pch= 25, cex= 2.0, 
         col= "blue")


points3D(test[, 2], test[, 1], post.mean, add= TRUE, pch= 8, cex= 2.0, 
         col= "red")


points3D(D.reg[, 2], D.reg[, 1], y.reg, add= TRUE, pch= 20, cex= 2.0, 
         col= "black")


legend("topleft", legend= c("pkGP", "Pseudo-Kriging"),
       col= c("red", "blue"), 
       pch= c(8, 25), cex= 4, bty= "n")
dev.off()




################################################################################
################################################################################


# Do ordinary kriging with an increasing number of pseudo-observations and plot
# results vs. pkGP

# No. of pseudo-observations is in reg.vals

r.errs <- rep(0, 10)
r.errs.full <- rep(0, 10)
n.ps <- seq(6, 55, length= 10)  # Number of pseudo-observations


x.vals <- seq(-1, 1, length= 30)
y.vals <- seq(-1, 1, length= 30)
test.grid <- as.matrix(expand.grid(x.vals, y.vals))
f.test.grid <- apply(test.grid, 1, f)  # True values

for (i in 1:10) {
  
  reg.bayes.it <- seq(-1, 1, length.out= n.ps[i])# Sequence of points
  reg.vals.it <- cbind(reg.bayes.it, reg.bayes.it)  # Convert to 2D boundary points
  D.reg.it <- rbind(D, reg.vals.it)
  
  y.reg.it <- apply(D.reg.it, 1, f)
  
  
  # E.reg <- eigen(K2(D.reg, D.reg))
  # K.reg.inv <- E.reg$vectors%*%diag(1/abs(E.reg$values))%*%t(E.reg$vectors)
  # #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
  # post.mean.reg <- K2(test,D.reg)%*%K.reg.inv%*%y.reg  # The usual predictor
  
  
  # Try with and without nugget
  E.reg.nug <- eigen(K2(D.reg.it, D.reg.it) + 1e-6*diag(dim(D.reg.it)[1]))
  #E.reg <- eigen(K2(D.reg, D.reg))
  K.reg.inv.nug <- E.reg.nug$vectors%*%diag(1/abs(E.reg.nug$values))%*%
    t(E.reg.nug$vectors)
  #post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
  post.mean.reg.nug <- K2(test, D.reg.it)%*%K.reg.inv.nug%*%y.reg.it
  
  post.mean.reg.full <- K2(test.grid, D.reg.it)%*%K.reg.inv.nug%*%y.reg.it
  
  
  #r.errs[i] <- max(abs(f.test-post.mean.reg.nug))/max(f.test)
  r.errs[i] <- sqrt(sum(((f.test-post.mean.reg.nug)^2)))/sqrt(sum(f.test^2))
  #r.errs.full[i] <- max(abs(f.test.grid-post.mean.reg.full))/max(f.test.grid)
  r.errs.full[i] <- sqrt(sum(((f.test.grid-post.mean.reg.full)^2)))/sqrt(sum(f.test.grid^2))
  
}



mu.post.grid <- mean.post(test.grid, f)  # Computing eq. (13) for t = test points. This, 
# and eq (15) determine the prior GP: GP(\mu_0, k_0)

# mu.post and K.post determine the (projected) GP "PRIOR" distribution. Which
# is updated to obtain the posterior predictive distribuiton. Below is the 
# mean of this distribution == the linear predictor using pGP
# post.mean <- mu.post + 
#     K.post(test, D)%*%solve(K.post(D, D) + 0.05*diag(dim(D)[1]))%*%
#     (y-mean.post(D,f))

post.mean.grid <- mu.post.grid + K.post(test.grid, D)%*%solve(K.post(D, D))%*%
  (y - mean.post(D,f)) # We only need the k_0(.,.) at the interior points
# because the boundary has already been absorbed
# into K.post (the basis functions and the i.p.
# on H(T_0)



# Figure 9
pdf("diagonalExPseudos.pdf", width= 25, height= 25)
#par(mfrow= c(1, 2), mar = c(5, 5, 4, 2) + 0.1)
par(mar = c(5, 5, 4, 2) + 0.1)



plot(n.ps, r.errs.full, type= "b", pch= 20, cex= 3, ylim= c(0, 0.8), lwd= 3,
     xlab= "Number of Pseudo-observations", ylab= "Relative L2 Error",
     main= "", 
     cex.lab = 3,    # Increases size of x, y, z labels
     cex.axis = 2.6,   # Increases size of axis tick labels
     cex.main = 4)
abline(h= sqrt(sum(((f.test.grid-post.mean.grid)^2)))/sqrt(sum(f.test.grid^2)), lty= 2, lwd= 3)
legend("topright",
       legend= c("pkGP", "Pseudo-Kriging"),
       lty= c(2,1), lwd= 3, bty= "n", cex= 4)
dev.off()
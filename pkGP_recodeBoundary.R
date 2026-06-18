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

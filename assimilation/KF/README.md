# lorenz96.jl
EKF data assimilation system for the nonlinear model in Lorenz (1996)

## Basic equations of the Lorenz96 model
\\[a+b=C\\]

# lorenz96_hybrid.jl
Hybrid EKF data assimilation system for the nonlinear model in Lorenz (1996)

## Required modules
lorenz96_module.jl

## Required files
You need to run lorenz96_NMC.jl to create a static background covariance error matrix file (Pstat.bin). 

# lorenz96_NMC.jl
Script to create a binary file of Pstat.bin (provides a static background covariance error matrix based on the NMC method). 

# lorenz96_module.jl
Module of functions to solve the nonlinear model in Lorenz (1996)


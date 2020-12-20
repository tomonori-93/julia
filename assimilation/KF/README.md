[**Formulation**](https://tomonori-93.github.io/julia/assimilation/KF/formulation.md)

# EnKF_module.jl
Module script: Basic functions for EnKF (Ensemble Kalman Filter)

# lorenz96_module.jl
Module script: Functions of functions to solve the nonlinear model in Lorenz (1996)

# lorenz96_NMC.jl
Script to create a binary file of Pstat.bin (provides a static background covariance error matrix based on the NMC method). 

# lorenz96.jl
EKF data assimilation system for the nonlinear model in Lorenz (1996)

## Basic equations of the Lorenz96 model
\\[a+b=C\\]
$$a+b+C$$

# lorenz96_hybrid.jl
Hybrid EKF data assimilation system for the nonlinear model in Lorenz (1996)

## Required modules
lorenz96_module.jl

## Required files
You need to run lorenz96_NMC.jl to create a static background covariance error matrix file (Pstat.bin). 

## Basic equations of the Lorenz96 model
\\[a+b=C\\]


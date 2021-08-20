# Documents
[**Formulation**](https://tomonori-93.github.io/julia/assimilation/KF/docs/formulation.html)

# Images
![X1-time_short](https://user-images.githubusercontent.com/11016219/130178988-d8dd749e-bce2-46ae-a74f-0b74d0c6297b.png)
![B-end_2d](https://user-images.githubusercontent.com/11016219/130179113-39ffeaaa-b72c-45d9-8fa3-95f21789ed1d.png)

# Modules
## EnKF_module.jl
Module script: Basic functions for EnKF (Ensemble Kalman Filter)

## lorenz96_module.jl
Module script: Functions of functions to solve the nonlinear model in Lorenz (1996)

# Scripts
**Note**: Tests of these scripts have been done on Jupyter-notebook with Julia-1.1.0. 

## lorenz96.jl
EKF (Extended Kalman Filter) data assimilation system for the nonlinear model in Lorenz (1996)

## lorenz96_NMC.jl
Script to create a binary file of Pstat.bin (provides a static background covariance error matrix based on the NMC method). 

## lorenz96_hybrid.jl
Hybrid EKF data assimilation system for the nonlinear model in Lorenz (1996)

- **Required modules**: [lorenz96_module.jl](lorenz96_module.jl)
- **Required files**: You need to run [lorenz96_NMC.jl](lorenz96_NMC.jl) to create a static background covariance error matrix file (Pstat.bin). 

## lorenz96_SEEK.jl
SEEK Filter data assimilation system for the nonlinear model in Lorenz (1996)

- **Required modules**: [lorenz96_module.jl](lorenz96_module.jl)

## lorenz96_LETKF.jl
LETKF data assimilation system for the nonlinear model in Lorenz (1996)

- **Required modules**: [lorenz96_module.jl](lorenz96_module.jl), [EnKF_module.jl](EnKF_module.jl)

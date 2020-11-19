# Making static background covariance error matrix by a NMC method in the Lorenz (1996) model
# Author: Satoki Tsujino (satoki_at_gfd-dennou.org)
# Date: 2020/11/19
# License: LGPL2.1
###########################
# calculation (main) part #
###########################

include("./KF/lorenz96_module.jl")

using .Lorenz96_functions
using LinearAlgebra
using Random
using Statistics

rng = MersenneTwister(1234)

nt = 15000
nx = 40
no = 40
n_ascyc = 5  # assimilation cycle interval for "t"
n_nmcyc = 20  # NMC cycle for "t"
rand_flag = true  # observation including random noise
nto = div(nt,n_ascyc) + 1
inf_fact = 1.05  # inflation factor
nspin = 10000  # spin-up steps
n_nmc = 365  # sample number for getting static background covariance error

# Allocating
y_o = reshape(zeros(no),no)
x_t = reshape(zeros(nx),nx)
x_spin = reshape(zeros(nx,nspin),nx,nspin)
x_nmc = reshape(zeros(nx,2),nx,2)
dd = reshape(zeros(nx,nx,n_nmc),nx,nx,n_nmc)
x_an = reshape(zeros(nx,nto),nx,nto)  # nx 行 nto 列行列へ変換
x_f = reshape(zeros(nx,nt),nx,nt)
Mop_linf = reshape(zeros(nx,nx),nx,nx)
Pf = reshape(zeros(nx,nx),nx,nx)
Pstat = reshape(zeros(nx,nx),nx,nx)
Pa = reshape(zeros(nx,nx),nx,nx)
Ro = reshape(zeros(no,no),no,no)
Hop = reshape(zeros(no,nx),no,nx)
Kg = reshape(zeros(nx,no),nx,no)
sigma_R = reshape(zeros(no,1),no,1)
I_mat = Matrix{Float64}(I,nx,nx)  # 単位行列の利用

# Setting parameters
dt = 0.01  # Time step for integration
F = 8.0  # default: 8
sigma_const_R = 1.0
Pf = (1.0 .* I_mat) + fill(20.0,nx,nx)
Ro = sigma_const_R * I_mat[1:no,1:no]
for i in 1:no
    Hop[i,i] = 1.0
end

xinit = fill(1.0,nx)
xinit[1] = xinit[1] + 0.1

# spin-up simulation and making initial states
x_spin[1:nx,1] = xinit
for i in 1:nspin-1
    x_spin[1:nx,i+1] = L96_RK4(x_spin[1:nx,i],dt,F,nx)
end

x_t[1:nx] = x_spin[1:nx,nspin]
for i in 1:nx
    x_f[i,1] = mean(x_spin[i,1:nspin])
end

# Main assimilation-prediction cycle

for i in 1:nt-1
    if mod(i-1,n_ascyc) == 0  # Entering the analysis processes
        #-- Analysis
        L_inv = inv(Hop * Pf * Hop' + Ro)
        Kg = (Pf * Hop') * L_inv
        Pa = Pf - (Kg * Hop) * Pf
        Pa = inf_fact .* Pa  # covariance inflation
        if rand_flag == true
            sigma_R = randn(rng, Float64) * map( x->sqrt(x), diag(Ro) )  # Ro の対角成分を抽出し, 平方根をとる.
            y_o = x_t[1:no] + sigma_R[1:no,1]
        else
            y_o = x_t[1:no]
        end
        d_innov = y_o - x_f[1:no,i]
        x_an[1:nx,div(i-1,n_ascyc)+1] = x_f[1:nx,i] + Kg * d_innov
        x_ftmp = x_f[1:nx,i] + Kg * d_innov
    else
        Pa = Pf
        x_ftmp = x_f[1:nx,i]
    end
    
    #-- Forecast
    Mop_linf = L96_RK4_tangent(x_ftmp,dt,F,nx)
    x_t = L96_RK4(x_t,dt,F,nx)
    x_f[1:nx,i+1] = L96_RK4(x_ftmp,dt,F,nx)
    Pf = (Mop_linf * Pa) * Mop_linf'
end

#-- NMC procedure
# nto: analysis cycle number, div(n_nmcyc,n_ascyc): NMC cycle on analysis cycle coordinate
# n_nmc: NMC cycle number
i_nmc = 1
for i in nto - div(n_nmcyc,n_ascyc) * n_nmc - 1:nto
    if mod(i - (nto - div(n_nmcyc,n_ascyc) * n_nmc - 1) + 1,div(n_nmcyc,n_ascyc)) == 0
        x_nmc[1:nx,1] = x_an[1:nx,i]  # One analysis step
        x_nmc[1:nx,2] = x_an[1:nx,i+1]  # The next analysis step
        for j in 1:n_nmcyc
            x_nmc[1:nx,1] = L96_RK4(x_nmc[1:nx,1],dt,F,nx)
        end
        for j in 1:n_nmcyc - n_ascyc
            x_nmc[1:nx,2] = L96_RK4(x_nmc[1:nx,2],dt,F,nx)
        end
        dx_nmc = x_nmc[1:nx,1] - x_nmc[1:nx,2]
        dd[1:nx,1:nx,i_nmc] = dx_nmc * dx_nmc'
        i_nmc = i_nmc + 1
    end
end

for i in 1:nx
    for j in 1:nx
        Pstat[i,j] = mean(dd[i,j,1:n_nmc])
    end
end

open("Pstat.bin","w") do io
    write(io,Pstat)
end

println("Output file Pstat.bin.")

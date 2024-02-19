#include("./KF/EnKF_module.jl")
#include("./KF/lorenz96_module.jl")

#using .EnKF_functions
#using .Lorenz96_functions
using LinearAlgebra
#using Random
#using Statistics

#rng = MersenneTwister(1234)

# set parameters
α_min = 0.0  # minimum alpha
α_max = 1.0  # maximum alpha
Re = 500.0  # Raynolds number
ε_min = 0.0  # minimum angle
ε_max = 1.0  # maximum angle
N = 400  # grid number
d2r = π / 180.0  # degree to radian

# initialize variables and matrix
ntot = 2*N-1  # matrix dimension
I_mat = Matrix{Complex{Float64}}(I,ntot,ntot)
X = reshape(zeros(ntot,1),ntot,1)
### --- From here, supposing matrix element as row, col as in mathematical notation
A = reshape(zeros(Complex{Float64},ntot,ntot),ntot,ntot)  # linear operator in differential equations
B = reshape(zeros(Complex{Float64},ntot,ntot),ntot,ntot)  # linear operator in differential equations
A_φφ = reshape(zeros(Complex{Float64},N-1,N-1),N-1,N-1)  # partial matrix (phi operator on phi EQ.) in A
A_φμ = reshape(zeros(Complex{Float64},N-1,N),N-1,N)      # partial matrix (mu operator on phi EQ.) in A
A_μφ = reshape(zeros(Complex{Float64},N,N-1),N,N-1)      # partial matrix (phi operator on mu EQ.) in A
A_μμ = reshape(zeros(Complex{Float64},N,N),N,N)          # partial matrix (mu operator on mu EQ.) in A
B_φφ = reshape(zeros(Complex{Float64},N-1,N-1),N-1,N-1)  # partial matrix (phi operator on phi EQ.) in B
B_φμ = reshape(zeros(Complex{Float64},N-1,N),N-1,N)      # partial matrix (mu operator on phi EQ.) in B == 0
B_μφ = reshape(zeros(Complex{Float64},N,N-1),N,N-1)      # partial matrix (phi operator on mu EQ.) in B == 0
B_μμ = reshape(zeros(Complex{Float64},N,N),N,N)          # partial matrix (mu operator on mu EQ.) in B
### --- From here, supposing matrix element as col, row as in julia notation (only using these arrays for matrix calc)
A_J = reshape(zeros(Complex{Float64},ntot,ntot),ntot,ntot)  # linear operator in differential equations
B_J = reshape(zeros(Complex{Float64},ntot,ntot),ntot,ntot)  # linear operator in differential equations

α = 0.45  # 0.5*(α_max + α_min)
ε = 0.0 * d2r  # 0.5*(ε_max + ε_min)
Δ = 1.0/sqrt(α * N)  # grid spacing
z_p   = [i*Δ for i in 1:N-1]  # non-dimensional vertical coordinate for φ
zh_p  = [(i-0.5)*Δ for i in 1:N]  # non-dimensional vertical coordinate for μ
U_bar = reshape(zeros(N,1),N,1)          # Mean flow in x-direction defined at i+1/2
V_bar = reshape(zeros(N-1,1),N-1,1)      # Mean flow in y-direction defined at i
Vh_bar = reshape(zeros(N,1),N,1)      # Mean flow in y-direction defined at i+1/2
dUdz  = reshape(zeros(N,1),N,1)          # Mean flow shear in x-direction defined at i+1/2
dVdz  = reshape(zeros(N-1,1),N-1,1)      # Mean flow shear in y-direction defined at i
d2Vdz2  = reshape(zeros(N-1,1),N-1,1)      # Mean flow shear in y-direction defined at i

z_p = reshape(z_p,N-1,1)
zh_p = reshape(zh_p,N,1)
U_bar[1:N,1] .= cos(ε) .- exp.(-zh_p[1:N,1]) .* cos.(zh_p[1:N,1] .+ ε)
V_bar[1:N-1,1] .= -sin(ε) .+ exp.(-z_p[1:N-1,1]) .* sin.(z_p[1:N-1,1] .+ ε)
Vh_bar[1:N,1] .= -sin(ε) .+ exp.(-zh_p[1:N,1]) .* sin.(zh_p[1:N,1] .+ ε)
dUdz[1:N,1] .= exp.(-zh_p[1:N,1]) .* cos.(zh_p[1:N,1] .+ ε) .+ exp.(-zh_p[1:N,1]) .* sin.(zh_p[1:N,1] .+ ε)
dVdz[1:N-1,1] .= -exp.(-z_p[1:N-1,1]) .* sin.(z_p[1:N-1,1] .+ ε) .+ exp.(-z_p[1:N-1,1]) .* cos.(z_p[1:N-1,1] .+ ε)
d2Vdz2[1:N-1,1] .= .- 2.0 .* exp.(-z_p[1:N-1,1]) .* cos.(z_p[1:N-1,1] .+ ε)

# setting matrix elements in A and B
Δ2 = Δ*Δ
Δ4 = Δ2*Δ2
αΔ = α*Δ
V_infl = -sin(ε) + exp(-(0.5*π-ε))

#      μ[1]        μ[2]       μ[3]        μ[4]
# ----- 0 --------- 1 --------- 2 --------- 3 -- i for μ = μ[i+1]
# ------|-----------|-----------|-----------|---   for μ = μ[i+1]
# 0---------- 1 --------- 2 --------- 3 -------- i for φ = φ[i]
# |-----------|-----------|-----------|---------   for φ = φ[i]
# For φ EQ.
for i in 3:N-3
#    println("$i")
   A_φφ[i,i+2] = complex(1.0,α*Re*V_bar[i,1]*Δ2/12.0)
   A_φφ[i,i-2] = A_φφ[i,i+2]
   A_φφ[i,i+1] = complex(-4.0-2.0*((αΔ)^2),-4.0*α*Re*V_bar[i,1]*Δ2/3.0)
   A_φφ[i,i-1] = A_φφ[i,i+1]
   A_φφ[i,i]   = complex(6.0+4.0*((αΔ)^2)+(α^4-d2Vdz2[i,1])*Δ4,α*Re*V_bar[i,1]*(2.5+(αΔ)^2)*Δ2)
   A_φμ[i,i]   = -2.0*(Δ^3)  # A^{φ,μ}_{i,i-1/2}
   A_φμ[i,i+1] = +2.0*(Δ^3)  # A^{φ,μ}_{i,i+1/2}
   B_φφ[i,i+2] = complex(0.0,α*Re*Δ2/12.0)
   B_φφ[i,i-2] = B_φφ[i,i+2]
   B_φφ[i,i+1] = -complex(0.0,4.0*α*Re*Δ2/3.0)
   B_φφ[i,i-1] = B_φφ[i,i+1]
   B_φφ[i,i]   = complex(0.0,α*Re*(2.5+(αΔ)^2)*Δ2)
end

# For μ EQ.
for i in 3:N-2
   A_μφ[i,i] = -complex(2.0*Δ,9.0*α*Re*dUdz[i,1]*Δ2/16.0)  # A^{μ,φ}_{i+1/2,i+1}
   A_μφ[i,i-2] = complex(0.0,α*Re*dUdz[i,1]*Δ2/16.0)  # A^{μ,φ}_{i+1/2,i-1}
   A_μφ[i,i+1] = A_μφ[i,i-2]  # A^{μ,φ}_{i+1/2,i+2}
   A_μφ[i,i-1] = complex(2.0*Δ,-9.0*α*Re*dUdz[i,1]*Δ2/16.0)  # A^{μ,φ}_{i+1/2,i}
   A_μμ[i,i+1] = 1.0  # A^{μ,μ}_{i+1/2,i-1/2}
   A_μμ[i,i-1] = 1.0  # A^{μ,μ}_{i+1/2,i+3/2}
   A_μμ[i,i]   = -complex(2.0+(αΔ)^2,α*Re*Vh_bar[i,1]*Δ2)
   B_μμ[i,i]   = -complex(0.0,α*Re*Δ2)
end

# Adjacent of boundaries
# For φ EQ. (at z=z_1)
A_φφ[1,3] = complex(1.0,α*Re*V_bar[1,1]*Δ2/12.0)
A_φφ[1,2] = complex(-4.0-2.0*((αΔ)^2),-4.0*α*Re*V_bar[1,1]*Δ2/3.0)
A_φφ[1,1] = complex(7.0+4.0*((αΔ)^2)+(α^4-d2Vdz2[1,1])*Δ4,α*Re*V_bar[1,1]*(31.0/12.0+(αΔ)^2)*Δ2)
A_φμ[1,1] = -2.0*(Δ^3)
A_φμ[1,2] = +2.0*(Δ^3)
B_φφ[1,3] = complex(0.0,α*Re*Δ2/12.0)
B_φφ[1,2] = -complex(0.0,4.0*α*Re*Δ2/3.0)
B_φφ[1,1] = complex(0.0,α*Re*(31.0/12.0+(αΔ)^2)*Δ2)
# For μ EQ. (at z=z_{1/2})
A_μφ[1,1] = -complex(2.0*Δ,0.5*α*Re*dUdz[1,1]*Δ2)
A_μφ[1,2] = complex(0.0,α*Re*dUdz[1,1]*Δ2/16.0)
A_μμ[1,2] = 1.0
A_μμ[1,1] = -complex(3.0+(αΔ)^2,α*Re*Vh_bar[1,1]*Δ2)
B_μμ[1,1] = -complex(0.0,α*Re*Δ2)

# For φ EQ. (at z=z_2)
A_φφ[2,4] = complex(1.0,α*Re*V_bar[2,1]*Δ2/12.0)
A_φφ[2,3] = complex(-4.0-2.0*((αΔ)^2),-4.0*α*Re*V_bar[2,1]*Δ2/3.0)
A_φφ[2,1] = A_φφ[2,3]
A_φφ[2,2] = complex(6.0+4.0*((αΔ)^2)+(α^4-d2Vdz2[2,1])*Δ4,α*Re*V_bar[2,1]*(2.5+(αΔ)^2)*Δ2)
A_φμ[2,2] = -2.0*(Δ^3)
A_φμ[2,3] = +2.0*(Δ^3)
B_φφ[2,4] = complex(0.0,α*Re*Δ2/12.0)
B_φφ[2,3] = -complex(0.0,4.0*α*Re*Δ2/3.0)
B_φφ[2,1] = B_φφ[2,3]
B_φφ[2,2] = complex(0.0,α*Re*(2.5+(αΔ)^2)*Δ2)
# For μ EQ. (at z=z_{1+1/2})
A_μφ[2,2] = -complex(2.0*Δ,9.0*α*Re*dUdz[2,1]*Δ2/16.0)
A_μφ[2,3] = complex(0.0,α*Re*dUdz[2,1]*Δ2/16.0)
A_μφ[2,1] = complex(2.0*Δ,-9.0*α*Re*dUdz[2,1]*Δ2/16.0)
A_μμ[2,3] = 1.0
A_μμ[2,1] = 1.0
A_μμ[2,2] = -complex(2.0+(αΔ)^2,α*Re*Vh_bar[2,1]*Δ2)
B_μμ[2,2] = -complex(0.0,α*Re*Δ2)

# For φ EQ. (at z=z_{N-2})
A_φφ[N-2,N-4] = complex(1.0,α*Re*V_bar[N-2,1]*Δ2/12.0)
A_φφ[N-2,N-1] = complex(-4.0-2.0*((αΔ)^2),-4.0*α*Re*V_bar[N-2,1]*Δ2/3.0)
A_φφ[N-2,N-3] = A_φφ[N-2,N-1]
A_φφ[N-2,N-2] = complex(6.0+4.0*((αΔ)^2)+(α^4-d2Vdz2[N-2,1])*Δ4,α*Re*V_bar[N-2,1]*(2.5+(αΔ)^2)*Δ2)
A_φμ[N-2,N-2] = -2.0*(Δ^3)
A_φμ[N-2,N-1] = +2.0*(Δ^3)
B_φφ[N-2,N-4] = complex(0.0,α*Re*Δ2/12.0)
B_φφ[N-2,N-1] = -complex(0.0,4.0*α*Re*Δ2/3.0)
B_φφ[N-2,N-3] = B_φφ[N-2,N-1]
B_φφ[N-2,N-2] = complex(0.0,α*Re*(2.5+(αΔ)^2)*Δ2)
# For μ EQ. (at z=z_{N-1-1/2})
A_μφ[N-1,N-1] = -complex(2.0*Δ,9.0*α*Re*dUdz[N-1,1]*Δ2/16.0)
A_μφ[N-1,N-3] = complex(0.0,α*Re*dUdz[N-1,1]*Δ2/16.0)
A_μφ[N-1,N-2] = complex(2.0*Δ,-9.0*α*Re*dUdz[N-1,1]*Δ2/16.0)
A_μμ[N-1,N-2] = 1.0
A_μμ[N-1,N]   = 1.0
A_μμ[N-1,N-1] = -complex(2.0+(αΔ)^2,α*Re*Vh_bar[N-1,1]*Δ2)
B_μμ[N-1,N-1] = -complex(0.0,α*Re*Δ2)

# For φ EQ. (at z=z_{N-1})
A_φφ[N-1,N-3]   = complex(1.0,α*Re*V_bar[N-1,1]*Δ2/12.0)
A_φφ[N-1,N-2]   = complex(-4.0-2.0*((αΔ)^2),-4.0*α*Re*V_bar[N-1,1]*Δ2/3.0)
A_φφ[N-1,N-1]   = complex(5.0+4.0*((αΔ)^2)+(α^4-d2Vdz2[N-1,1])*Δ4,α*Re*V_bar[N-1,1]*(29.0/12.0+(αΔ)^2)*Δ2)
A_φμ[N-1,N-1]   = -2.0*(Δ^3)
A_φμ[N-1,N]     = +2.0*(Δ^3)
B_φφ[N-1,N-3]   = complex(0.0,α*Re*Δ2/12.0)
B_φφ[N-1,N-2]   = -complex(0.0,4.0*α*Re*Δ2/3.0)
B_φφ[N-1,N-1]   = complex(0.0,α*Re*(29.0/12.0+(αΔ)^2)*Δ2)
# For μ EQ. (at z=z_{N-1/2})
A_μφ[N,N-2]     = complex(0.0,α*Re*dUdz[N]*Δ2/16.0)
A_μφ[N,N-1]     = complex(2.0*Δ,-5.0*α*Re*dUdz[N]*Δ2/8.0)
A_μμ[N,N-1]     = 1.0
A_μμ[N,N]       = -complex(1.0+(αΔ)^2,α*Re*Vh_bar[N]*Δ2)
B_μμ[N,N]       = -complex(0.0,α*Re*Δ2)

# assigned partial matrices A^ and B^ to A and B
A[1:N-1,1:N-1] = A_φφ[1:N-1,1:N-1]
A[1:N-1,N:2*N-1] = A_φμ[1:N-1,1:N]
A[N:2*N-1,1:N-1] = A_μφ[1:N,1:N-1]
A[N:2*N-1,N:2*N-1] = A_μμ[1:N,1:N]
# non-diagonal parts in B are zero
B[1:N-1,1:N-1] = B_φφ[1:N-1,1:N-1]
B[N:2*N-1,N:2*N-1] = B_μμ[1:N,1:N]

#for i in 3:ntot-3
#   println(A[i-1:i+1,i], B[i-1:i+1,i])
#end

A_J = transpose(A)  # mathematical A -> Julia A_J
B_J = transpose(B)  # mathematical A -> Julia A_J

C_J = inv(B_J) * A_J
D_J = eigvals(C_J)
F_J = imag.(D_J)

E_J = inv(B_J) * B_J

println(D_J[findmax(F_J)[2]])
println(F_J)
println(V_infl)

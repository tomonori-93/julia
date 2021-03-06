# Data assimilation simulation with the SEEK filter in Lorenz (1996) model.
# Author: Satoki Tsujino (satoki_at_gfd-dennou.org)
# Date: 2020/12/13
# Modification: 
# License: LGPL2.1
###########################
# calculation (main) part #
###########################

include("./lorenz96_module.jl")

using .Lorenz96_functions
using LinearAlgebra
using Random
using Statistics

rng = MersenneTwister(1234)

nt = 20000
nx = 40
no = 40
mx = 20
n_ascyc = 5  # assimilation cycle interval for "t"
fact_1day = 0.2  # Time for 1 day
rand_flag = true  # observation including random noise
nto = div(nt,n_ascyc) + 1
inf_fact = 5.0  # inflation factor
eps = 0.00001  # small value for tangential linear matrix of the model operator
nspin = 10000  # spin-up steps

# Allocating
t = zeros(nt)
t_o = zeros(nto)
y_o = reshape(zeros(no,nto),no,nto)
x_t = reshape(zeros(nx,nt),nx,nt)
x_spin = reshape(zeros(nx,nspin),nx,nspin)
x_an = reshape(zeros(nx),nx,1)
x_f = reshape(zeros(nx,nt),nx,nt)
x_e = reshape(zeros(nx,nt),nx,nt)
x_inc = reshape(zeros(nx,nto),nx,nto)
Mop_linf = reshape(zeros(nx,nx),nx,nx)
Pf = reshape(zeros(nx,nx),nx,nx)
Ro = reshape(zeros(no,no),no,no)
Df = reshape(zeros(mx,mx),mx,mx)
Uf = reshape(zeros(nx,mx),nx,mx)
Ufd = reshape(zeros(nx,mx),nx,mx)
Ua = reshape(zeros(nx,mx),nx,mx)
Da = reshape(zeros(mx,mx),mx,mx)
LD = reshape(zeros(mx,mx),mx,mx)
VEVT = reshape(zeros(mx,mx),mx,mx)
E12 = reshape(zeros(mx,mx),mx,mx)
Vevec = reshape(zeros(mx,mx),mx,mx)
Hop = reshape(zeros(no,nx),no,nx)
Kg = reshape(zeros(nx,no),nx,no)
d_innov = reshape(zeros(no,1),no,1)
Pftime = reshape(zeros(nx,nx,2),nx,nx,2)
Egtime = reshape(zeros(nt,1),nt,1)
Egftime = reshape(zeros(nt,1),nt,1)
L2Ntime = reshape(zeros(nt,1),nt,1)
L2Natime = reshape(zeros(nto,1),nto,1)
sigma_R = reshape(zeros(no,1),no,1)
delta_x = reshape(zeros(nx,1),nx,1)
I_mat = Matrix{Float64}(I,nx,nx)

# Setting parameters
dt = 0.01  # Time step for integration
n_1day = Int(fact_1day / dt)
delx = 0.1  # small departure for the TL check
F = 8.0  # default: 8
sigma_const_R = 1.0
Pf = (1.0 .* I_mat) + fill(1.0,nx,nx)
for i in 1:nx
    Pf[i,i] = Pf[i,i] + i
end
Ro = sigma_const_R * I_mat[1:no,1:no]
for i in 1:no
    Hop[i,i] = 1.0
end

xinit = fill(1.0,nx)
xinit[1] = xinit[1] + 0.1

eps_inv = 1.0 / eps

# spin-up simulation and making initial states
x_spin[1:nx,1] = xinit
for i in 1:nspin-1
    x_spin[1:nx,i+1] = L96_RK4(x_spin[1:nx,i],dt,F,nx)
end

x_t[1:nx,1] = x_spin[1:nx,nspin]
for i in 1:nx
    x_f[i,1] = mean(x_spin[i,1:nspin])
end
x_e[1:nx,1] = x_f[1:nx,1]

# TL check (If you need to check the Mop, please activate the next line.)
#L96_TL_check(x_t[1:nx,1],delx.*fill(1.0,nx,1),dt,1000,F,nx,"RK4")

# Main assimilation-prediction cycle

Pftime[1:nx,1:nx,1] = Pf

# Decomposition of Pf
lambda = eigen(Symmetric(Pf),1:nx)
for i in 1:mx
    Da[i,i] = lambda.values[nx - i + 1]  # set eigenvalues at the maximum to the mx th in the diagonal components
    Uf[1:nx,i] = lambda.vectors[1:nx,nx - i + 1]  # set Uf based on eigenvectors at the maximum to the mx th
end

for i in 1:nt-1
    if mod(i-1,n_ascyc) == 0  # Entering the analysis processes
        #println("Enter hear")
        #-- Analysis
        L_inv = inv(Hop * Uf * Da * Uf' * Hop' + Ro)
        Df = Da - Da * Uf' * Hop' * L_inv * Hop * Uf * Da
        Kg = Uf * Df * Uf' * Hop' * inv(Ro)
        
        LLT = cholesky(Hermitian(Df))
        LD = LLT.L
        VEVT = LD' * Uf' * Uf * LD
        lambda = eigen(Symmetric(VEVT),1:mx)
        Da = fill(0.0,mx,mx)
        E12 = fill(0.0,mx,mx)
        for j in 1:mx
            Da[j,j] = lambda.values[mx - j + 1]
            E12[j,j] = 1.0 / sqrt(Da[j,j])
            Vevec[1:mx,j] = lambda.vectors[1:mx,mx - j + 1]
        end
        Ua = Uf * LD * Vevec * E12
        Da = inf_fact .* Da  # covariance inflation
        
        if rand_flag == true
            sigma_R = randn(rng, Float64) * map( x->sqrt(x), diag(Ro) )
            y_o[1:no,div(i-1,n_ascyc)+1] = x_t[1:no,i] + sigma_R[1:no,1]
        else
            y_o[1:no,div(i-1,n_ascyc)+1] = x_t[1:no,i]
        end
        d_innov = y_o[1:no,div(i-1,n_ascyc)+1] - x_f[1:no,i]
        x_an = x_f[1:nx,i] + Kg * d_innov
        x_inc[1:nx,div(i-1,n_ascyc)+1] = x_an - x_f[1:nx,i]
        t_o[div(i-1,n_ascyc)+1] = dt*(i-1)
        
        L2Natime[div(i-1,n_ascyc)+1,1] = sqrt((x_an - x_t[1:nx,i])' * (x_an - x_t[1:nx,i]) / nx)
    else
        Ua = Uf
        x_an = x_f[1:nx,i]
    end
    
    #-- Forecast
    Mop_linf = L96_RK4_tangent(x_an,dt,F,nx)
    t[i] = dt*(i-1)
    
    x_t[1:nx,i+1] = L96_RK4(x_t[1:nx,i],dt,F,nx)
    x_e[1:nx,i+1] = L96_RK4(x_e[1:nx,i],dt,F,nx)
    x_f[1:nx,i+1] = L96_RK4(x_an,dt,F,nx)
    Uf = Mop_linf * Ua
    #for j in 1:mx
    #    Uf[1:nx,j] = eps_inv .* (L96_RK4(x_an + eps .* Ua[1:nx,j],dt,F,nx) - L96_RK4(x_an,dt,F,nx))
    #end

    Egftime[i+1] = tr(Uf * Df * Uf') / nx
    Egtime[i+1] = tr(Ua * Da * Ua') / nx
    delta_x = x_f[1:nx,i+1] - x_t[1:nx,i+1]
    L2Ntime[i+1] = sqrt((delta_x' * delta_x) / nx)
    
end
t[nt] = t[nt-1] + dt

#Pftime[1:nx,1:nx,2] = Pf

#-- Drawing
##########
#  Plot  #
##########
using PyPlot

fig = figure("pyplot_majorminor",figsize=(7,5))

draw_num = 3
if draw_num == 1  # X-time
    p = plot(5.0.*t[1:nt],x_t[1,1:nt],color="black",label="Perfect")
    p = plot(5.0.*t[1:nt],x_e[1,1:nt],color="blue",label="False")
    p = plot(5.0.*t[1:nt],x_f[1,1:nt],color="red",label="Assim")
    p = plot(5.0.*t[1:nt],Egtime[1:nt],color="green",label="Tr(Pa)")
    p = plot(5.0.*t[1:nt],L2Ntime[1:nt],color="orange",label="L2x-xt")
    p = plot(5.0.*t_o[1:div(nt-1,n_ascyc)+1],y_o[1,1:div(nt-1,n_ascyc)+1],linestyle="",marker="o",color="black",label="Obs")
    
elseif draw_num == 2  # X-Z phase
    p = plot(x_t[1,1:nt],x_t[3,1:nt],color="black",label="Perfect")
    p = plot(x_e[1,1:nt],x_e[3,1:nt],color="blue",label="False")
    p = plot(x_f[1,1:nt],x_f[3,1:nt],color="red",label="Assim")
    p = plot(y_o[1,1:div(nt-1,n_ascyc)+1],y_o[3,1:div(nt-1,n_ascyc)+1],linestyle="",marker="o",color="black",label="Obs")
    
elseif draw_num == 3  # error-time
    p = plot(5.0.*t[1:nt],Egftime[1:nt],color="black",label="Tr(Pf)")
#    p = plot(5.0.*t[1:nt],Egtime[1:nt],color="blue",label="Tr(Pa)")
    p = plot(5.0.*t[1:nt],L2Ntime[1:nt],color="red",label="RMSE(x-xt)")
#    p = plot(5.0.*t_o[1:div(nt-1,n_ascyc)+1],L2Natime[1:div(nt-1,n_ascyc)+1,1],color="blue",label="RMSE(xa-xt)")
    p = plot(5.0.*t_o[1:div(nt-1,n_ascyc)+1],fill(0.0,div(nt-1,n_ascyc)+1),linestyle="",marker="o",color="black",label="Obs")
#    p = plot(t[1:nt],lam_max[1:nt],color="green",label="Lambda_max")
        
end

legend()
ax = gca()

if draw_num == 1  # X-time
    xlabel("Time (days)")
    ylabel("X^T")
    
elseif draw_num == 2  # X-Z phase
    xlabel("X")
    ylabel("Z")
    
elseif draw_num == 3  # error-time
    xlabel("Time (days)")
    ylabel("ΔX")

end
    
grid("on")
PyPlot.title("Lorenz96")

###########################
#  Set the tick interval  #
###########################
###Mx = matplotlib.ticker.MultipleLocator(0.1) # Define interval of major ticks
##Mx = matplotlib.ticker.MultipleLocator(1.0) # Define interval of major ticks
#Mx = matplotlib.ticker.MultipleLocator(4.0) # Define interval of major ticks
#f = matplotlib.ticker.FormatStrFormatter("%1.2f") # Define format of tick labels
##ax.xaxis.set_major_locator(Mx) # Set interval of major ticks
#ax.xaxis.set_major_formatter(f) # Set format of tick labels

###mx = matplotlib.ticker.MultipleLocator(0.02) # Define interval of minor ticks
##mx = matplotlib.ticker.MultipleLocator(0.1) # Define interval of minor ticks
#mx = matplotlib.ticker.MultipleLocator(1.0) # Define interval of minor ticks
##ax.xaxis.set_minor_locator(mx) # Set interval of minor ticks

###My = matplotlib.ticker.MultipleLocator(2.0) # Define interval of major ticks
##My = matplotlib.ticker.MultipleLocator(2.0) # Define interval of major ticks
#My = matplotlib.ticker.MultipleLocator(4.0) # Define interval of major ticks
##ax.yaxis.set_major_locator(My) # Set interval of major ticks

###my = matplotlib.ticker.MultipleLocator(0.5) # Define interval of minor ticks
##my = matplotlib.ticker.MultipleLocator(0.5) # Define interval of minor ticks
#my = matplotlib.ticker.MultipleLocator(1.0) # Define interval of minor ticks
##ax.yaxis.set_minor_locator(my) # Set interval of minor ticks

#########################
#  Set tick dimensions  #
#########################
#ax.xaxis.set_tick_params(which="major",length=5,width=2,labelsize=10)
#ax.xaxis.set_tick_params(which="minor",length=5,width=2)

fig.canvas.draw() # Update the figure
gcf() # Needed for IJulia to plot inline
savefig("Error-time.pdf")
##savefig("X-time.pdf")
#savefig("X-Z_lorenz.pdf")


xax = reshape(zeros(nx,nt),nx,nt)
tax = reshape(zeros(nx,nt),nx,nt)
xoax = reshape(zeros(nx,nto),nx,nto)
toax = reshape(zeros(nx,nto),nx,nto)

for i in 1:nx
    xax[i,1:nt] .= i
    xoax[i,1:nto] .= i
end
for i in 1:nt
    tax[1:nx,i] .= t[i] * 5.0
end
for i in 1:nto
    toax[1:nx,i] .= t_o[i] * 5.0
end


#-- Drawing partII
##########
#  Plot  #
##########
rc("font", family="IPAPGothic")
fig = figure("pyplot_majorminor",figsize=(7,5))

draw_num2 = 3
if draw_num2 == 1
    #cp = contourf(xax[1:nx,1:nx], xax'[1:nx,1:nx], Pftime[1:nx,1:nx,1], levels=[-15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0], origin="image", cmap=ColorMap("viridis"), extend="both")
    cp = contourf(xax[1:nx,1:nx], xax'[1:nx,1:nx], Pftime[1:nx,1:nx,2], origin="image", cmap=ColorMap("viridis"), extend="both")
elseif draw_num2 == 2
    cp = contourf(xax, tax, x_f, 10, levels=[-15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0], cmap=ColorMap("viridis"), extend="both")
elseif draw_num2 == 3
    println(evec_max[1:nx,nto-1])
    cp = contourf(xoax[1:nx,1:nto-1], toax[1:nx,1:nto-1], evec_max[1:nx,1:nto-1], 10, cmap=ColorMap("viridis"), extend="both")
end
#ax.label(cp, inline=1, fontsize=10)
#legend()
ax = gca()

if draw_num2 == 1 || draw_num2 == 2
    xlabel("変数 (X_k)")
    ylabel("変数 (X_k)")
elseif draw_num2 == 3
    xlabel("変数 (X_k)")
    ylabel("時間 (days)")
end
plt.colorbar(cp)
grid("on")

if draw_num2 == 1 || draw_num2 == 2
    PyPlot.title("Lorenz96 (Pf)")
elseif draw_num2 == 3
    PyPlot.title("Lorenz96 (Eigenvec_λmax)")
end

#########################
#  Set tick dimensions  #
#########################
#ax.xaxis.set_tick_params(which="major",length=5,width=2,labelsize=10)
#ax.xaxis.set_tick_params(which="minor",length=5,width=2)

fig.canvas.draw() # Update the figure
gcf() # Needed for IJulia to plot inline
#savefig("B-init_2d.pdf")
savefig("B-lammax_2d.pdf")

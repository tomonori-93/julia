# KF assimilation simulation in Lorenz (1963) model.
# Author: Satoki Tsujino (satoki_at_gfd-dennou.org)
# Date: 2020/11/01
# Modification: 2020/11/20
# License: LGPL2.1
#############################
# functions for calculation #
#############################

include("./KF/lorenz63_module.jl")

using .Lorenz63_functions
using Random
using LinearAlgebra

rng = MersenneTwister(1234)

nt = 1000000
nx = 3
no = 3
n_ascyc = 250000  # assimilation cycle interval for "t"
rand_flag = false  # observation including random noise
nto = div(nt,n_ascyc) + 1
inf_fact = 1.0  # inflation factor
eps = 0.0005  # small value for tangential linear matrix of the model operator

# Allocating
t = zeros(nt)
t_o = zeros(nto)
y_o = reshape(zeros(no,nto),no,nto)
x_t = reshape(zeros(nx,nt),nx,nt)
x_an = reshape(zeros(nx),nx,1)  # nx 行 1 列行列へ変換
x_f = reshape(zeros(nx,nt),nx,nt)
x_e = reshape(zeros(nx,nt),nx,nt)
Mop_linf = reshape(zeros(nx,nx),nx,nx)
Pf = reshape(zeros(nx,nx),nx,nx)
Pa = reshape(zeros(nx,nx),nx,nx)
Ro = reshape(zeros(no,no),no,no)
Hop = reshape(zeros(no,nx),no,nx)
Kg = reshape(zeros(nx,no),nx,no)
Pftime = reshape(zeros(nx,nt),nx,nt)
sigma_R = reshape(zeros(no,1),no,1)
evec_max = reshape(zeros(nx,nt),nx,nt)
lam_max = reshape(zeros(nt),nt)
Egtime = reshape(zeros(nt,1),nt,1)
L2Ntime = reshape(zeros(nt,1),nt,1)
delta_x = reshape(zeros(nx,1),nx,1)
I_mat = Matrix{Float64}(I,nx,nx)  # 単位行列の利用

# Setting parameters
dt = 0.00001
delx = 0.1  # small departure for the TL check
k = 10.0
r = 23.0
b = 8.0/3.0
Pf = [20.0 20.0 20.0
      20.0 20.0 20.0
      20.0 20.0 20.0]
Ro = [1.0 0.0 0.0
      0.0 1.0 0.0
      0.0 0.0 1.0]  # 行ベクトル [a b], 列ベクトル [a, b]

xinit = [10.0, 20.0, 30.0]

x_t[1:nx,1] = xinit[1:nx]
x_f[1:nx,1] = x_t[1:nx,1] + 0.1 * xinit[1:nx]
x_e[1:nx,1] = x_f[1:nx,1]
Hop = I_mat
lam_max[1] = 1.0
evec_max = fill(1.0,nx,nt)
eps_inv = 1.0 / eps

# TL check (If you need to check the Mop, please activate the next line.)
# L63_TL_check(x_t[1:nx,1],delx.*fill(1.0,nx,1),dt,100,k,r,b)

# Main assimilation-prediction cycle
for i in 1:nt-1
    if mod(i-1,n_ascyc) == 0  # Entering the analysis processes
        println("Enter hear")
        #-- Analysis
        L_inv = inv(Hop * Pf * Hop' + Ro)
        Kg = (Pf * Hop') * L_inv
        Pa = Pf - (Kg * Hop) * Pf
        if rand_flag == true
            sigma_R = randn(rng, Float64) * map( x->sqrt(x), diag(Ro) )  # Ro の対角成分を抽出し, 平方根をとる.
            y_o[1:no,div(i-1,n_ascyc)+1] = x_t[1:no,i] + sigma_R[1:no,1]
        else
            y_o[1:no,div(i-1,n_ascyc)+1] = x_t[1:no,i]
        end
        d_innov = y_o[1:no,div(i-1,n_ascyc)+1] - x_f[1:no,i]        
        x_an = x_f[1:nx,i] + Kg * d_innov
        t_o[div(i-1,n_ascyc)+1] = dt*(i-1)
    else
        Pa = Pf
        x_an = x_f[1:nx,i]
    end

    #-- Forecast
    t[i] = dt*(i-1)
    
    Mop_linf = L63_EU1_tangent(x_an,dt,k,r,b)
    #for j in 1:nx
    #    Mop_linf[1:nx,j] = eps_inv .* (L63_RK4(x_an_lin + eps .* I_mat[1:nx,j],dt,k,r,b) - L63_RK4(x_an_lin,dt,k,r,b))
    #end

    x_t[1:nx,i+1] = L63_EU1(x_t[1:nx,i],dt,k,r,b)
    x_e[1:nx,i+1] = L63_EU1(x_e[1:nx,i],dt,k,r,b)
    x_f[1:nx,i+1] = L63_EU1(x_an,dt,k,r,b)
    Pf = inf_fact * Mop_linf * Pa * Mop_linf'
    Pftime[1:nx,i+1] = diag(Pf)
    Egtime[i+1] = tr(Pf) / nx
    delta_x = x_f[1:nx,i+1] - x_t[1:nx,i+1]
    L2Ntime[i+1] = sqrt(delta_x' * delta_x ./ nx)
    
    if mod(i-1,n_ascyc) == 0
        lambda = eigen(Symmetric(Mop_linf' * Mop_linf),nx:nx)
        lam_max[i+1] = lambda.values[1]  # 最大固有値
        evec_max[1:nx,i+1] = lambda.vectors[1:nx,1]  # 最大固有値の固有ベクトル
    end

end

t[nt] = t[nt-1] + dt


#-- Drawing
##########
#  Plot  #
##########
using PyPlot

fig = figure("pyplot_majorminor",figsize=(7,5))
draw_num = 3
if draw_num == 1
    p = plot(t[1:nt],x_t[1,1:nt],color="black",label="Perfect")
    p = plot(t[1:nt],x_e[1,1:nt],color="blue",label="False")
    p = plot(t[1:nt],x_f[1,1:nt],color="red",label="Assim")
    p = plot(t_o[1:div(nt-1,n_ascyc)+1],y_o[1,1:div(nt-1,n_ascyc)+1],linestyle="",marker="o",color="black",label="Obs")
elseif draw_num == 2
    p = plot(x_t[1,1:nt],x_t[3,1:nt],color="black",label="Perfect")
    p = plot(x_e[1,1:nt],x_e[3,1:nt],color="blue",label="False")
    p = plot(x_f[1,1:nt],x_f[3,1:nt],color="red",label="Assim")
    p = plot(y_o[1,1:div(nt-1,n_ascyc)+1],y_o[3,1:div(nt-1,n_ascyc)+1],linestyle="",marker="o",color="black",label="Obs")
elseif draw_num == 3
    p = plot(t[1:nt],Egtime[1:nt]*0.01,color="black",label="Tr(Pf)")
    p = plot(t[1:nt],L2Ntime[1:nt],color="red",label="RMSE(x-xt)")
    p = plot(t_o[1:div(nt-1,n_ascyc)+1],fill(0.0,div(nt-1,n_ascyc)+1),linestyle="",marker="o",color="black",label="Obs")
    #p = plot(t[1:nt],lam_max[1:nt],color="green",label="Lambda_max")
end
legend()
ax = gca()

if draw_num == 1
    xlabel("Time")
    ylabel("X^T")
elseif draw_num == 2
    xlabel("Time")
    ylabel("ΔX")
elseif draw_num == 3
    xlabel("X")
    ylabel("Z")
end

grid("on")
PyPlot.title("Lorenz63")

###########################
#  Set the tick interval  #
###########################
##Mx = matplotlib.ticker.MultipleLocator(1.0) # Define interval of major ticks
#Mx = matplotlib.ticker.MultipleLocator(4.0) # Define interval of major ticks
#f = matplotlib.ticker.FormatStrFormatter("%1.2f") # Define format of tick labels
##ax.xaxis.set_major_locator(Mx) # Set interval of major ticks
#ax.xaxis.set_major_formatter(f) # Set format of tick labels

##mx = matplotlib.ticker.MultipleLocator(0.1) # Define interval of minor ticks
#mx = matplotlib.ticker.MultipleLocator(1.0) # Define interval of minor ticks
##ax.xaxis.set_minor_locator(mx) # Set interval of minor ticks

##My = matplotlib.ticker.MultipleLocator(2.0) # Define interval of major ticks
#My = matplotlib.ticker.MultipleLocator(4.0) # Define interval of major ticks
##ax.yaxis.set_major_locator(My) # Set interval of major ticks

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

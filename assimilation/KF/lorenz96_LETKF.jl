# Data assimilation simulation with the LETKF in Lorenz (1996) model.
# Author: Satoki Tsujino (satoki_at_gfd-dennou.org)
# Date: 2020/12/20
# Modification: 
# License: LGPL2.1
###########################
# calculation (main) part #
###########################

include("./EnKF_module.jl")
include("./lorenz96_module.jl")

using .EnKF_functions
using .Lorenz96_functions
using LinearAlgebra
using Random
using Statistics

rng = MersenneTwister(1234)

rng = MersenneTwister(1234)

nt = 20000
nx = 40
no = 40
nex = 78  # ensemble member
n_ascyc = 5  # assimilation cycle interval for "t"
fact_1day = 0.2  # Time for 1 day
rand_flag = true  # observation including random noise
nto = div(nt,n_ascyc) + 1
inf_fact = 1.5  # inflation factor
nspin = 10000  # spin-up steps

# Allocating
t = zeros(nt)
t_o = zeros(nto)
y_o = reshape(zeros(no,nto),no,nto)
x_t = reshape(zeros(nx,nt),nx,nt)
x_spin = reshape(zeros(nx,nspin),nx,nspin)
x_f = reshape(zeros(nx,nt),nx,nt)
Xen = reshape(zeros(nx,nex),nx,nex)
X_an = reshape(zeros(nx,nex),nx,nex)
dXen = reshape(zeros(nx,nex),nx,nex)
Yo = reshape(zeros(no,nex),no,nex)
x_inc = reshape(zeros(nx,nto),nx,nto)
draw_Xen = reshape(zeros(nex,nt),nex,nt)
Pf = reshape(zeros(nx,nx),nx,nx)
Ro = reshape(zeros(no,no),no,no)
Rinv = reshape(zeros(no,no),no,no)
Da = reshape(zeros(nex,nex),nex,nex)
UDUT = reshape(zeros(nex,nex),nex,nex)
D12 = reshape(zeros(nex,nex),nex,nex)
Uevec = reshape(zeros(nex,nex),nex,nex)
Mop = reshape(zeros(nx,nx),nx,nx)
Hop = reshape(zeros(no,nx),no,nx)
Kg = reshape(zeros(nx,no),nx,no)
d_innov = reshape(zeros(no,nex),no,nex)
Xfmean = reshape(zeros(nx,nex),nx,nex)
HXmean = reshape(zeros(no,nex),no,nex)
Pftime = reshape(zeros(nx,nx,2),nx,nx,2)
Egtime = reshape(zeros(nt,1),nt,1)
Egftime = reshape(zeros(nt,1),nt,1)
L2Ntime = reshape(zeros(nt,1),nt,1)
L2Natime = reshape(zeros(nto,1),nto,1)
sigma_R = reshape(zeros(no,1),no,1)
delta_x = reshape(zeros(nx,1),nx,1)
I_mat = Matrix{Float64}(I,nx,nx)
I_mat_en = Matrix{Float64}(I,nex,nex)

# Setting parameters
dt = 0.01  # Time step for integration
n_1day = Int(fact_1day / dt)
delx = 0.1  # small departure for the TL check
F = 8.0  # default: 8
sigma_const_R = 1.0
#Pf = (1.0 .* I_mat) + fill(1.0,nx,nx)
#for i in 1:nx
#    Pf[i,i] = Pf[i,i] + i
#end
Ro = sigma_const_R * I_mat[1:no,1:no]
Rinv = I_mat[1:no,1:no] ./ sigma_const_R
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

x_t[1:nx,1] = x_spin[1:nx,nspin]
for i in 1:nx
    x_f[i,1] = mean(x_spin[i,1:nspin])
end

# make ensemble members
Mop = L96_RK4_tangent(x_f[1:nx,1],dt,F,nx)
lambda = eigen(Symmetric(Mop*Mop'),1:nx)
for i in 1:div(nex,2)
    Xen[1:nx,i] = x_f[1:nx,1] + lambda.vectors[1:nx,nx - i + 1]
    Xen[1:nx,nex - i + 1] = x_f[1:nx,1] - lambda.vectors[1:nx,nx - i + 1]
end

draw_Xen[1:nex,1] = Xen[1,1:nex]
# Main assimilation-prediction cycle

Pftime[1:nx,1:nx,1] = Pf

for i in 1:nt-1
    if mod(i-1,n_ascyc) == 0  # Entering the analysis processes
        #println("Enter hear")
        #-- Analysis
        #-- Calculation of ensemble mean Xf
        Xfmean = En_mean(Xen,nx,nex)
        HXmean = En_mean(Hop*Xen,no,nex)
        dXen = Xen - Xfmean
        dY = Hop * dXen
        UDUT = (nex - 1.0) .* I_mat_en + dY' * Rinv * dY
        lambda = eigen(Symmetric(UDUT),1:nex)
        for j in 1:nex
            Da[j,j] = lambda.values[nex - j + 1]
            D12[j,j] = 1.0 / sqrt(Da[j,j])
            Uevec[1:nex,j] = lambda.vectors[1:nex,nex - j + 1]
        end
        Tm = sqrt(nex - 1.0) .* Uevec * D12 * Uevec'
        Kg = dXen * Uevec * (D12 * D12) * Uevec' * dY' * Rinv
        
        if rand_flag == true
            sigma_R = randn(rng, Float64) * map( x->sqrt(x), diag(Ro) )  # Ro の対角成分を抽出し, 平方根をとる.
            y_o[1:no,div(i-1,n_ascyc)+1] = x_t[1:no,i] + sigma_R[1:no,1]
        else
            y_o[1:no,div(i-1,n_ascyc)+1] = x_t[1:no,i]
        end
        for j in 1:no
            Yo[j,1:nex] = fill(y_o[j,div(i-1,n_ascyc)+1],1,nex)
        end
        d_innov = Yo - HXmean
        X_an = Xfmean + inf_fact .* ( Kg * d_innov + dXen * Tm )
        #x_inc[1:nx,div(i-1,n_ascyc)+1] = x_an - x_f[1:nx,i]
        t_o[div(i-1,n_ascyc)+1] = dt*(i-1)
        
        #L2Natime[div(i-1,n_ascyc)+1,1] = sqrt((x_an - x_t[1:nx,i])' * (x_an - x_t[1:nx,i]) / nx)
    else
        X_an = Xen
    end
    
    #-- Forecast
    t[i] = dt*(i-1)
    
    x_t[1:nx,i+1] = L96_RK4(x_t[1:nx,i],dt,F,nx)
    for j in 1:nex
        Xen[1:nx,j] = L96_RK4(X_an[1:nx,j],dt,F,nx)
    end
    x_f[1:nx,i+1] = En_mean(Xen,nx,nex)[1:nx,1]
    draw_Xen[1:nex,i+1] = Xen[1,1:nex]

    #Egftime[i+1] = tr(Uf * Df * Uf') / nx
    #Egtime[i+1] = tr(Ua * Da * Ua') / nx
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

draw_num = 1
if draw_num == 1  # X-time
    p = plot(5.0.*t[1:200],draw_Xen[1,1:200],color="red",label="Ensemble")
    for i in 2:nex
        p = plot(5.0.*t[1:200],draw_Xen[i,1:200],color="red")
    end
    p = plot(5.0.*t[1:200],x_t[1,1:200],color="black",label="Perfect")
    p = plot(5.0.*t_o[1:div(200-1,n_ascyc)+1],y_o[1,1:div(200-1,n_ascyc)+1],linestyle="",marker="o",color="black",label="Obs")
    #p = plot(5.0.*t[1:nt],x_f[1,1:nt],color="red",label="Assim")
    #p = plot(5.0.*t[1:nt],Egtime[1:nt],color="green",label="Tr(Pa)")
    #p = plot(5.0.*t[1:nt],L2Ntime[1:nt],color="orange",label="L2x-xt")
    #p = plot(5.0.*t_o[1:div(nt-1,n_ascyc)+1],y_o[1,1:div(nt-1,n_ascyc)+1],linestyle="",marker="o",color="black",label="Obs")
    
elseif draw_num == 2  # X-Z phase
    p = plot(x_t[1,1:nt],x_t[3,1:nt],color="black",label="Perfect")
    p = plot(x_e[1,1:nt],x_e[3,1:nt],color="blue",label="False")
    p = plot(x_f[1,1:nt],x_f[3,1:nt],color="red",label="Assim")
    p = plot(y_o[1,1:div(nt-1,n_ascyc)+1],y_o[3,1:div(nt-1,n_ascyc)+1],linestyle="",marker="o",color="black",label="Obs")
    
elseif draw_num == 3  # error-time
#    p = plot(5.0.*t[1:nt],Egftime[1:nt],color="black",label="Tr(Pf)")
    p = plot(5.0.*t[1:nt],Egtime[1:nt],color="blue",label="Tr(Pa)")
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
PyPlot.title("Lorenz96-LETKF")

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

#-- Drawing
##########
#  Plot  #
##########
rc("font", family="IPAPGothic")
fig = figure("pyplot_majorminor",figsize=(7,5))

draw_num2 = 2
if draw_num2 == 1
    #cp = contourf(xax[1:nx,1:nx], xax'[1:nx,1:nx], Pftime[1:nx,1:nx,1], levels=[-15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0], origin="image", cmap=ColorMap("viridis"), extend="both")
    #cp = contourf(xax[1:nx,1:nx], xax'[1:nx,1:nx], Pftime[1:nx,1:nx,2], levels=[-2.0, -1.0, -0.5, -0.25, 0.25, 0.5, 1.0, 1.5, 2.0], origin="image", cmap=ColorMap("viridis"), extend="both")
    cp = contourf(xax[1:nx,1:no], xax'[1:nx,1:no], Pftime[1:nx,1:nx,2]*Hop'*inv(Hop*Pftime[1:nx,1:nx,2]*Hop'+Ro), origin="image", cmap=ColorMap("viridis"), extend="both")
    #cp = contourf(xax[1:nx,1:nx], xax'[1:nx,1:nx], Pstat[1:nx,1:nx], origin="image", cmap=ColorMap("viridis"), extend="both")
elseif draw_num2 == 2
    cp = contourf(xax[1:nx,1:nt], tax[1:nx,1:nt], x_f[1:nx,1:nt], 10, levels=[-15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0], cmap=ColorMap("viridis"), extend="both")
elseif draw_num2 == 3
    cp = contourf(xoax[1:nx,1:nto-1], toax[1:nx,1:nto-1], evec_max[1:nx,1:nto-1], 10, cmap=ColorMap("viridis"), extend="both")
elseif draw_num2 == 4
    cp = contourf(xoax[1:nx,1:nto], toax[1:nx,1:nto], x_inc[1:nx,1:nto], cmap=ColorMap("viridis"), extend="both")
elseif draw_num2 == 5
    cp = contourf(xoax[1:no,1:40], toax[1:no,1:40], y_innov[1:no,1:40], levels=[-2.0, -1.0, -0.5, -0.25, 0.25, 0.5, 1.0, 1.5, 2.0], cmap=ColorMap("viridis"), extend="both")
end
#ax.label(cp, inline=1, fontsize=10)
#legend()
ax = gca()

if draw_num2 == 1
    xlabel("変数 (X_k)")
    ylabel("変数 (X_k)")
elseif draw_num2 == 5 || draw_num2 == 4 || draw_num2 == 3 || draw_num2 == 2
    xlabel("変数 (X_k)")
    ylabel("時間 (days)")
end
plt.colorbar(cp)
grid("on")

if draw_num2 == 1
    PyPlot.title("Lorenz96 (Pf)")
elseif draw_num2 == 2
    PyPlot.title("Lorenz96 (Xf)")
elseif draw_num2 == 3
    PyPlot.title("Lorenz96 (Eigenvec_λmax)")
elseif draw_num2 == 4
    PyPlot.title("Lorenz96 (δx)")
elseif draw_num2 == 5
    PyPlot.title("Lorenz96 (D_innov)")
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

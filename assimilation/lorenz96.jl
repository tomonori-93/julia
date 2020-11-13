# Lorenz (1996) モデルで同化・解析処理を行うシミュレーション
# Author: Satoki Tsujino (satoki_at_gfd-dennou.org)
# Date: 2020/11/01
# License: LGPL2.1
#############################
# functions for calculation #
#############################

using LinearAlgebra

function L96_tendency(x,F,N)  # calculation of tendency in Lorenz (1996)
    # x: 状態変数
    # F: 強制項の値 (スカラー)
    # N: x の個数
    m_res = fill(0.0,N)
    
    for i in 3:N-1
        m_res[i] = - x[i] + (x[i+1] - x[i-2]) * x[i-1] + F
    end
    # i=1
    m_res[1] = - x[1] + (x[2] - x[N-1]) * x[N] + F
    # i=2
    m_res[2] = - x[2] + (x[3] - x[N]) * x[1] + F
    # i=N
    m_res[N] = - x[N] + (x[1] - x[N-2]) * x[N-1] + F
    
    return m_res
end

function L96_EU1(x,dt,F,N)  # Lorenz (1996) with explicit Euler scheme
    # x: 状態変数
    # dt: 時間ステップ
    # F: 強制項の値 (スカラー)
    # N: x の個数
    m_res = fill(0.0,N)
    
    m_res = x + dt * L96_tendency(x,F,N)
    
    return m_res
end

function L96_RK4(x,dt,F,N)  # Lorenz (1996) with the 4th-order Runge-Kutta scheme
    # x: 状態変数
    # dt: 時間ステップ
    # F: 強制項の値 (スカラー)
    # N: x の個数
    m_res = fill(0.0,N)
    
    k1 = L96_tendency(x,F,N)
    k2 = L96_tendency(x+0.5*dt*k1,F,N)
    k3 = L96_tendency(x+0.5*dt*k2,F,N)
    k4 = L96_tendency(x+dt*k3,F,N)
    
    m_res = x + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    
    return m_res
end

function L96_EU1_tangent(x,dt,N)  # The tangential operator for L96_EU1
    # x: 状態変数
    # dt: 時間ステップ
    # N: x の個数
    mt_res = fill(0.0,N,N)
    
    # 対角成分
    for i in 1:N
        mt_res[i,i] = (1.0 - dt)
    end
    # k-2
    for i in 3:N
        mt_res[i,i-2] = - dt * x[i-1]
    end
    # k-1
    for i in 3:N-1
        mt_res[i,i-1] = dt * (x[i+1] - x[i-2])
    end
    # k+1
    for i in 2:N-1
        mt_res[i,i+1] = dt * x[i-1]
    end
    # other terms
    mt_res[1,2] = dt * x[N]
    mt_res[N,1] = dt * x[N-1]
    mt_res[1,N] = dt * (x[2] - x[N-1])
    mt_res[2,1] = dt * (x[3] - x[N])
    mt_res[N,N-1] = dt * (x[1] - x[N-2])
    mt_res[1,N-1] = - dt * x[N]
    mt_res[2,N] = - dt * x[1]
    
    return mt_res
end

function L96_RK4_Fdx(x,N)  # The gradient for the forcing term (i.e., dx/dt = F <- this) in L96_RK4
    # x: 状態変数
    # dt: 時間ステップ
    # N: x の個数
    mt_res = fill(0.0,N,N)
    
    # 対角成分
    for i in 1:N
        mt_res[i,i] = -1.0
    end
    # k-2
    for i in 3:N
        mt_res[i,i-2] = -x[i-1]
    end
    # k-1
    for i in 3:N-1
        mt_res[i,i-1] = x[i+1] - x[i-2]
    end
    # k+1
    for i in 2:N-1
        mt_res[i,i+1] = x[i-1]
    end
    # other terms
    mt_res[1,2] = x[N]
    mt_res[N,1] = x[N-1]
    mt_res[1,N] = x[2] - x[N-1]
    mt_res[2,1] = x[3] - x[N]
    mt_res[N,N-1] = x[1] - x[N-2]
    mt_res[1,N-1] = -x[N]
    mt_res[2,N] = -x[1]
    
    return mt_res
end

function L96_RK4_tangent(x,dt,F,N)  # The tangential operator for L96_RK4
    # x: 状態変数
    # dt: 時間ステップ
    # F: 強制項の値 (スカラー)
    # N: x の個数
    mt_res = fill(0.0,N,N)
    dt6 = dt / 6.0
    dt2 = 0.5 * dt
    
    I_mat = Matrix{Float64}(I,nx,nx)  # 単位行列の利用

    k1 = L96_tendency(x,F,N)
    x2 = x + dt2 .* k1
    k2 = L96_tendency(x2,F,N)
    x3 = x + dt2 .* k2
    k3 = L96_tendency(x3,F,N)
    x4 = x + dt .* k3
    k4 = L96_tendency(x4,F,N)
    
    F1 = L96_RK4_Fdx(x,N)
    F2 = L96_RK4_Fdx(x2,N)
    F3 = L96_RK4_Fdx(x3,N)
    F4 = L96_RK4_Fdx(x4,N)
    
    F12 = (I_mat + dt2 .* F1) * F2
    F123 = (I_mat + dt2 .* F12) * F3
    F1234 = (I_mat + dt .* F123) * F4
    
    mt_res = I_mat + dt6 .* (F1 + 2.0 .* F12 + 2.0 .* F123 + F1234)

    return mt_res
end

function L96_TL_check(x,deltax,dt,nt,F,N,scheme)  # The TL check for L96_{EU1,RK4}_tangent
    # x: 状態変数
    # deltax: x に対する微小変数
    # dt: 時間ステップ
    # nt: 検証する時間ステップ数
    # F: 強制項の値 (スカラー)
    # N: x の個数
    # scheme: 検証するスキームの種類, "EU1" or "RK4"
    xt = x
    delx = deltax
    println("delx ",delx)
    
    if scheme == "EU1"
        println("Start TL check for EU1")
    elseif scheme == "RK4"
        println("Start TL check for RK4")
    end
    
    for i in 1:nt
        if scheme == "EU1"
            Lxt = L96_EU1(xt,dt,F,N)  # Time integration for x
            Lxtd = L96_EU1(xt+delx,dt,F,N)  # Time integration for x + deltax
            Mop = L96_EU1_tangent(xt,dt,N)  # TLM for x
        elseif scheme == "RK4"
            Lxt = L96_RK4(xt,dt,F,N)  # Time integration for x
            Lxtd = L96_RK4(xt+delx,dt,F,N)  # Time integration for x + deltax
            Mop = L96_RK4_tangent(xt,dt,F,N)  # TLM for x
        end
        Ldx = Lxtd - Lxt
        Mdx = Mop * delx
        println("TL check:i=",i,", (dL/dx)δx=",Ldx'*Ldx,", Mδx=",Mdx'*Mdx)
        
        # Updating x and deltax
        xt = Lxt
        delx = Mdx
        println("dx check ",delx' * delx)
    end
    
end

###########################
# calculation (main) part #
###########################

using Random
using Statistics

rng = MersenneTwister(1234)

nt = 15000
nx = 40
no = 40
n_ascyc = 5  # assimilation cycle interval for "t"
rand_flag = true  # observation including random noise
nto = div(nt,n_ascyc) + 1
inf_fact = 1.05  # inflation factor
eps = 0.00001  # small value for tangential linear matrix of the model operator
nspin = 3000

# Allocating
t = zeros(nt)
t_o = zeros(nto)
y_o = reshape(zeros(no,nto),no,nto)
x_t = reshape(zeros(nx,nt),nx,nt)
x_spin = reshape(zeros(nx,nspin),nx,nspin)
x_an = reshape(zeros(nx),nx,1)  # nx 行 1 列行列へ変換
x_f = reshape(zeros(nx,nt),nx,nt)
x_e = reshape(zeros(nx,nt),nx,nt)
errors = reshape(zeros(nx,nt),nx,nt)
Mop_linf = reshape(zeros(nx,nx),nx,nx)
Mop_fulf = reshape(zeros(nx,1),nx,1)
Mop_fult = reshape(zeros(nx,1),nx,1)
Mop_fule = reshape(zeros(nx,1),nx,1)
Pf = reshape(zeros(nx,nx),nx,nx)
Pa = reshape(zeros(nx,nx),nx,nx)
Ro = reshape(zeros(no,no),no,no)
Hop = reshape(zeros(no,nx),no,nx)
Kg = reshape(zeros(nx,no),nx,no)
d_innov = reshape(zeros(no,1),no,1)
Pftime = reshape(zeros(nx,nx,2),nx,nx,2)
Egtime = reshape(zeros(nt,1),nt,1)
Egftime = reshape(zeros(nt,1),nt,1)
L2Ntime = reshape(zeros(nt,1),nt,1)
L2Natime = reshape(zeros(nto,1),nto,1)
sigma_R = reshape(zeros(no,1),no,1)
evec_max = reshape(zeros(nx,nto),nx,nto)
lam_max = reshape(zeros(nto),nto)
delta_x = reshape(zeros(nx,1),nx,1)

# Setting parameters
xinit = fill(1.0,nx)
xinit[1] = xinit[1] + 0.1
dt = 0.01
delx = 0.1  # small departure for the TL check
F = 8.0  # default: 8
sigma_const_R = 1.0

I_mat = Matrix{Float64}(I,nx,nx)  # 単位行列の利用
Pf = (1.0 .* I_mat) + fill(20.0,nx,nx)
Ro = sigma_const_R * I_mat[1:no,1:no]

#Hop = I_mat
Hop = fill(0.0,no,nx)
for i in 1:no
    Hop[i,i] = 1.0
end
lam_max[1] = 1.0
evec_max = fill(0.0,nx,nto)
eps_inv = 1.0 / eps

# TL check (If you need to check the Mop, please activate the next line.)
#L96_TL_check(x_t[1:nx,1],delx.*fill(1.0,nx,1),dt,1000,F,nx,"RK4")

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


Pftime[1:nx,1:nx,1] = Pf


for i in 1:nt-1
    if mod(i-1,n_ascyc) == 0  # Entering the analysis processes
        #println("Enter hear")
        #-- Analysis
        L_inv = inv(Hop * Pf * Hop' + Ro)
        Kg = (Pf * Hop') * L_inv
        Pa = Pf - (Kg * Hop) * Pf
        Pa = inf_fact .* Pa  # covariance inflation
        if rand_flag == true
            sigma_R = randn(rng, Float64) * map( x->sqrt(x), diag(Ro) )  # Ro の対角成分を抽出し, 平方根をとる.
            y_o[1:no,div(i-1,n_ascyc)+1] = x_t[1:no,i] + sigma_R[1:no,1]
        else
            y_o[1:no,div(i-1,n_ascyc)+1] = x_t[1:no,i]
        end
        d_innov = y_o[1:no,div(i-1,n_ascyc)+1] - x_f[1:no,i]        
        x_an = x_f[1:nx,i] + Kg * d_innov
        t_o[div(i-1,n_ascyc)+1] = dt*(i-1)
        
        L2Natime[div(i-1,n_ascyc)+1,1] = sqrt((x_an - x_t[1:nx,i])' * (x_an - x_t[1:nx,i]) / nx)
    else
        Pa = Pf
        x_an = x_f[1:nx,i]
    end
    Mop_fult = L96_RK4(x_t[1:nx,i],dt,F,nx)
    Mop_fulf = L96_RK4(x_an,dt,F,nx)
    Mop_fule = L96_RK4(x_e[1:nx,i],dt,F,nx)
    Mop_linf = L96_RK4_tangent(x_an,dt,F,nx)
    #for j in 1:nx
    #    Mop_linf[1:nx,j] = eps_inv .* (L96_RK4(x_an + eps .* I_mat[1:nx,j],dt,F,nx) - L96_RK4(x_an,dt,F,nx))
    #end
    #-- Forecast
    t[i] = dt*(i-1)
    
    x_t[1:nx,i+1] = Mop_fult
    x_e[1:nx,i+1] = Mop_fule
    x_f[1:nx,i+1] = Mop_fulf
    Pf = (Mop_linf * Pa) * Mop_linf'

    Egftime[i+1] = tr(Pf) / nx
    Egtime[i+1] = tr(Pa) / nx
    delta_x = x_f[1:nx,i+1] - x_t[1:nx,i+1]
    L2Ntime[i+1] = sqrt((delta_x' * delta_x) / nx)
    
    if mod(i-1,n_ascyc) == 0  # Entering the analysis processes
        lambda = eigen(Symmetric(Mop_linf' * Mop_linf),nx:nx)
        lam_max[div(i-1,n_ascyc)+1] = lambda.values[1]  # 最大固有値
        evec_max[1:nx,div(i-1,n_ascyc)+1] = lambda.vectors[1:nx,1]  # 最大固有値の固有ベクトル
    end
    
end
t[nt] = t[nt-1] + dt

errors = map(x->abs(x), x_f - x_t)
Pftime[1:nx,1:nx,2] = Pf

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

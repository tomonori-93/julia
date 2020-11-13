# Lorenz (1963) モデルで同化・解析処理を行うシミュレーション
# Author: Satoki Tsujino (satoki_at_gfd-dennou.org)
# Date: 2020/11/01
# License: LGPL2.1
#############################
# functions for calculation #
#############################

function L63_tendency(x,cx,cy,cz)  # calculation of tendency in Lorenz (1963)
    # x: 状態変数
    # c{x,y,z}: 各時間発展方程式のパラメータ (定数) 
    m_res = fill(0.0,3)
    
    # i=1
    m_res[1] = - cx * (x[1] - x[2])
    # i=2
    m_res[2] = - x[2] - x[1] * x[3] + cy * x[1]
    # i=3
    m_res[3] = - cz * x[3] + x[1] * x[2]
    
    return m_res
end

function L63_EU1(x,dt,cx,cy,cz)  # Lorenz (1963) with explicit Euler scheme
    # x: 状態変数
    # dt: 時間ステップ
    # c{x,y,z}: 各時間発展方程式のパラメータ (定数) 
    m_res = fill(0.0,3)
    
    m_res = x + dt * L63_tendency(x,cx,cy,cz)
    
    return m_res
end

function L63_RK4(x,dt,cx,cy,cz)  # Lorenz (1963) with the 4th-order Runge-Kutta scheme
    # x: 状態変数
    # dt: 時間ステップ
    # c{x,y,z}: 各時間発展方程式のパラメータ (定数) 
    m_res = fill(0.0,3)
    
    k1 = L63_tendency(x,cx,cy,cz)
    k2 = L63_tendency(x+0.5*dt*k1,cx,cy,cz)
    k3 = L63_tendency(x+0.5*dt*k2,cx,cy,cz)
    k4 = L63_tendency(x+dt*k3,cx,cy,cz)
    
    m_res = x + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    
    return m_res
end

function L63_tangent(x,dt,cx,cy,cz)  # The tangential operator for L63_full
    # x: 状態変数
    # dt: 時間ステップ
    # c{x,y,z}: 各時間発展方程式のパラメータ (定数) 
    mt_res = fill(0.0,3,3)
    
    mt_res[1,1] = 1.0 - cx * dt
    mt_res[2,2] = 1.0 - dt
    mt_res[3,3] = 1.0 - cz * dt
    mt_res[1,2] = cx * dt
    mt_res[2,1] = (cy - x[3]) * dt
    mt_res[2,3] = - x[1] * dt
    mt_res[3,1] = x[2] * dt
    mt_res[3,2] = x[1] * dt
    
    return mt_res
end

function L63_TL_check(x,deltax,dt,nt,cx,cy,cz)  # The TL check for L63_tangent
    # x: 状態変数
    # deltax: x に対する微小変数
    # dt: 時間ステップ
    # nt: 検証する時間ステップ数
    # c{x,y,z}: 各時間発展方程式のパラメータ (定数) 
    xt = x
    delx = deltax
    println("delx ",delx)
    
    for i in 1:nt
        Lxt = L63_EU1(xt,dt,cx,cy,cz)  # Time integration for x
        Lxtd = L63_EU1(xt+delx,dt,cx,cy,cz)  # Time integration for x + deltax
        Ldx = Lxtd - Lxt
        Mop = L63_tangent(xt,dt,cx,cy,cz)  # TLM for x
        Mdx = Mop * delx
        println("TL check:i=",i,", (dL/dx)δx=",Ldx'*Ldx,", Mδx=",Mdx'*Mdx)
        
        # Updating x and deltax
        xt = Lxt
        delx = Mdx
        println("dx check ",delx' * delx)
    end
    
end

# Lorenz (1963) モデルで各時刻に同化・解析処理を挟むシミュレーション
using PyPlot
using Random
using LinearAlgebra

rng = MersenneTwister(1234)

nt = 1000000
nx = 3
no = 3
n_ascyc = 250000  # assimilation cycle interval for "t"
rand_flag = false # observation including random noise
nto = div(nt,n_ascyc) + 1
inf_fact = 1.0  # inflation factor
eps = 0.0005  # small value for tangential linear matrix of the model operator

# Allocating
t = zeros(nt)
t_o = zeros(nto)
y_o = reshape(zeros(no,nto),no,nto)
x_t = reshape(zeros(nx,nt),nx,nt)
x_an = reshape(zeros(nx),nx,1)  # nx 行 1 列行列へ変換
x_an_lin = reshape(zeros(nx),nx,1)  # nx 行 1 列行列へ変換
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
Pftime = reshape(zeros(nx,nt),nx,nt)
sigma_R = reshape(zeros(no,1),no,1)
evec_max = reshape(zeros(nx,nt),nx,nt)
lam_max = reshape(zeros(nt),nt)
Egtime = reshape(zeros(nt,1),nt,1)
L2Ntime = reshape(zeros(nt,1),nt,1)
delta_x = reshape(zeros(nx,1),nx,1)

# Setting parameters
xinit = [10.0, 20.0, 30.0]
dt = 0.00001
delx = 0.1  # small departure for the TL check
k = 10.0
r = 23.0
b = 8.0/3.0
I_mat = Matrix{Float64}(I,nx,nx)  # 単位行列の利用
Pf = [20.0 20.0 20.0
      20.0 20.0 20.0
      20.0 20.0 20.0]
Ro = [1.0 0.0 0.0
      0.0 1.0 0.0
      0.0 0.0 1.0]  # 行ベクトル [a b], 列ベクトル [a, b]
x_t[1:nx,1] = xinit[1:nx]
x_f[1:nx,1] = x_t[1:nx,1] + 0.1 * xinit[1:nx]
x_e[1:nx,1] = x_f[1:nx,1]
Hop = I_mat
lam_max[1] = 1.0
evec_max = fill(1.0,nx,nt)
eps_inv = 1.0 / eps

# TL check (If you need to check the Mop, please activate the next line.)
# L63_TL_check(x_t[1:nx,1],delx.*fill(1.0,nx,1),dt,100,k,r,b)

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
        x_an_lin = x_an
        t_o[div(i-1,n_ascyc)+1] = dt*(i-1)
    else
        Pa = Pf
        x_an = x_f[1:nx,i]
        x_an_lin = x_an
    end

    Mop_fult = L63_EU1(x_t[1:nx,i],dt,k,r,b)
    Mop_fulf = L63_EU1(x_an,dt,k,r,b)
    Mop_fule = L63_EU1(x_e[1:nx,i],dt,k,r,b)
    Mop_linf = L63_tangent(x_an_lin,dt,k,r,b)
    #for j in 1:nx
    #    Mop_linf[1:nx,j] = eps_inv .* (L63_RK4(x_an_lin + eps .* I_mat[1:nx,j],dt,k,r,b) - L63_RK4(x_an_lin,dt,k,r,b))
    #end

    #-- Forecast
    t[i] = dt*(i-1)
    
    x_t[1:nx,i+1] = Mop_fult
    x_e[1:nx,i+1] = Mop_fule
    x_f[1:nx,i+1] = Mop_fulf
    Pf = inf_fact * Mop_linf * Pa * Mop_linf'
    #println(Pf[1,1], ",", Pf[2,2], ",", Pf[3,3])
    Pftime[1:nx,i+1] = diag(Pf)
    Egtime[i+1] = tr(Pf) / nx
    delta_x = x_f[1:nx,i+1] - x_t[1:nx,i+1]
    L2Ntime[i+1] = sqrt(delta_x' * delta_x) / nx
    lambda = eigen(Symmetric(Mop_linf' * Mop_linf),nx:nx)
    lam_max[i+1] = lambda.values[1]  # 最大固有値
    evec_max[1:nx,i+1] = lambda.vectors[1:nx,1]  # 最大固有値の固有ベクトル
    
    if mod(i-1,n_ascyc) == 0
        println("Max Lam = ", i, ", ", lam_max[i+1])
        println("Max Evec = ", evec_max[1:nx,i+1])
    end

end
t[nt] = t[nt-1] + dt

errors = map(x->abs(x), x_f - x_t)

#-- Drawing
##########
#  Plot  #
##########
fig = figure("pyplot_majorminor",figsize=(7,5))
##p = plot(t[1:nt],x_t[1,1:nt],color="black",label="Perfect")
##p = plot(t[1:nt],x_e[1,1:nt],color="blue",label="False")
##p = plot(t[1:nt],x_f[1,1:nt],color="red",label="Assim")
##p = plot(t_o[1:div(nt-1,n_ascyc)+1],y_o[1,1:div(nt-1,n_ascyc)+1],linestyle="",marker="o",color="black",label="Obs")
#p = plot(x_t[1,1:nt],x_t[3,1:nt],color="black",label="Perfect")
#p = plot(x_e[1,1:nt],x_e[3,1:nt],color="blue",label="False")
#p = plot(x_f[1,1:nt],x_f[3,1:nt],color="red",label="Assim")
#p = plot(y_o[1,1:div(nt-1,n_ascyc)+1],y_o[3,1:div(nt-1,n_ascyc)+1],linestyle="",marker="o",color="black",label="Obs")
p = plot(t[1:nt],Egtime[1:nt]*0.01,color="black",label="Tr(Pf)")
p = plot(t[1:nt],L2Ntime[1:nt],color="red",label="L2norm")
p = plot(t_o[1:div(nt-1,n_ascyc)+1],fill(0.0,div(nt-1,n_ascyc)+1),linestyle="",marker="o",color="black",label="Obs")
###p = plot(t[1:nt],lam_max[1:nt],color="green",label="Lambda_max")
legend()
ax = gca()
xlabel("Time")
ylabel("ΔX")
##ylabel("X^T")
#xlabel("X")
#ylabel("Z")
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

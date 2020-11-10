# Lorenz (1963) モデルで各時刻に同化・解析処理を挟むシミュレーション
using PyPlot
using Random
using LinearAlgebra

rng = MersenneTwister(1234)

nt = 10000
nx = 3
no = 3
n_ascyc = 2500  # assimilation cycle interval for "t"
rand_flag = false # observation including random noise
nto = div(nt,n_ascyc) + 1
inf_fact = 1.0  # inflation factor

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
x_prev = reshape(zeros(nx),nx,1)
x_eprev = reshape(zeros(nx),nx,1)
x_anprev = reshape(zeros(nx),nx,1)
Hop = reshape(zeros(no,nx),no,nx)
Kg = reshape(zeros(nx,no),nx,no)
Pftime = reshape(zeros(nx,nt),nx,nt)
sigma_R = reshape(zeros(no,1),no,1)
evec_max = reshape(zeros(nx,nt),nx,nt)
lam_max = reshape(zeros(nt),nt)

# Setting parameters
xinit = [10.0, 20.0, 30.0]
dt = 0.001
k = 10.0
r = 23.0
b = 8.0/3.0
Pf = [20.0 20.0 20.0
      20.0 20.0 20.0
      20.0 20.0 20.0]
Ro = [1.0 0.0 0.0
      0.0 1.0 0.0
      0.0 0.0 1.0]  # 行ベクトル [a b], 列ベクトル [a, b]
x_t[1:nx,1] = xinit[1:nx]
x_f[1:nx,1] = x_t[1:nx,1] + 0.1 * xinit[1:nx]
x_e[1:nx,1] = x_f[1:nx,1]
Hop = [1.0 0.0 0.0
       0.0 1.0 0.0
       0.0 0.0 1.0]  # 行ベクトル [a b], 列ベクトル [a, b]
lam_max[1] = 1.0
evec_max = fill(1.0,nx,nt)

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
        #println("x_t = ", x_t[1:nx,i])
        #println("x_f = ", x_f[1:nx,i])
        #println("x_an = ", x_an)
        println("d_innov = ", d_innov)
        #println("Kg = ", Kg)
        #println("Pf = ", Pf)
        #println("Pa = ", Pa)
        #println("L_inv = ", L_inv)
        println("incre = ", Kg * d_innov)
        t_o[div(i-1,n_ascyc)+1] = dt*(i-1)
    else
        Pa = Pf
        x_an = x_f[1:nx,i]
        x_an_lin = x_an
    end
    Mop_fult = [ (1.0 - k * dt) * x_t[1,i] + dt * k * x_t[2,i]
                 dt * r * x_t[1,i] + (1.0 - dt) * x_t[2,i] - dt * x_t[1,i] * x_t[3,i]
                 dt * x_t[1,i] * x_t[2,i] + (1.0 - b * dt) * x_t[3,i] ]
    Mop_fulf = [ (1.0 - k * dt) * x_an[1] + dt * k * x_an[2]
                 dt * r * x_an[1] + (1.0 - dt) * x_an[2] - dt * x_an[1] * x_an[3]
                 dt * x_an[1] * x_an[2] + (1.0 - b * dt) * x_an[3] ]
    Mop_fule = [ (1.0 - k * dt) * x_e[1,i] + dt * k * x_e[2,i]
                 dt * r * x_e[1,i] + (1.0 - dt) * x_e[2,i] - dt * x_e[1,i] * x_e[3,i]
                 dt * x_e[1,i] * x_e[2,i] + (1.0 - b * dt) * x_e[3,i] ]
    Mop_linf = [ (1.0 - k * dt) (dt * k) 0.0
                 (dt * r - dt * x_an_lin[3]) (1.0 - dt) (- dt * x_an_lin[1])
                 (dt * x_an_lin[2]) (dt * x_an_lin[1]) (1.0 - b * dt) ]
    
    #-- Forecast
    t[i] = dt*(i-1)
    
    x_t[1:nx,i+1] = Mop_fult
    x_e[1:nx,i+1] = Mop_fule
    x_f[1:nx,i+1] = Mop_fulf
    Pf = inf_fact * Mop_linf * Pa * Mop_linf'
    #println(Pf[1,1], ",", Pf[2,2], ",", Pf[3,3])
    Pftime[1:nx,i+1] = diag(Pf)
    lambda = eigen(Symmetric(Mop_linf' * Mop_linf),nx:nx)
    lam_max[i+1] = lambda.values[1]  # 最大固有値
    evec_max[1:nx,i+1] = lambda.vectors[1:nx,1]  # 最大固有値の固有ベクトル
    
    if mod(i-1,n_ascyc) == 0
        println("Max Lam = ", i, ", ", lam_max[i+1])
        println("Max Evec = ", evec_max[1:nx,i+1])
#    println("lambda", lambda[nx,i+1], evec_max[1:nx,i+1])
#    println(diag(Pf)[1])
#    println("test",eigmax(Mop_linf' * Mop_linf))
#    println("test",eigvals(Mop_linf' * Mop_linf))
    end
    if (i-1)%100 == 0
    println("Pf, ",i,", ",Pf)
    #println("M",Mop_linf * Pa * Mop_linf')
    end

end
t[nt] = t[nt-1] + dt

errors = map(x->abs(x), x_f - x_t)

#-- Drawing
##########
#  Plot  #
##########
#nt = 20
fig = figure("pyplot_majorminor",figsize=(7,5))
p = plot(t[1:nt],x_t[1,1:nt],color="black",label="Perfect")
p = plot(t[1:nt],x_e[1,1:nt],color="blue",label="False")
p = plot(t[1:nt],x_f[1,1:nt],color="red",label="Assim")
p = plot(t_o[1:div(nt-1,n_ascyc)+1],y_o[1,1:div(nt-1,n_ascyc)+1],linestyle="",marker="o",color="black",label="Obs")
#p = plot(x_t[1,1:nt],x_t[3,1:nt],color="black",label="Perfect")
#p = plot(x_e[1,1:nt],x_e[3,1:nt],color="blue",label="False")
#p = plot(x_f[1,1:nt],x_f[3,1:nt],color="red",label="Assim")
#p = plot(y_o[1,1:div(nt-1,n_ascyc)+1],y_o[3,1:div(nt-1,n_ascyc)+1],linestyle="",marker="o",color="black",label="Obs")
###p = plot(t[n_ascyc-20:n_ascyc+20],errors[1,n_ascyc-20:n_ascyc+20],color="black",label="V_x")
###p = plot(t[n_ascyc-20:n_ascyc+20],errors[2,n_ascyc-20:n_ascyc+20],color="blue",label="V_y")
###p = plot(t[n_ascyc-20:n_ascyc+20],errors[3,n_ascyc-20:n_ascyc+20],color="red",label="V_z")
###p = plot(t_o[2],fill(0.0,1),linestyle="",marker="o",color="black",label="Obs")
###p = plot(t[1:nt],lam_max[1:nt],color="green",label="Lambda_max")
legend()
ax = gca()
xlabel("Time")
###ylabel("ΔX")
ylabel("X^T")
#xlabel("X")
#ylabel("Z")
grid("on")
PyPlot.title("Lorenz63")
#PyPlot.title("Phase")

###########################
#  Set the tick interval  #
###########################
###Mx = matplotlib.ticker.MultipleLocator(0.1) # Define interval of major ticks
Mx = matplotlib.ticker.MultipleLocator(1.0) # Define interval of major ticks
#Mx = matplotlib.ticker.MultipleLocator(4.0) # Define interval of major ticks
#f = matplotlib.ticker.FormatStrFormatter("%1.2f") # Define format of tick labels
ax.xaxis.set_major_locator(Mx) # Set interval of major ticks
#ax.xaxis.set_major_formatter(f) # Set format of tick labels

###mx = matplotlib.ticker.MultipleLocator(0.02) # Define interval of minor ticks
mx = matplotlib.ticker.MultipleLocator(0.1) # Define interval of minor ticks
#mx = matplotlib.ticker.MultipleLocator(1.0) # Define interval of minor ticks
ax.xaxis.set_minor_locator(mx) # Set interval of minor ticks

###My = matplotlib.ticker.MultipleLocator(2.0) # Define interval of major ticks
My = matplotlib.ticker.MultipleLocator(2.0) # Define interval of major ticks
#My = matplotlib.ticker.MultipleLocator(4.0) # Define interval of major ticks
ax.yaxis.set_major_locator(My) # Set interval of major ticks

###my = matplotlib.ticker.MultipleLocator(0.5) # Define interval of minor ticks
my = matplotlib.ticker.MultipleLocator(0.5) # Define interval of minor ticks
#my = matplotlib.ticker.MultipleLocator(1.0) # Define interval of minor ticks
ax.yaxis.set_minor_locator(my) # Set interval of minor ticks

#########################
#  Set tick dimensions  #
#########################
#ax.xaxis.set_tick_params(which="major",length=5,width=2,labelsize=10)
#ax.xaxis.set_tick_params(which="minor",length=5,width=2)

fig.canvas.draw() # Update the figure
gcf() # Needed for IJulia to plot inline
###savefig("Error-time.pdf")
savefig("X-time.pdf")
#savefig("X-Z_lorenz.pdf")

using PyPlot
using Random
using LinearAlgebra

rng = MersenneTwister(1234)

nt = 100000
nx = 2
no = 1
n_ascyc = 45000  # assimilation cycle interval for "t"
rand_flag = false  # observation including random noise
nto = div(nt,n_ascyc) + 1

# Initializing
xinit = [1.0, 1.0]
t = zeros(nt)
t_o = zeros(nto)
y_o = reshape(zeros(no,nto),no,nto)
x_t = reshape(zeros(nx,nt),nx,nt)
x_a = reshape(zeros(nx,nt),nx,nt)  # nx 行 nt 列で定義: x_a[行,列]
x_an = reshape(zeros(nx),nx,1)  # nx 行 1 列行列へ変換
x_f = reshape(zeros(nx,nt),nx,nt)
x_e = reshape(zeros(nx,nt),nx,nt)
Mop = reshape(zeros(nx,nx),nx,nx)
Pf = reshape(zeros(nx,nx),nx,nx)
Pa = reshape(zeros(nx,nx),nx,nx)
#Ro = reshape(zeros(no,no),no,no)
x_prev = reshape(zeros(nx),nx,1)
x_eprev = reshape(zeros(nx),nx,1)
x_anprev = reshape(zeros(nx),nx,1)
Hop = reshape(zeros(no,nx),no,nx)
Kg = reshape(zeros(nx,no),nx,no)
Pftime = reshape(zeros(nx,nt),nx,nt)

# Setting parameters
dt = 0.0001
k = 40.0
ksqrt = sqrt(k)
Pf = fill(1.0,nx,nx)
Ro = 0.001
x_t[1:nx,1] = xinit[1:nx]
x_a[1:nx,1] = xinit[1:nx]
x_f[1:nx,1] = x_t[1:nx,1] + [-0.1, 0.2]
x_e[1:nx,1] = x_f[1:nx,1]
Mop = [1.0 dt
       -k*dt 1.0]
Hop = [1.0 0.0]  # 行ベクトル [a b], 列ベクトル [a, b]
Pftime[1:nx,1] = diag(Pf)  # 対角成分のみ抽出

for i in 1:nt-1
    if mod(i-1,n_ascyc) == 0  # Entering the analysis processes
        #-- Analysis
        L_inv = 1.0 / (Pf[1,1]+Ro)
        Kg = (Pf * Hop') * L_inv
        Pa = Pf - (Kg * Hop) * Pf
        if rand_flag == true
            sigma_R = randn(rng, Float64) * sqrt(Ro)
            y_o[1,div(i-1,n_ascyc)+1] = x_t[1,i] + sigma_R
        else
            y_o[1,div(i-1,n_ascyc)+1] = x_t[1,i]
        end
        d_innov = y_o[1,div(i-1,n_ascyc)+1] - x_f[1,i]
        x_an = x_f[1:nx,i] + Kg[1:nx,1] * d_innov
        t_o[div(i-1,n_ascyc)+1] = dt*(i-1)
    else
        Pa = Pf
        x_an = x_f[1:nx,i]
    end
    
    #-- Forecast
    t[i] = dt*(i-1)
    x_prev = x_t[1:nx,i]
    x_eprev = x_e[1:nx,i]
    x_anprev = x_an[1:nx]
    
    x_t[1:nx,i+1] = Mop * x_prev
    x_e[1:nx,i+1] = Mop * x_eprev
    x_f[1:nx,i+1] = Mop * x_anprev
    Pf = Mop * Pa * Mop'
    Pftime[1:nx,i+1] = diag(Pf)

    #-- Analytical solution
    x_a[1,i] = xinit[1] * cos(ksqrt*t[i]) + (xinit[2] / ksqrt) * sin(ksqrt*t[i])
    x_a[2,i] = - ksqrt * xinit[2] * sin(ksqrt*t[i]) + xinit[2] * cos(ksqrt*t[i])
end
t[nt] = t[nt-1] + dt
x_a[1,nt] = xinit[1] * cos(ksqrt*t[nt]) + (xinit[2] / ksqrt) * sin(ksqrt*t[nt])
x_a[2,nt] = - ksqrt * xinit[2] * sin(ksqrt*t[nt]) + xinit[2] * cos(ksqrt*t[nt])


#-- Drawing
##########
#  Plot  #
##########
#nt = 20
fig = figure("pyplot_majorminor",figsize=(7,5))
p = plot(t[1:nt],x_t[1,1:nt],color="black",label="Perfect")
#p = plot(t[1:nt],x_a[1,1:nt],label="Exact")
p = plot(t[1:nt],x_e[1,1:nt],color="blue",label="False")
p = plot(t[1:nt],x_f[1,1:nt],color="red",label="Assim")
p = plot(t_o[1:div(nt-1,n_ascyc)+1],y_o[1,1:div(nt-1,n_ascyc)+1],linestyle="",marker="o",color="black",label="Obs")
legend()
ax = gca()
xlabel("Time")
ylabel("X^T")
grid("on")
PyPlot.title("Linear-Osci.")

###########################
#  Set the tick interval  #
###########################
Mx = matplotlib.ticker.MultipleLocator(1) # Define interval of major ticks
f = matplotlib.ticker.FormatStrFormatter("%1.2f") # Define format of tick labels
ax.xaxis.set_major_locator(Mx) # Set interval of major ticks
ax.xaxis.set_major_formatter(f) # Set format of tick labels

mx = matplotlib.ticker.MultipleLocator(0.1) # Define interval of minor ticks
ax.xaxis.set_minor_locator(mx) # Set interval of minor ticks

My = matplotlib.ticker.MultipleLocator(0.5) # Define interval of major ticks
ax.yaxis.set_major_locator(My) # Set interval of major ticks

my = matplotlib.ticker.MultipleLocator(0.1) # Define interval of minor ticks
ax.yaxis.set_minor_locator(my) # Set interval of minor ticks

#########################
#  Set tick dimensions  #
#########################
#ax.xaxis.set_tick_params(which="major",length=5,width=2,labelsize=10)
#ax.xaxis.set_tick_params(which="minor",length=5,width=2)

fig.canvas.draw() # Update the figure
gcf() # Needed for IJulia to plot inline
savefig("X-time_linear.pdf")

#-- Drawing
##########
#  Plot  #
##########
#nt = 20
fig = figure("pyplot_majorminor",figsize=(7,5))
q = plot(t[1:nt],Pftime[2,1:nt],color="black",label="AssimB")
#p = plot(t[1:nt],x_a[1,1:nt],label="Exact")
#p = plot(t[1:nt],x_e[1,1:nt],color="blue",label="False")
#p = plot(t[1:nt],x_f[1,1:nt],color="red",label="Assim")
#p = plot(t_o[1:div(nt-1,n_ascyc)+1],y_o[1,1:div(nt-1,n_ascyc)+1],linestyle="",marker="o",color="black",label="Obs")
legend()
ax = gca()
xlabel("Time")
ylabel("X^T")
grid("on")
PyPlot.title("Time series")

###########################
#  Set the tick interval  #
###########################
Mx = matplotlib.ticker.MultipleLocator(1) # Define interval of major ticks
f = matplotlib.ticker.FormatStrFormatter("%1.2f") # Define format of tick labels
ax.xaxis.set_major_locator(Mx) # Set interval of major ticks
ax.xaxis.set_major_formatter(f) # Set format of tick labels

mx = matplotlib.ticker.MultipleLocator(0.1) # Define interval of minor ticks
ax.xaxis.set_minor_locator(mx) # Set interval of minor ticks

My = matplotlib.ticker.MultipleLocator(0.5) # Define interval of major ticks
ax.yaxis.set_major_locator(My) # Set interval of major ticks

my = matplotlib.ticker.MultipleLocator(0.1) # Define interval of minor ticks
ax.yaxis.set_minor_locator(my) # Set interval of minor ticks

#########################
#  Set tick dimensions  #
#########################
#ax.xaxis.set_tick_params(which="major",length=5,width=2,labelsize=10)
#ax.xaxis.set_tick_params(which="minor",length=5,width=2)

fig.canvas.draw() # Update the figure
#gcf() # Needed for IJulia to plot inline
savefig("X-time.pdf")

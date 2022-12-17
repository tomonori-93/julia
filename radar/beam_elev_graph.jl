#-- Constants
rad_e = 6.371e6  # Earth radius (in meridional) [m]
d2r = π / 180.0  # degree to radian
r2d = 180.0 / π  # radian to degree

#--- Setting radar parameters
nr = 100  # sampling number in radius
dr = 2.0e3  # sampling grid interval in radius [m]
na = 14   # number of the elevation angle
phi = reshape([ 0.0 0.1 0.3 0.5 0.7 1.0 2.0 3.0 4.0 5.0 10.0 15.0 20.0 25.0 ],na,1)  # elevation angles [degree]
lidx = [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ]  # line indices for each radar beam
h_radar = 350.0  # the above surface level [m] of the radar

#--- Initialize
rl = reshape(collect(Float32, 0.0:dr:dr*(nr-1)),nr,1)  # radar line
I_nr = reshape(fill(1.0,nr),nr,1)  # dummy vector for nr
I_na = reshape(fill(1.0,na),na,1)  # dummy vector for na
rl_mat = rl * I_na'  # = rl(nr,na)
phi_mat = I_nr * phi'   # = phi(nr,na)

#println(size(rl_mat),size(phi_mat))
#--- Calculation of radial and vertical points for (r,z)

x .= rad_e .+ h_radar .+ rl_mat .* cos.(π ./ 2 .- phi_mat .* d2r)
y .= rl_mat .* sin.(π ./ 2 .- phi_mat .* d2r)
theta .= atan.( y ./ x)

r .= x .* theta ./ cos.(theta)
z .= x ./ cos.(theta) .- rad_e

#-- Drawing
##########
#  Plot  #
##########
using PyPlot

#--- Setting parameters
dfact = 1.0e-3

rc("font", family="IPAPGothic")
fig = figure("pyplot_majorminor",figsize=(7,5))

ax = gca()

xlabel("Radius (km)")
ylabel("Height (km)")
#plt.colorbar(cp)
grid("on")

xlim(0.0,dr*(nr-1)*dfact)
ylim(0.0,10.0e3*dfact)

PyPlot.title("Radar beams")

plot(r.*dfact,z.*dfact)

fig.canvas.draw() # Update the figure
gcf() # Needed for IJulia to plot inline
fname = "test.png"
savefig(fname)

# Pull in data
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather.csv"
weather <- read.csv(file_path, header = TRUE)

# Simulations using Dew Temp
weather$sim1 = -2.901738 + weather$VPmax..mbar.*0.947056
weather$sim2 = -2.274677  + weather$VPmax..mbar.*0.815350 + weather$Tdew..degC.*0.236818
weather$sim3 = 45.09 + 1.354*(243.5*log(weather$VPmax..mbar./100))/(17.62-log(weather$VPmax..mbar./100))
weather$sim4 = 44.91  + 1.348*(243.5*log(weather$VPmax..mbar./100))/(17.62-log(weather$VPmax..mbar./100)) + 2.793e-04*weather$VPmax..mbar.*weather$Tdew..degC.

# Plots
plotLim = 2000
plot(weather$T..degC.[1:plotLim], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Tdew Lin")
points(weather$sim1[1:plotLim], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[1:plotLim], col = "green", pch = 18, cex = 1.5)
legend("topleft", legend = c("Temp", "VPMax", "VPMax+Tdew"), 
       col = c("red", "blue", "green"), lwd = 2)
plotStart=30000
plotLim = 2000
ymin = min(weather$T..degC.[(plotStart):(plotStart+plotLim-1)],weather$sim1[1:plotLim],weather$sim1[(plotStart):(plotStart+plotLim-1)],weather$sim2[(plotStart):(plotStart+plotLim-1)])
ymax = max(weather$T..degC.[(plotStart):(plotStart+plotLim-1)],weather$sim1[1:plotLim],weather$sim1[(plotStart):(plotStart+plotLim-1)],weather$sim2[(plotStart):(plotStart+plotLim-1)])
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Tdew Lin",ylim=c(ymin,ymax))
points(weather$sim1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 18, cex = 1.5)
legend("topleft", legend = c("Temp", "VPMax", "VPMax+Tdew"), 
       col = c("red", "blue", "green"), lwd = 2)

plotLim = 20
plot(weather$T..degC.[1:plotLim], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Tdew NonLin")
points(weather$sim3[1:plotLim], col = "blue", pch = 17, cex = 1.5)
points(weather$sim4[1:plotLim], col = "green", pch = 18, cex = 1.5)
legend("topleft", legend = c("Temp", "VPMax", "VPMax+Tdew"), 
       col = c("red", "blue", "green"), lwd = 2)

################################################################################

# Common variables
library(tswge)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
f1 = 1/(6*24)
f2 = 1/(6*24*366)
t = 1:(6*24*366)
n = 6*24*366
freq=c(f1,f2)

b0_inc = 5
A1_mult = 0.1
A2_mult = 1.2
var_mult = 0.6
sn = 2
b0 = mean(weather$VPmax..mbar.) + b0_inc
A1 = -3.552730 * A1_mult
phi1 = 5.356089 
A2 = -7.932237 * A2_mult
phi2 = -0.379505
psi=c(phi1,phi2+2*pi) # 5.356089 -0.379505
resid = weather$VPmax..mbar. - b0 - A1 * cos(2 * pi * f1 * t + phi1) - A2 * cos(2 * pi * f2 * t + phi2)
coef=c(A1,A2) 
est=est.ar.wge(resid,p=1,method='burg') # Trying p = 1
phi = est$phi
vara = est$avar * var_mult
y_gen_vpmax = gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=sn)

b0_inc = 0
A1_mult = 1.0
A2_mult = 1.0
var_mult = 1.0
sn = sn
b0 = mean(weather$Tdew..degC.) + b0_inc
A1 = 0.124165 * A1_mult
A2 = 6.557257 * A2_mult
phi1 = -1.181011 + pi # 1.960582
phi2 = 15.013296 - 4*pi # 2.446925
resid = weather$Tdew..degC. - b0 - A1 * cos(2 * pi * f1 * t + phi1) - A2 * cos(2 * pi * f2 * t + phi2)
plotLim = 6*24*366 
plotStart = 1
y_pred=A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2) + mean(weather$Tdew..degC.)
est=est.ar.wge(resid,p=1,method='burg') # trying p = 1
coef=c(A1,A2) 
psi=c(phi1,phi2+2*pi) 
phi = est$phi
vara = est$avar * var_mult
y_gen_tdew = gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=sn)

cor(weather$VPmax..mbar.,weather$Tdew..degC.) # 0.7126456
cor(y_gen_vpmax,y_gen_tdew) 
diff = cor(weather$VPmax..mbar.,weather$Tdew..degC.) - cor(y_gen_vpmax,y_gen_tdew)
alpha <- diff * (1-abs(diff))
y_gen_tdew_adj <- (1 - alpha) * y_gen_tdew + alpha * y_gen_vpmax
y_gen_vpmax_adj <- (1 - alpha) * y_gen_vpmax + alpha * y_gen_tdew_adj
y_gen_vpmax_adj[y_gen_vpmax_adj<0] = 0
cor(y_gen_vpmax_adj, y_gen_tdew_adj) 
# TODO: Maybe use a while loop to modify alpha

plotLim = 6*24*366 
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'VPMax')
points(y_gen_vpmax_adj[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 19, cex = 1.5)
plotLim = 2000
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'VPMax', ylim=c(0,15))
points(y_gen_vpmax_adj[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 19, cex = 1.5)
plotLim = 2000
plotStart = 30000
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'VPMax', ylim=c(10,45))
points(y_gen_vpmax_adj[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 19, cex = 1.5)

plotLim = 6*24*366 
plotStart = 1
plot(weather$Tdew..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'Tdew')
points(y_gen_tdew_adj[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 19, cex = 1.5)
plotLim = 2000 
plotStart = 1
plot(weather$Tdew..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'Tdew')
points(y_gen_tdew_adj[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 19, cex = 1.5)
# Initial stuff
library(tswge)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
f1 = 1/(6*24)
f2 = 1/(6*24*366)
t = 1:(6*24*366)
n = 6*24*366
freq=c(f1,f2)

b0 = mean(weather$VPmax..mbar.) + 10
A1_mult = 0.3
A2_mult = 1.6
A1 = -3.552730 * A1_mult
phi1 = 5.356089 
A2 = -7.932237 * A2_mult
phi2 = -0.379505
psi=c(phi1,phi2+2*pi) # 5.356089 -0.379505
resid = weather$VPmax..mbar. - b0 - A1 * cos(2 * pi * f1 * t + phi1) - A2 * cos(2 * pi * f2 * t + phi2)
coef=c(A1,A2) # -3.552730 -7.932237
est=est.ar.wge(resid,p=1,method='burg') # Trying p = 1
phi = est$phi
vara = est$avar
y_gen1 = gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=0)

plotLim = 6*24*366 
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'VPMax')
points(y_gen1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 19, cex = 1.5)
plotLim = 2000
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'VPMax', ylim=c(0,15))
points(y_gen1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 19, cex = 1.5)
plotts.wge(y_gen1[(plotStart):(plotStart+plotLim-1)], main = 'VPMax')
# Is there a way to increase daily variance without increasing 10 min variance?
# I'll probably need to do max(0,val), so no negs

###########################################################################

b0 = mean(log(weather$VPdef..mbar.+1)) + 0.0
A1_mult = 0.05
A2_mult = 0.4
var_mult = 0.6
A1 = 0.55900 * A1_mult
A2 = 0.642318 * A2_mult
phi1 = -0.942353 + pi # 2.199243
phi2 = 0.170724 + pi # 3.312313
y_pred=A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2) + b0 
resid = log(weather$VPdef..mbar.+1) - y_pred
est=est.ar.wge(resid,p=1,method='burg')
coef=c(A1,A2) 
psi=c(phi1,phi2+2*pi) 
phi = est$phi
vara = est$avar * var_mult
y_gen2 = pmax(gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=0),0)
plotLim = 6*24*366 
plotStart = 1
plot(weather$VPdef..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'VPDef')
points((exp(y_gen2)-1)[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 19, cex = 1.5)
plotLim = 2000 
plotStart = 1
plot(weather$VPdef..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'VPDef')
points((exp(y_gen2)-1)[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 19, cex = 1.5)
plotts.wge((exp(y_gen2)-1)[(plotStart):(plotStart+plotLim-1)], main = 'VPDef')
# I feel like this is a little closer

################################################################################

# Look at correlation
cor(weather$VPmax..mbar.,weather$VPdef..mbar.) # 0.860394
cor(y_gen1,y_gen2) # 0.4590828

alpha <- abs(cor(y_gen1,y_gen2)-cor(weather$VPmax..mbar.,weather$VPdef..mbar.))*0.10
y_gen2_adj <- (1 - alpha) * y_gen2 + alpha * y_gen1
y_gen1_adj <- (1 - alpha) * y_gen1 + alpha * y_gen2_adj
cor(y_gen1_adj, y_gen2_adj) # 0.8488848

plotLim = 6*24*366 
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'VPMax')
points(y_gen1_adj[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 19, cex = 1.5)
plotLim = 2000
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'VPMax', ylim=c(0,15))
points(y_gen1_adj[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 19, cex = 1.5)
plotts.wge(y_gen1_adj[(plotStart):(plotStart+plotLim-1)], main = 'VPMax')

plotLim = 6*24*366 
plotStart = 1
plot(weather$VPdef..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'VPDef')
points((exp(y_gen2_adj)-1)[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 19, cex = 1.5)
plotLim = 2000 
plotStart = 1
plot(weather$VPdef..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'VPDef')
points((exp(y_gen2_adj)-1)[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 19, cex = 1.5)
plotts.wge((exp(y_gen2_adj)-1)[(plotStart):(plotStart+plotLim-1)], main = 'VPDef')

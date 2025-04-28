# Plot VPMax for full year
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
t = 1:(6*24*366)
t = t/(6*24)
len = length(weather$date)
plot(x=t[1:len],y=weather$VPmax..mbar.[1:len], type = "p", col = "blue", pch = 16, cex = 0.5,
     xlab = "Days since 1/1/2020", ylab = "Vapor Pressure (mbar)", 
     main = "Maximum Vapor Pressure - Full Year")

# Plot a week at the beginning of the year
startDay = 0
startNum = startDay*6*24
len = 6*24*7
plot(x=t[startNum:(startNum+len)],y=weather$VPmax..mbar.[startNum:(startNum+len)], 
     type = "p", col = "blue", pch = 16, cex = 0.8,
     xlab = "Days since 1/1/2020", ylab = "Vapor Pressure (mbar)", 
     main = "Maximum Vapor Pressure - Week in Winter")

# Plot a week in the middle of the year
startDay = 180
startNum = startDay*6*24
len = 6*24*7
plot(x=t[startNum:(startNum+len)],y=weather$VPmax..mbar.[startNum:(startNum+len)], 
     type = "p", col = "blue", pch = 16, cex = 0.8,
     xlab = "Days since 1/1/2020", ylab = "Vapor Pressure (mbar)", 
     main = "Maximum Vapor Pressure - Week in Summer")

# Plot daily mean
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/vectors_VPMax.csv"
vecsVPMax = read.csv(file_path, header = TRUE)
t=1:366
plot(x=t,y=vecsVPMax$mean, type = "p", col = "blue", pch = 16, cex = 1,
     xlab = "Days since 1/1/2020", ylab = "Vapor Pressure (mbar)", 
     main = "Daily Mean of Maximum Vapor Pressure")

# Plot daily amplitude
AVec = abs(vecsVPMax$A1)
plot(x=t,y=AVec, type = "p", col = "blue", pch = 16, cex = 1,
     xlab = "Days since 1/1/2020", ylab = "Vapor Pressure (mbar)", 
     main = "Daily Amplitude of Maximum Vapor Pressure")

# Plot daily log amplitude
plot(x=t,y=log(AVec), type = "p", col = "blue", pch = 16, cex = 1,
     xlab = "Days since 1/1/2020", ylab = "Log Vapor Pressure (mbar)", 
     main = "Daily Log Amplitude of Maximum Vapor Pressure")

# Plot daily variance
plot(x=t,y=vecsVPMax$var, type = "p", col = "blue", pch = 16, cex = 1,
     xlab = "Days since 1/1/2020", ylab = "Variance", 
     main = "Daily Variance of Maximum Vapor Pressure")

# Plot daily log variance
plot(x=t,y=log(vecsVPMax$var), type = "p", col = "blue", pch = 16, cex = 1,
     xlab = "Days since 1/1/2020", ylab = "Log Variance", 
     main = "Daily Log Variance of Maximum Vapor Pressure")

# Sim mean
library(tswge)
set.seed(1)
meanVec = vecsVPMax$mean
f=1/366
ampShift = 1.5
b0A = abs(-7.9308) + ampShift
b0psi = -0.3872 + pi
b0phi = 0.7
b0var = 6
meanShift = 3
b0sin_true = b0A * cos(2 * pi * f * t + b0psi) + mean(meanVec) + meanShift
b0sin = b0sin_true + gen.arma.wge(n=366,phi=b0phi,theta = 0, vara = b0var, sn = 0, plot=FALSE) # use gen.arma since errors had relationship to prior values
b0sin[b0sin < (min(meanVec)*0.9)] = min(meanVec)*0.9
plot(x=t,y=vecsVPMax$mean, type = "p", col = "blue", pch = 16, cex = 1,
     xlab = "Days since 1/1/2020", ylab = "Vapor Pressure (mbar)", 
     main = "Single Realization of Simulated Daily Vapor Pressure Mean")
points(x=t,y=b0sin_true,col = "red", pch = 16, cex = 1)
points(x=t,y=b0sin,col = "green", pch = 16, cex = 1)
legend('topleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("blue", "red", "green"),pch=16)

# Sim amp
set.seed(1)
AVec = abs(vecsVPMax$A1)
lAvec = log(AVec)
Aamp = abs(-1.00422)
psiAmp = -0.07425 + pi
ampPhi = 0.4
ampVar = 0.4
Asin_true = exp(Aamp * cos(2 * pi * f * t + psiAmp) + mean(lAvec)) # use gen.arma since errors had relationship to prior values
Asin_error = exp(gen.arma.wge(n=366,phi=ampPhi,theta = 0, vara = ampVar, sn = 0, plot=FALSE))
ASin = Asin_true * Asin_error
ymin = min(ASin,AVec)
ymax = max(ASin,AVec)
plot(x=t,y=AVec, type = "p", col = "blue", pch = 16, cex = 1, ylim=c(ymin,ymax),
     xlab = "Days since 1/1/2020", ylab = "Vapor Pressure (mbar)", 
     main = "Single Realization of Simulated Daily Vapor Pressure Amplitude")
points(x=t,y=Asin_true,col = "red", pch = 16, cex = 1)
points(x=t,y=ASin,col = "green", pch = 16, cex = 1)
legend('topleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("blue", "red", "green"),pch=16)

# Sim var
set.seed(1)
varAmp = abs(-1.63190)
varPsi = -0.12543 + pi
varB = -3.83464
varPhi = 0.4
varVar = 0.7
varSin_true = exp(varAmp * cos(2 * pi * f * t + varPsi) + varB) # use gen.arma since errors had relationship to prior values
varSin_error = exp(gen.arma.wge(n=366,phi=varPhi,theta = 0, vara = varVar, sn = 0, plot=FALSE))
varSin = varSin_true * varSin_error
plot(x=t,y=vecsVPMax$var, type = "p", col = "blue", pch = 16, cex = 1,
     xlab = "Days since 1/1/2020", ylab = "Variance", 
     main = "Single Realization of Simulated Daily Vapor Pressure Variance")
points(x=t,y=varSin_true,col = "red", pch = 16, cex = 1)
points(x=t,y=varSin,col = "green", pch = 16, cex = 1)
legend('topleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("blue", "red", "green"),pch=16)

# Example of realization for 10 days
set.seed(1)
meanVec = vecsVPMax$mean
varVec = vecsVPMax$var
psi1vec = vecsVPMax$psi1
psi1vec = psi1vec %% (2*pi)
psi1vec2 = numeric(366)
for (day in 1:366){
  if (psi1vec[day] < pi){
    psi1vec2[day] = psi1vec[day]
  } else { # psi1vec[day] > pi
    psi1vec2[day] = psi1vec[day] - pi
  }
}
psi1mean = mean(psi1vec2) 
A1vec = vecsVPMax$A1
AVec = abs(vecsVPMax$A1)
lAvec = log(AVec)
pvec = vecsVPMax$p
phicoeff1vec = vecsVPMax$phi1
phicoeff2vec = vecsVPMax$phi2
phi1mean = mean(phicoeff1vec[pvec==2]) 
phi2mean = mean(phicoeff2vec[pvec==2])
phi = c(phi1mean,phi2mean)
ampShift = 1.5
b0A = abs(-7.9308) + ampShift
b0psi = -0.3872 + pi
b0phi = 0.7
b0var = 6
f=1/366
meanShift = 3
Aamp = abs(-1.00422)
psiAmp = -0.07425 + pi
ampPhi = 0.4
ampVar = 0.4
varAmp = abs(-1.63190)
varPsi = -0.12543 + pi
varB = -3.83464
varPhi = 0.4
varVar = 0.7
freq1 = 1/144
dayStart = 1
dayStop = 10
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  t = ((day-1)*6*24+1):(day*6*24+1)
  A1 = ASin[day]
  b0 = b0sin[day]
  psi1 = psi1mean
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  zt = gen.arma.wge(n = (6*24), phi = phi, theta = 0, vara = varSin[day], sn = 0, plot=FALSE)
  y_gen_vpmax[t] = st + zt
}
numDays = dayStop - dayStart + 1
ymin = min(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "blue",ylim=c(ymin,ymax),
     main='Realization of Maximum Vapor Pressure - No Transition', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)  # get default x-axis tick positions
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_vpmax[1:(numDays*6*24)],col = "red")
legend('topleft',legend = c("Real Data", "Simulated Data"), 
       col = c("blue", "red"),pch=1)

# Adding overlap
set.seed(1)
meanVec = vecsVPMax$mean
varVec = vecsVPMax$var
psi1vec = vecsVPMax$psi1
psi1vec = psi1vec %% (2*pi)
psi1vec2 = numeric(366)
for (day in 1:366){
  if (psi1vec[day] < pi){
    psi1vec2[day] = psi1vec[day]
  } else { # psi1vec[day] > pi
    psi1vec2[day] = psi1vec[day] - pi
  }
}
psi1mean = mean(psi1vec2) 
A1vec = vecsVPMax$A1
AVec = abs(vecsVPMax$A1)
lAvec = log(AVec)
pvec = vecsVPMax$p
phicoeff1vec = vecsVPMax$phi1
phicoeff2vec = vecsVPMax$phi2
phi1mean = mean(phicoeff1vec[pvec==2]) 
phi2mean = mean(phicoeff2vec[pvec==2])
phi = c(phi1mean,phi2mean)
ampShift = 1.5
b0A = abs(-7.9308) + ampShift
b0psi = -0.3872 + pi
b0phi = 0.7
b0var = 6
f=1/366
meanShift = 3
Aamp = abs(-1.00422)
psiAmp = -0.07425 + pi
ampPhi = 0.4
ampVar = 0.4
varAmp = abs(-1.63190)
varPsi = -0.12543 + pi
varB = -3.83464
varPhi = 0.4
varVar = 0.7
freq1 = 1/144
numInDay = 6*24
overlap = 20
dayStart = 1
dayStop = 10
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  t = ((day-1)*6*24+1):(day*6*24+1)
  A1 = ASin[day]
  b0 = b0sin[day]
  psi1 = psi1mean
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax[t] = st
}
lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
}
zt_int = data.frame(matrix(0, nrow = (6*24+1), ncol = numDays))  # add an extra value, because we want to see what the difference is
colnames(zt_int) <- paste0("day", 1:numDays)
zt_d = numeric(numDays*6*24)
for (day in 1:numDays){
  zt_int[, day] = gen.arma.wge(n = (6*24*1+1), phi = phi, theta = 0, vara = varSin[day], sn = 0, plot=FALSE)
}
zt_d[1:(6*24)] = zt_int[1:(6*24), 1] # Set the first day
for (day in 2:numDays){
  t = ((day-1)*6*24+1):(day*6*24)
  diff = zt_int[(6*24*1+1),(day-1)] - zt_int[(6*24*1),(day-1)]
  sim = zt_int[1:(6*24), day]
  zt_d[t] = sim + diff
}
y_gen_vpmax_noise = y_gen_vpmax[1:(6*24*366)] + zt_d
y_gen_vpmax_noise[y_gen_vpmax_noise < (min(weather$VPmax..mbar.)*0.5)] = min(weather$VPmax..mbar.)*0.5
numDays = dayStop - dayStart + 1
ymin = min(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_noise[1:(numDays*6*24)])
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_noise[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "blue",ylim=c(ymin,ymax),
     main='Realization of Maximum Vapor Pressure - With Transition', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)  # get default x-axis tick positions
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_vpmax_noise[1:(numDays*6*24)],col = "red")
legend('topleft',legend = c("Real Data", "Simulated Data"), 
       col = c("blue", "red"),pch=1)

# Full year
set.seed(1)
meanVec = vecsVPMax$mean
varVec = vecsVPMax$var
psi1vec = vecsVPMax$psi1
psi1vec = psi1vec %% (2*pi)
psi1vec2 = numeric(366)
for (day in 1:366){
  if (psi1vec[day] < pi){
    psi1vec2[day] = psi1vec[day]
  } else { # psi1vec[day] > pi
    psi1vec2[day] = psi1vec[day] - pi
  }
}
psi1mean = mean(psi1vec2) 
A1vec = vecsVPMax$A1
AVec = abs(vecsVPMax$A1)
lAvec = log(AVec)
pvec = vecsVPMax$p
phicoeff1vec = vecsVPMax$phi1
phicoeff2vec = vecsVPMax$phi2
phi1mean = mean(phicoeff1vec[pvec==2]) 
phi2mean = mean(phicoeff2vec[pvec==2])
phi = c(phi1mean,phi2mean)
ampShift = 1.5
b0A = abs(-7.9308) + ampShift
b0psi = -0.3872 + pi
b0phi = 0.7
b0var = 6
f=1/366
meanShift = 3
Aamp = abs(-1.00422)
psiAmp = -0.07425 + pi
ampPhi = 0.4
ampVar = 0.4
varAmp = abs(-1.63190)
varPsi = -0.12543 + pi
varB = -3.83464
varPhi = 0.4
varVar = 0.7
freq1 = 1/144
numInDay = 6*24
overlap = 20
dayStart = 1
dayStop = 366
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  t = ((day-1)*6*24+1):(day*6*24+1)
  A1 = ASin[day]
  b0 = b0sin[day]
  psi1 = psi1mean
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax[t] = st
}
lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
}
zt_int = data.frame(matrix(0, nrow = (6*24+1), ncol = numDays))  # add an extra value, because we want to see what the difference is
colnames(zt_int) <- paste0("day", 1:numDays)
zt_d = numeric(numDays*6*24)
for (day in 1:numDays){
  zt_int[, day] = gen.arma.wge(n = (6*24*1+1), phi = phi, theta = 0, vara = varSin[day], sn = 0, plot=FALSE)
}
zt_d[1:(6*24)] = zt_int[1:(6*24), 1] # Set the first day
for (day in 2:numDays){
  t = ((day-1)*6*24+1):(day*6*24)
  diff = zt_int[(6*24*1+1),(day-1)] - zt_int[(6*24*1),(day-1)]
  sim = zt_int[1:(6*24), day]
  zt_d[t] = sim + diff
}
y_gen_vpmax_noise = y_gen_vpmax[1:(6*24*366)] + zt_d
y_gen_vpmax_noise[y_gen_vpmax_noise < (min(weather$VPmax..mbar.)*0.5)] = min(weather$VPmax..mbar.)*0.5
numDays = dayStop - dayStart + 1
ymin = min(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_noise[1:(numDays*6*24)])
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_noise[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "blue",ylim=c(ymin,ymax),
     main='Realization of Maximum Vapor Pressure - Full Year', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)  # get default x-axis tick positions
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_vpmax_noise[1:(numDays*6*24)],col = "red")
legend('topleft',legend = c("Real Data", "Simulated Data"), 
       col = c("blue", "red"),pch=1)

# Plot of TDew for full year
t = 1:(6*24*366)
t = t/(6*24)
len = length(weather$date)
plot(x=t[1:len],y=weather$Tdew..degC.[1:len], type = "p", col = "blue", pch = 16, cex = 0.5,
     xlab = "Days since 1/1/2020", ylab = "Dew Temperature (C)", 
     main = "Dew Temperature - Full Year")

# Plot a week of TDew at the beginning of the year
startDay = 0
startNum = startDay*6*24
len = 6*24*7
plot(x=t[startNum:(startNum+len)],y=weather$Tdew..degC.[startNum:(startNum+len)], 
     type = "p", col = "blue", pch = 16, cex = 0.8,
     xlab = "Days since 1/1/2020", ylab = "Dew Temperature (C)", 
     main = "Dew Temperature - Week in Winter")

# Plot a week of TDew in the middle of the year
startDay = 180
startNum = startDay*6*24
len = 6*24*7
plot(x=t[startNum:(startNum+len)],y=weather$Tdew..degC.[startNum:(startNum+len)], 
     type = "p", col = "blue", pch = 16, cex = 0.8,
     xlab = "Days since 1/1/2020", ylab = "Dew Temperature (C)", 
     main = "Dew Temperature - Week in Summer")

# Plot daily amplitude for TDew
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/vectors_TDew.csv"
vecs_TDew = read.csv(file_path, header = TRUE)
AVec = abs(vecs_TDew$A1)
t=1:366
plot(x=t,y=AVec, type = "p", col = "blue", pch = 16, cex = 1,
     xlab = "Days since 1/1/2020", ylab = "Dew Temperature (C)", 
     main = "Daily Amplitude of Dew Temperature")

# Plot 1 week of a TDew realization
set.seed(1)
psi1vec = vecs_TDew$psi1
psi1vec = psi1vec %% (2*pi)
psi1vec2 = numeric(366)
for (day in 1:366){
  if (psi1vec[day] < pi){
    psi1vec2[day] = psi1vec[day]
  } else { # psi1vec[day] > pi
    psi1vec2[day] = psi1vec[day] - pi
  }
}
psi1mean = mean(psi1vec2) 
phi = 0.9
f=1/366
b0A = -6.55756
b0Psi = 5.57984
b0B = 5.40932
# b0_phi = 0.7
b0_phi = 0.9
b0_avar = 5.48
AA = -0.29247
APsi = 1.15203
AB = 1.54916
# Amp_phi = 0.1603404
Amp_phi = 0.7
Amp_avar = 0.9459506
varA = -0.89346    
varPsi = 0.18136
varB = -3.75244
# Var_phi = 0.2812496
Var_phi = 0.7
Var_avar = 0.4903821
freq1 = 1/144
psi1 = psi1mean
dayStart = 1
dayStop = 366
numInDay = 6*24
overlap = 20
t=1:366
b0sin_true = b0A * cos(2 * pi * f * t + b0Psi) + b0B
b0sin = b0sin_true + gen.arma.wge(n=366, phi=b0_phi, theta=0, sn=0, vara =b0_avar, plot = FALSE)
Asin_true = AA * cos(2 * pi * f * t + APsi) + AB
Asin = Asin_true + gen.arma.wge(n=366, phi=Amp_phi, theta=0, sn=0, vara =Amp_avar, plot = FALSE)
logVar_sin_true = varA * cos(2 * pi * f * t + varPsi) + varB
logVar_sin = logVar_sin_true + gen.arma.wge(n=366, phi=Var_phi, theta=0, sn=0, vara =Var_avar, plot = FALSE)
var_sin = exp(logVar_sin)
var_sin_true = exp(logVar_sin_true)
y_gen_tdew_true = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop) {
  t_day = ((day-1)*6*24+1):(day*6*24)
  A1 = Asin_true[day]
  b0 = b0sin_true[day]
  st = b0 + A1 * cos(2 * pi * freq1 * (1:length(t_day)) + psi1)
  y_gen_tdew_true[t_day] = st
}
lastOfDay = seq(numInDay, length(y_gen_tdew_true)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_tdew_true)-numInDay, by = numInDay)
for (ii in 1:overlap) {
  beforeDay = seq(numInDay+1-ii, length(y_gen_tdew_true)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_tdew_true)-numInDay+overlap, by = numInDay)
  factor = overlap + 1 - ii
  y_gen_tdew_true[beforeDay] = ((y_gen_tdew_true[lastOfDay] + y_gen_tdew_true[firstOfDay]) * factor/2 + y_gen_tdew_true[beforeDay] * ii) / (overlap + 1)
  y_gen_tdew_true[afterDay] = ((y_gen_tdew_true[lastOfDay] + y_gen_tdew_true[firstOfDay]) * factor/2 + y_gen_tdew_true[afterDay] * ii) / (overlap + 1)
}
y_gen_tdew = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop) {
  t_day = ((day-1)*6*24+1):(day*6*24)
  A1 = Asin[day]
  b0 = b0sin[day]
  st = b0 + A1 * cos(2 * pi * freq1 * (1:length(t_day)) + psi1)
  y_gen_tdew[t_day] = st
}
lastOfDay = seq(numInDay, length(y_gen_tdew)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_tdew)-numInDay, by = numInDay)
for (ii in 1:overlap) {
  beforeDay = seq(numInDay+1-ii, length(y_gen_tdew)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_tdew)-numInDay+overlap, by = numInDay)
  factor = overlap + 1 - ii
  y_gen_tdew[beforeDay] = ((y_gen_tdew[lastOfDay] + y_gen_tdew[firstOfDay]) * factor/2 + y_gen_tdew[beforeDay] * ii) / (overlap + 1)
  y_gen_tdew[afterDay] = ((y_gen_tdew[lastOfDay] + y_gen_tdew[firstOfDay]) * factor/2 + y_gen_tdew[afterDay] * ii) / (overlap + 1)
}
zt_int = data.frame(matrix(0, nrow = (6*24+1), ncol = numDays))  # plus one because of connecting days
colnames(zt_int) <- paste0("day", 1:numDays)
zt_d = numeric(numDays*6*24)
for (day in 1:numDays) {
  zt_int[, day] = gen.arma.wge(n = (6*24+1), phi = phi, theta = 0, vara = var_sin[day], sn = 0, plot=FALSE)
}
zt_d[1:(6*24)] = zt_int[1:(6*24), 1]
for (day in 2:numDays) {
  t_day = ((day-1)*6*24+1):(day*6*24)
  diff = zt_int[(6*24+1), (day-1)] - zt_int[(6*24), (day-1)]
  sim = zt_int[1:(6*24), day]
  zt_d[t_day] = sim + diff
}
y_gen_tdew_noise = y_gen_tdew[1:(6*24*366)] + zt_d
numDays = 7
ymin = min(weather$Tdew..degC.[1:(numDays*6*24)],y_gen_tdew_noise[1:(numDays*6*24)])
ymax = max(weather$Tdew..degC.[1:(numDays*6*24)],y_gen_tdew_noise[1:(numDays*6*24)])
plot(weather$Tdew..degC.[1:(numDays*6*24)],col = "blue",ylim=c(ymin,ymax),
     main='Realization of Dew Temperature - First Week', xlab = 'Days',ylab = 'Dew Temperature (C)',xaxt='n')
ticks <- axTicks(1)  # get default x-axis tick positions
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_tdew_noise[1:(numDays*6*24)],col = "red")
legend('topleft',legend = c("Real Data", "Simulated Data"), 
       col = c("blue", "red"),pch=1)

# Plot the full year
numDays = 366
ymin = min(weather$Tdew..degC.[1:(numDays*6*24)],y_gen_tdew_noise[1:(numDays*6*24)])
ymax = max(weather$Tdew..degC.[1:(numDays*6*24)],y_gen_tdew_noise[1:(numDays*6*24)])
plot(weather$Tdew..degC.[1:(numDays*6*24)],col = "blue",ylim=c(ymin,ymax),
     main='Realization of Dew Temperature - Full Year', xlab = 'Days',ylab = 'Dew Temperature (C)',xaxt='n')
ticks <- axTicks(1)  # get default x-axis tick positions
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_tdew_noise[1:(numDays*6*24)],col = "red")
legend('topleft',legend = c("Real Data", "Simulated Data"), 
       col = c("blue", "red"),pch=1)

# ACF of VPMax
psi1mean = mean(psi1vec2) 
A1vec = vecsVPMax$A1
AVec = abs(vecsVPMax$A1)
lAvec = log(AVec)
pvec = vecsVPMax$p
phicoeff1vec = vecsVPMax$phi1
phicoeff2vec = vecsVPMax$phi2
phi1mean = mean(phicoeff1vec[pvec==2]) 
phi2mean = mean(phicoeff2vec[pvec==2])
phi = c(phi1mean,phi2mean)
# phiShift = 0
phiShift = 0.1
ampShift = 1.5
b0A = abs(-7.9308) + ampShift
b0psi = -0.3872 + pi
b0phi = 0.7
if (phiShift > 0){
  b0phi = b0phi + (1-b0phi)*phiShift
} else {
  b0phi = b0phi - b0phi*phiShift
}
b0var = 6
f=1/366
meanShift = 3
Aamp = abs(-1.00422)
psiAmp = -0.07425 + pi
ampPhi = 0.4
if (phiShift > 0){
  ampPhi = ampPhi + (1-ampPhi)*phiShift
} else {
  ampPhi = ampPhi - ampPhi*phiShift
}
ampVar = 0.4
varAmp = abs(-1.63190)
varPsi = -0.12543 + pi
varB = -3.83464
varPhi = 0.4
if (phiShift > 0){
  varPhi = varPhi + (1-varPhi)*phiShift
} else {
  varPhi = varPhi - varPhi*phiShift
}
varVar = 0.7
minVPMax = 0.6
t=1:366
numDays = 366
b0sin_true = b0A * cos(2 * pi * f * t + b0psi) + mean(meanVec) + meanShift
Asin_true = exp(Aamp * cos(2 * pi * f * t + psiAmp) + mean(lAvec))
y_gen_vpmax_true = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  t = ((day-1)*6*24+1):(day*6*24+1)
  A1 = Asin_true[day]
  b0 = b0sin_true[day]
  psi1 = psi1mean
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax_true[t] = st
}
lastOfDay = seq(numInDay, length(y_gen_vpmax_true)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax_true)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax_true)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax_true)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax_true[beforeDay] = ((y_gen_vpmax_true[lastOfDay]+y_gen_vpmax_true[firstOfDay])*factor/2 + y_gen_vpmax_true[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax_true[afterDay] = ((y_gen_vpmax_true[lastOfDay]+y_gen_vpmax_true[firstOfDay])*factor/2 + y_gen_vpmax_true[afterDay]*ii)/(overlap+1)
}
minVPMax = 0.6
# numSims = 20
numSims = 100
y_gen_vpmax_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_mult) <- paste0("sim", 1:numSims)
for (seed in 1:numSims){
  if (seed %% 10 == 0){
    print(seed)
  }
  set.seed(seed)
  t = 1:366
  b0sin = b0A * cos(2 * pi * f * t + b0psi) + mean(meanVec) + meanShift + gen.arma.wge(n=366,phi=b0phi,theta = 0, vara = b0var, sn = 0, plot=FALSE)
  b0sin[b0sin < (min(meanVec)*0.9)] = min(meanVec)*0.9
  ASin = exp(Aamp * cos(2 * pi * f * t + psiAmp) + mean(lAvec) + gen.arma.wge(n=366,phi=ampPhi,theta = 0, vara = ampVar, sn = 0, plot=FALSE))
  varSin = exp(varAmp * cos(2 * pi * f * t + varPsi) + varB + gen.arma.wge(n=366,phi=varPhi,theta = 0, vara = varVar, sn = 0, plot=FALSE))
  y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
  for (day in dayStart:dayStop){
    t = ((day-1)*6*24+1):(day*6*24)
    A1 = ASin[day] 
    b0 = b0sin[day] 
    psi1 = psi1mean
    st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
    y_gen_vpmax[t] = st
  }
  lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
  firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
  for (ii in 1:overlap){
    beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
    afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
    factor = overlap+1-ii
    y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
    y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
  }
  numDays = 366
  zt_int = data.frame(matrix(0, nrow = (6*24+1), ncol = numDays))  
  colnames(zt_int) <- paste0("day", 1:numDays)
  zt_d = numeric(numDays*6*24)
  for (day in 1:numDays){
    zt_int[, day] = gen.arma.wge(n = (6*24*1+1), phi = phi, theta = 0, vara = varSin[day], sn = 0, plot=FALSE)
  }
  zt_d[1:(6*24)] = zt_int[1:(6*24), 1] 
  for (day in 2:numDays){
    t = ((day-1)*6*24+1):(day*6*24)
    diff = zt_int[(6*24*1+1),(day-1)] - zt_int[(6*24*1),(day-1)]
    sim = zt_int[1:(6*24), day]
    zt_d[t] = sim + diff
  }
  y_gen_vpmax_noise = y_gen_vpmax[1:(6*24*366)] + zt_d
  y_gen_vpmax_noise[y_gen_vpmax_noise < (min(weather$VPmax..mbar.)*minVPMax)] = min(weather$VPmax..mbar.)*minVPMax
  y_gen_vpmax_mult[, seed] = y_gen_vpmax_noise
}
numSims = 20
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6,
     main='ACF Comparison - Maximum Vapor Pressure', xlab = 'Lag',ylab = 'ACF')
for (sim in 1:numSims){
  ACF2 = acf(y_gen_vpmax_mult[,sim], plot = "FALSE")
  lines(ACF2$lag, ACF2$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1)

# ACF of TDew
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/vectors_TDew.csv"
vecs_TDew = read.csv(file_path, header = TRUE)
meanVec = vecs_TDew$mean
varVec = vecs_TDew$var
psi1vec = vecs_TDew$psi1
psi1vec = psi1vec %% (2*pi)
psi1vec2 = numeric(366)
for (day in 1:366){
  if (psi1vec[day] < pi){
    psi1vec2[day] = psi1vec[day]
  } else { # psi1vec[day] > pi
    psi1vec2[day] = psi1vec[day] - pi
  }
}
psi1mean = mean(psi1vec2) 
A1vec = vecs_TDew$A1
AVec = abs(vecs_TDew$A1)
lAvec = log(AVec)
pvec = vecs_TDew$p
phicoeff1vec = vecs_TDew$phi1
phicoeff2vec = vecs_TDew$phi2
phi1mean = mean(phicoeff1vec[pvec==2]) 
phi2mean = mean(phicoeff2vec[pvec==2])
# phi = c(phi1mean,phi2mean)
# phiShift = 0
phiShift = 0.4
phi = 0.7
# phi = 0.9
if (phiShift > 0){
  phi = phi + (1-phi)*phiShift
} else {
  phi = phi - phi*phiShift
}
f=1/366
b0A = -6.55756
b0Psi = 5.57984
b0B = 5.40932
b0_phi = 0.7
# b0_phi = 0.9
if (phiShift > 0){
  b0_phi = b0_phi + (1-b0_phi)*phiShift
} else {
  b0_phi = b0_phi - b0_phi*phiShift
}
b0_avar = 5.48
AA = -0.29247
APsi = 1.15203
AB = 1.54916
Amp_phi = 0.1603404
# Amp_phi = 0.7
if (phiShift > 0){
  Amp_phi = Amp_phi + (1-Amp_phi)*phiShift
} else {
  Amp_phi = Amp_phi - Amp_phi*phiShift
}
Amp_avar = 0.9459506
varA = -0.89346    
varPsi = 0.18136
varB = -3.75244
Var_phi = 0.2812496
# Var_phi = 0.7
if (phiShift > 0){
  Var_phi = Var_phi + (1-Var_phi)*phiShift
} else {
  Var_phi = Var_phi - Var_phi*phiShift
}
Var_avar = 0.4903821
t=1:366
b0sin_true = b0A * cos(2 * pi * f * t + b0Psi) + b0B
Asin_true = AA * cos(2 * pi * f * t + APsi) + AB
y_gen_tdew_true = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop) {
  t_day = ((day-1)*6*24+1):(day*6*24)
  A1 = Asin_true[day]
  b0 = b0sin_true[day]
  st = b0 + A1 * cos(2 * pi * freq1 * (1:length(t_day)) + psi1)
  y_gen_tdew_true[t_day] = st
}
lastOfDay = seq(numInDay, length(y_gen_tdew_true)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_tdew_true)-numInDay, by = numInDay)
for (ii in 1:overlap) {
  beforeDay = seq(numInDay+1-ii, length(y_gen_tdew_true)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_tdew_true)-numInDay+overlap, by = numInDay)
  factor = overlap + 1 - ii
  y_gen_tdew_true[beforeDay] = ((y_gen_tdew_true[lastOfDay] + y_gen_tdew_true[firstOfDay]) * factor/2 + y_gen_tdew_true[beforeDay] * ii) / (overlap + 1)
  y_gen_tdew_true[afterDay] = ((y_gen_tdew_true[lastOfDay] + y_gen_tdew_true[firstOfDay]) * factor/2 + y_gen_tdew_true[afterDay] * ii) / (overlap + 1)
}
# numSims = 20
numSims = 100
y_gen_tdew_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_mult) <- paste0("sim", 1:numSims)
for (seed in 1:numSims){
  if (seed %% 10 == 0){
    print(seed)
  }
  set.seed(seed)
  t = 1:366
  b0sin_true = b0A * cos(2 * pi * f * t + b0Psi) + b0B
  b0sin = b0sin_true + gen.arma.wge(n=366, phi=0.7, theta=0, sn=0, vara =5.48, plot = FALSE)
  Asin_true = AA * cos(2 * pi * f * t + APsi) + AB
  Asin = Asin_true + gen.arma.wge(n=366, phi=Amp_phi, theta=0, sn=0, vara =Amp_avar, plot = FALSE)
  logVar_sin_true = varA * cos(2 * pi * f * t + varPsi) + varB
  logVar_sin = logVar_sin_true + gen.arma.wge(n=366, phi=Var_phi, theta=0, sn=0, vara =Var_avar, plot = FALSE)
  var_sin = exp(logVar_sin)
  var_sin_true = exp(logVar_sin_true)
  y_gen_tdew = numeric(6*24*(dayStop-dayStart+1))
  for (day in dayStart:dayStop){
    t = ((day-1)*6*24+1):(day*6*24)
    A1 = Asin[day] 
    b0 = b0sin[day] 
    st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
    y_gen_tdew[t] = st
  }
  lastOfDay = seq(numInDay, length(y_gen_tdew)-numInDay, by = numInDay)
  firstOfDay = seq(numInDay+1, length(y_gen_tdew)-numInDay, by = numInDay)
  for (ii in 1:overlap){
    beforeDay = seq(numInDay+1-ii, length(y_gen_tdew)-numInDay, by = numInDay)
    afterDay = seq(numInDay+ii, length(y_gen_tdew)-numInDay+overlap, by = numInDay)
    factor = overlap+1-ii
    y_gen_tdew[beforeDay] = ((y_gen_tdew[lastOfDay]+y_gen_tdew[firstOfDay])*factor/2 + y_gen_tdew[beforeDay]*ii)/(overlap+1)
    y_gen_tdew[afterDay] = ((y_gen_tdew[lastOfDay]+y_gen_tdew[firstOfDay])*factor/2 + y_gen_tdew[afterDay]*ii)/(overlap+1)
  }
  numDays = 366
  zt_int = data.frame(matrix(0, nrow = (6*24+1), ncol = numDays))  # plus one because of connecting days
  colnames(zt_int) <- paste0("day", 1:numDays)
  zt_d = numeric(numDays*6*24)
  for (day in 1:numDays) {
    zt_int[, day] = gen.arma.wge(n = (6*24+1), phi = phi, theta = 0, vara = var_sin[day], sn = 0, plot=FALSE)
  }
  zt_d[1:(6*24)] = zt_int[1:(6*24), 1]
  for (day in 2:numDays) {
    t_day = ((day-1)*6*24+1):(day*6*24)
    diff = zt_int[(6*24+1), (day-1)] - zt_int[(6*24), (day-1)]
    sim = zt_int[1:(6*24), day]
    zt_d[t_day] = sim + diff
  }
  y_gen_tdew_noise = y_gen_tdew[1:(6*24*366)] + zt_d
  # y_gen_tdew_noise[y_gen_tdew_noise < (min(weather$Tdew..degC., na.rm=TRUE)*0.5)] = min(weather$Tdew..degC., na.rm=TRUE)*0.5
  y_gen_tdew_mult[, seed] = y_gen_tdew_noise
}
numSims = 20
ACF = acf(weather$Tdew..degC., plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6,
     main='ACF Comparison - Dew Temperature', xlab = 'Lag',ylab = 'ACF')
for (sim in 1:numSims){
  ACF2 = acf(y_gen_tdew_mult[,sim], plot = "FALSE")
  lines(ACF2$lag, ACF2$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1)

# Spectral Density of VPMax, first 20 days
numDays = 20
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, 
     main='Spectral Density Comparison - First 20 Days - Maximum Vapor Pressure', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}
legend('topright',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1)

# Spectral Density of VPMax, first 20 days, daily frequency
trunc = 1000
lenPlot = 35
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE", trunc=trunc)
plot(SpecDen$freq[1:lenPlot],SpecDen$pzgram[1:lenPlot], type = "l", lwd = 6, 
     main='Spectral Density Comparison - First 20 Days - Daily Frequency - Maximum Vapor Pressure', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[1:(6*24*numDays),sim], plot = "FALSE", trunc=trunc)
  lines(SpecDen2$freq[1:lenPlot],SpecDen2$pzgram[1:lenPlot], lwd = 2, col = "red")
}
legend('topright',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1)

# Spectral Density of VPMax, middle 20 days
numDays = 20
startDay = 180
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, 
     main='Spectral Density Comparison - Middle 20 Days - Maximum Vapor Pressure', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}
legend('topright',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1)

# Spectral Density of TDew, first 20 days
numDays = 20
SpecDen = parzen.wge(weather$Tdew..degC.[1:(6*24*numDays)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, 
     main='Spectral Density Comparison - First 20 Days - Dew Temperature', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_tdew_mult[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}
legend('topright',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1)

# Spectral Density of TDew, first 20 days, daily frequency
lenPlot = 35
SpecDen = parzen.wge(weather$Tdew..degC.[1:(6*24*numDays)], plot = "FALSE", trunc=trunc)
plot(SpecDen$freq[1:lenPlot],SpecDen$pzgram[1:lenPlot], type = "l", lwd = 6, 
     main='Spectral Density Comparison - First 20 Days - Daily Frequency - Dew Temperature', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_tdew_mult[1:(6*24*numDays),sim], plot = "FALSE", trunc=trunc)
  lines(SpecDen2$freq[1:lenPlot],SpecDen2$pzgram[1:lenPlot], lwd = 2, col = "red")
}
legend('topright',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1)

# Calculating Correlation between VPMax and TDew
corVD = cor(weather$VPmax..mbar.,weather$Tdew..degC.)
numSims = 100
corVec = numeric(numSims)
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  diff = corVD - cor(y_gen_vpmax_mult[sim],y_gen_tdew_mult[sim])
  # alpha <- diff * (1-abs(diff))
  alpha <- diff * (1-abs(diff))^2
  y_gen_tdew_adj <- (1 - alpha) * y_gen_tdew_mult[sim] + alpha * y_gen_vpmax_mult[sim]
  y_gen_vpmax_adj <- (1 - alpha) * y_gen_vpmax_mult[sim] + alpha * y_gen_tdew_adj
  # y_gen_vpmax_adj[y_gen_vpmax_adj<0] = 0 # remove for now, but I'll need to deal with negatives for VPMax
  # corVec[sim] = cor(y_gen_vpmax_adj, y_gen_tdew_adj)
  corVec[sim] = cor(y_gen_vpmax_mult[sim], y_gen_tdew_mult[sim])
}
hist(corVec,breaks=50,xlab='Correlation',ylab='Count',
     main = 'Histogram of Simulation Correlations between\nMaximum Vapor Pressure and Dew Temperature')
abline(v=corVD,col='red')
legend('topright',legend='Actual Correlation = 0.713',col='red',pch=16)

# Realization of Act Temp using linear transformation
numSims = 1
y_gen_TLin_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_TLin_adj) <- paste0("sim", 1:numSims)
y_gen_TNonlin_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_TNonlin_adj) <- paste0("sim", 1:numSims)
varTemp = 0.1
minTemp = 0.7
for (sim in 1:numSims){
  # Oh shoot, I might need to add some error variance as well
  set.seed(sim)
  varLin = rnorm(n=length(y_gen_TLin_adj[sim]),mean=0,sd=sqrt(varTemp))
  varNonlin = rnorm(n=length(y_gen_TNonlin_adj[sim]),mean=0,sd=sqrt(varTemp))
  
  # y_gen_TLin_adj[sim] = -2.274795  + y_gen_vpmax_mult[sim]*0.815397 + y_gen_tdew_mult[sim]*0.236779 + varLin
  y_gen_TLin_adj[,sim] = -2.274795  + y_gen_vpmax_mult[,sim]*0.815397 + y_gen_tdew_mult[,sim]*0.236779 + varLin
  # y_gen_TNonlin_adj[sim] = 44.06+243.5*log(y_gen_vpmax_mult[sim]/100)/(17.62-log(y_gen_vpmax_mult[sim]/100))*1.326+y_gen_vpmax_mult[sim]*2.102e-02-y_gen_tdew_mult[sim]*1.000e-03+y_gen_vpmax_mult[sim]*y_gen_tdew_mult[sim]*6.250e-05 + varNonlin
  y_gen_TNonlin_adj[,sim] = 44.06+243.5*log(y_gen_vpmax_mult[,sim]/100)/(17.62-log(y_gen_vpmax_mult[,sim]/100))*1.326+y_gen_vpmax_mult[,sim]*2.102e-02-y_gen_tdew_mult[,sim]*1.000e-03+y_gen_vpmax_mult[,sim]*y_gen_tdew_mult[,sim]*6.250e-05 + varNonlin
  
}
t = 1:(6*24*366)
t = t/(6*24)
numDays = 366
len = length(weather$date)
plot(x=t[1:len],y=weather$T..degC.[1:len], type = "p", col = "blue", pch = 16, cex = 0.5,
     xlab = "Days since 1/1/2020", ylab = "Air Temperature (C)", 
     main = "Air Temperature Simulations")
points(x=t[1:len],y=y_gen_TLin_adj[1:len,1],col = "red",pch = 16, cex = 0.5)
points(x=t[1:len],y=y_gen_TNonlin_adj[1:len,1],col = "green",pch = 16, cex = 0.5)
legend('topleft',legend = c("Real Data", "Simulated Data - Linear", "Simulated Data - Non Linear"), 
       col = c("blue", "red", "green"),pch=1)
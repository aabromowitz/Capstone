# Demo

# Zoomed in spectral density
library(tswge)
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/vectors_VPMax.csv"
vecs = read.csv(file_path, header = TRUE)
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
meanVec = vecs$mean
varVec = vecs$var
psi1vec = vecs$psi1
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
A1vec = vecs$A1
AVec = abs(vecs$A1)
lAvec = log(AVec)
pvec = vecs$p
phicoeff1vec = vecs$phi1
phicoeff2vec = vecs$phi2
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
dayStop = 366
numInDay = 6*24
overlap = 20
set.seed(1)
t=1:366
numDays = 366
b0sin_true = b0A * cos(2 * pi * f * t + b0psi) + mean(meanVec) + meanShift
b0sin = b0sin_true + gen.arma.wge(n=366,phi=b0phi,theta = 0, vara = b0var, sn = 0, plot=FALSE)
b0sin[b0sin < (min(meanVec)*0.9)] = min(meanVec)*0.9
Asin_true = exp(Aamp * cos(2 * pi * f * t + psiAmp) + mean(lAvec))
Asin_error = exp(gen.arma.wge(n=366,phi=ampPhi,theta = 0, vara = ampVar, sn = 0, plot=FALSE))
ASin = Asin_true * Asin_error
varSin_true = exp(varAmp * cos(2 * pi * f * t + varPsi) + varB)
varSin_error = exp(gen.arma.wge(n=366,phi=varPhi,theta = 0, vara = varVar, sn = 0, plot=FALSE))
varSin = varSin_true * varSin_error
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
numSims = 20
y_gen_vpmax_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_mult) <- paste0("sim", 1:numSims)
for (seed in 1:numSims){
  print(seed)
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
  y_gen_vpmax_noise[y_gen_vpmax_noise < (min(weather$VPmax..mbar.)*0.5)] = min(weather$VPmax..mbar.)*0.5
  y_gen_vpmax_mult[, seed] = y_gen_vpmax_noise
}
numDays = 20
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE")
SpecDen2 = parzen.wge(y_gen_vpmax_true[1:(6*24*numDays)], plot = "FALSE")
ymax = max(SpecDen$pzgram,SpecDen2$pzgram)
ymin = min(SpecDen$pzgram,SpecDen2$pzgram)
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - First 20 Days', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_mult[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen3$freq,SpecDen3$pzgram, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
startDay = 180
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
SpecDen2 = parzen.wge(y_gen_vpmax_true[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
ymax = max(SpecDen$pzgram,SpecDen2$pzgram)
ymin = min(SpecDen$pzgram,SpecDen2$pzgram)
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - Middle 20 Days', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_mult[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE")
  lines(SpecDen3$freq,SpecDen3$pzgram, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
trunc = 1000
lenPlot = 35
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE", trunc=trunc)
SpecDen2 = parzen.wge(y_gen_vpmax_true[1:(6*24*numDays)], plot = "FALSE", trunc=trunc)
ymax = max(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
ymin = min(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
plot(SpecDen$freq[1:lenPlot],SpecDen$pzgram[1:lenPlot], type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - First 40 Days - First Peak', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq[1:lenPlot],SpecDen2$pzgram[1:lenPlot], lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_mult[1:(6*24*numDays),sim], plot = "FALSE", trunc=trunc)
  lines(SpecDen3$freq[1:lenPlot],SpecDen3$pzgram[1:lenPlot], lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE", trunc=trunc)
SpecDen2 = parzen.wge(y_gen_vpmax_true[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE", trunc=trunc)
ymax = max(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
ymin = min(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
plot(SpecDen$freq[1:lenPlot],SpecDen$pzgram[1:lenPlot], type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - First 40 Days - First Peak', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq[1:lenPlot],SpecDen2$pzgram[1:lenPlot], lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_mult[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE", trunc=trunc)
  lines(SpecDen3$freq[1:lenPlot],SpecDen3$pzgram[1:lenPlot], lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
lenPlot = 50
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE", trunc=trunc)
SpecDen2 = parzen.wge(y_gen_vpmax_true[1:(6*24*numDays)], plot = "FALSE", trunc=trunc)
ymax = max(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
ymin = min(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
plot(SpecDen$freq[1:lenPlot],SpecDen$pzgram[1:lenPlot], type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - First 40 Days - First Two Peaks', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq[1:lenPlot],SpecDen2$pzgram[1:lenPlot], lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_mult[1:(6*24*numDays),sim], plot = "FALSE", trunc=trunc)
  lines(SpecDen3$freq[1:lenPlot],SpecDen3$pzgram[1:lenPlot], lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE", trunc=trunc)
SpecDen2 = parzen.wge(y_gen_vpmax_true[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE", trunc=trunc)
ymax = max(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
ymin = min(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
plot(SpecDen$freq[1:lenPlot],SpecDen$pzgram[1:lenPlot], type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - First 40 Days - First Two Peaks', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq[1:lenPlot],SpecDen2$pzgram[1:lenPlot], lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_mult[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE", trunc=trunc)
  lines(SpecDen3$freq[1:lenPlot],SpecDen3$pzgram[1:lenPlot], lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)

# CCFs
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
  y_gen_vpmax_noise[y_gen_vpmax_noise < (min(weather$VPmax..mbar.)*0.5)] = min(weather$VPmax..mbar.)*0.5
  y_gen_vpmax_mult[, seed] = y_gen_vpmax_noise
}
y_gen_tdew_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_mult) <- paste0("sim", 1:numSims)
for (seed in 1:numSims){
  y_gen_tdew_mult[, seed] = weather$Tdew..degC.
}
corVD = cor(weather$VPmax..mbar.,weather$Tdew..degC.)
y_gen_vpmax_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_adj) <- paste0("sim", 1:numSims)
y_gen_tdew_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_adj) <- paste0("sim", 1:numSims)
corVec = numeric(numSims)
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  y_gen_vpmax_adj[sim] = y_gen_vpmax_mult[sim]
  y_gen_tdew_adj[sim] = y_gen_tdew_mult[sim]
  diff = corVD - cor(y_gen_vpmax_mult[sim],y_gen_tdew_mult[sim])
  alpha <- diff * (1-abs(diff))^2
  y_gen_tdew_adj[,sim] <- (1 - alpha) * y_gen_tdew_mult[sim] + alpha * y_gen_vpmax_adj[sim]
  y_gen_vpmax_adj[,sim] <- (1 - alpha) * y_gen_vpmax_adj[sim] + alpha * y_gen_tdew_adj[,sim]
  corVec[sim] = cor(y_gen_vpmax_adj[,sim], y_gen_tdew_adj[,sim])
}
hist(corVec, breaks=20,main='Histogram of Correlations\nbetween VPMax (Generated) and TDew',ylab='Correlation')
abline(v=corVD, col='red')
numDays = 366
ymin = min(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_mult[1:(numDays*6*24),1])
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_mult[1:(numDays*6*24),1])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of VPMax - 366 Days', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)  # get default x-axis tick positions
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_vpmax_mult[1:(numDays*6*24),1],col = "blue")
points(y_gen_vpmax_adj[1:(numDays*6*24),1],col = "green")
legend('topleft',legend = c("Real Data", "Simulated Data - Original", "Simulated Data - After Update"), 
       col = c("red", "blue", "green"),pch=1)
numDays = 10
ymin = min(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_mult[1:(numDays*6*24),1])
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_mult[1:(numDays*6*24),1])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of VPMax - First 10 Days', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_vpmax_mult[1:(numDays*6*24),1],col = "blue")
points(y_gen_vpmax_adj[1:(numDays*6*24),1],col = "green")
legend('top',legend = c("Real Data", "Simulated Data - Original", "Simulated Data - After Update"), 
       col = c("red", "blue", "green"),pch=1)
startDay = 366/2
numDays = 10
ymin = min(c(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax_mult[(1+6*24*startDay):(6*24*startDay+numDays*6*24),1]))
ymax = max(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax_mult[(1+6*24*startDay):(6*24*startDay+numDays*6*24),1])
plot(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of VPMax - Middle 10 Days', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)
axis(1, at = ticks, labels = round(ticks / 144,0) + startDay)
points(y_gen_vpmax_mult[(1+6*24*startDay):(6*24*startDay+numDays*6*24),1],col = "blue")
points(y_gen_vpmax_adj[(1+6*24*startDay):(6*24*startDay+numDays*6*24),1],col = "green")
legend('topleft',legend = c("Real Data", "Simulated Data - Original", "Simulated Data - After Update"), 
       col = c("red", "blue", "green"),pch=1)
numSims = 20
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
ACF2 = acf(y_gen_vpmax_true, plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6,
     main='ACF Comparison - After Update', xlab = 'Lag',ylab = 'ACF')
lines(ACF2$lag, ACF2$acf, lwd = 6, col = "blue")
for (sim in 1:numSims){
  ACF3 = acf(y_gen_vpmax_adj[,sim], plot = "FALSE")
  lines(ACF3$lag, ACF3$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
numDays = 20
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE")
SpecDen2 = parzen.wge(y_gen_vpmax_true[1:(6*24*numDays)], plot = "FALSE")
ymax = max(SpecDen$pzgram,SpecDen2$pzgram)
ymin = min(SpecDen$pzgram,SpecDen2$pzgram)
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - First 20 Days - After Update', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_adj[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen3$freq,SpecDen3$pzgram, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
startDay = 180
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
SpecDen2 = parzen.wge(y_gen_vpmax_true[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
ymax = max(SpecDen$pzgram,SpecDen2$pzgram)
ymin = min(SpecDen$pzgram,SpecDen2$pzgram)
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - Middle 20 Days - After Update', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_adj[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE")
  lines(SpecDen3$freq,SpecDen3$pzgram, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
trunc = 1000
lenPlot = 35
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE", trunc=trunc)
SpecDen2 = parzen.wge(y_gen_vpmax_true[1:(6*24*numDays)], plot = "FALSE", trunc=trunc)
ymax = max(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
ymin = min(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
plot(SpecDen$freq[1:lenPlot],SpecDen$pzgram[1:lenPlot], type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - First 40 Days - First Peak - After Update', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq[1:lenPlot],SpecDen2$pzgram[1:lenPlot], lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_adj[1:(6*24*numDays),sim], plot = "FALSE", trunc=trunc)
  lines(SpecDen3$freq[1:lenPlot],SpecDen3$pzgram[1:lenPlot], lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE", trunc=trunc)
SpecDen2 = parzen.wge(y_gen_vpmax_true[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE", trunc=trunc)
ymax = max(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
ymin = min(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
plot(SpecDen$freq[1:lenPlot],SpecDen$pzgram[1:lenPlot], type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - First 40 Days - First Peak - After Update', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq[1:lenPlot],SpecDen2$pzgram[1:lenPlot], lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_adj[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE", trunc=trunc)
  lines(SpecDen3$freq[1:lenPlot],SpecDen3$pzgram[1:lenPlot], lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
CCF = ccf(weather$VPmax..mbar.,weather$Tdew..degC.,plot=FALSE)
plot(CCF$lag, CCF$acf , type = "l", lwd = 6, ylim=c(0.6,0.73),
     main='CCF Comparison - VPMax and TDew', xlab = 'Lag',ylab = 'ACF')
for (sim in 1:numSims){
  CCF2 = ccf(y_gen_vpmax_adj[,sim], y_gen_tdew_adj[,sim],plot=FALSE)
  lines(CCF2$lag, CCF2$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1)

################################################################################

# What is the real cross correlation values between the variables
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
corVD = cor(weather$VPmax..mbar.,weather$Tdew..degC.) # 0.7126456
corTV = cor(weather$T..degC.,weather$VPmax..mbar.) # 0.9679004
corTD = cor(weather$T..degC.,weather$Tdew..degC.) # 0.7827042

# Also, dr sadler wanted us to look at CCF
ccf(weather$VPmax..mbar.,weather$Tdew..degC.) 
crosscf = ccf(weather$T..degC.,weather$VPmax..mbar.) # very high peak at 0
ccf(weather$T..degC.,weather$Tdew..degC.) 

# Generate 100 VPMax
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/vectors_VPMax.csv"
vecsVPMax = read.csv(file_path, header = TRUE)
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
dayStop = 366
numInDay = 6*24
overlap = 20
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
  y_gen_vpmax_noise[y_gen_vpmax_noise < (min(weather$VPmax..mbar.)*0.5)] = min(weather$VPmax..mbar.)*0.5
  y_gen_vpmax_mult[, seed] = y_gen_vpmax_noise
}

# Catherine hasn't gotten to generating a realization yet, so I'll just try with the real data for now
numSim = 100
y_gen_tdew_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_mult) <- paste0("sim", 1:numSims)
for (seed in 1:numSim){
  y_gen_tdew_mult[, seed] = weather$Tdew..degC.
}

# Now generate a vector of correlations
corVec = numeric(numSim)
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  corVec[sim] = cor(y_gen_vpmax_mult[sim],y_gen_tdew_mult[sim])
}
hist(corVec)
abline(corVD) # None of the correlations were as high

# Now try making them closer using that one trick I did before
corVec = numeric(numSim)
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  diff = corVD - cor(y_gen_vpmax_mult[sim],y_gen_tdew_mult[sim])
  alpha <- diff * (1-abs(diff))
  y_gen_tdew_adj <- (1 - alpha) * y_gen_tdew_mult[sim] + alpha * y_gen_vpmax_mult[sim]
  y_gen_vpmax_adj <- (1 - alpha) * y_gen_vpmax_mult[sim] + alpha * y_gen_tdew_adj
  # y_gen_vpmax_adj[y_gen_vpmax_adj<0] = 0 # remove for now, but I'll need to deal with negatives for VPMax
  corVec[sim] = cor(y_gen_vpmax_adj, y_gen_tdew_adj)
}
hist(corVec)
abline(corVD)

# That was too similar, maybe normal diff will work
corVec = numeric(numSim)
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  diff = corVD - cor(y_gen_vpmax_mult[sim],y_gen_tdew_mult[sim])
  # alpha <- diff * (1-abs(diff))
  alpha = diff
  y_gen_tdew_adj <- (1 - alpha) * y_gen_tdew_mult[sim] + alpha * y_gen_vpmax_mult[sim]
  y_gen_vpmax_adj <- (1 - alpha) * y_gen_vpmax_mult[sim] + alpha * y_gen_tdew_adj
  # y_gen_vpmax_adj[y_gen_vpmax_adj<0] = 0 # remove for now, but I'll need to deal with negatives for VPMax
  corVec[sim] = cor(y_gen_vpmax_adj, y_gen_tdew_adj)
}
hist(corVec)

# Maybe try squaring?
corVec = numeric(numSim)
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
  corVec[sim] = cor(y_gen_vpmax_adj, y_gen_tdew_adj)
}
hist(corVec)
abline(v=corVD) # That worked pretty well

# What does the realization look like compared to before?
corVec = numeric(numSim)
y_gen_vpmax_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_adj) <- paste0("sim", 1:numSims)
y_gen_tdew_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_adj) <- paste0("sim", 1:numSims)
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  diff = corVD - cor(y_gen_vpmax_mult[sim],y_gen_tdew_mult[sim])
  # alpha <- diff * (1-abs(diff))
  alpha <- diff * (1-abs(diff))^2
  y_gen_tdew_adj[,sim] <- (1 - alpha) * y_gen_tdew_mult[sim] + alpha * y_gen_vpmax_mult[sim]
  y_gen_vpmax_adj[,sim] <- (1 - alpha) * y_gen_vpmax_mult[sim] + alpha * y_gen_tdew_adj[,sim]
  # y_gen_vpmax_adj[y_gen_vpmax_adj<0] = 0 # remove for now, but I'll need to deal with negatives for VPMax
  corVec[sim] = cor(y_gen_vpmax_adj[,sim], y_gen_tdew_adj[,sim])
}

# Realization Plots
numDays = 366
ymin = min(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_mult[1:(numDays*6*24),1])
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_mult[1:(numDays*6*24),1])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of VPMax - 366 Days', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)  # get default x-axis tick positions
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_vpmax_mult[1:(numDays*6*24),1],col = "blue")
points(y_gen_vpmax_adj[1:(numDays*6*24),1],col = "green")
legend('topleft',legend = c("Real Data", "Simulated Data - Original", "Simulated Data - Nudged"), 
       col = c("red", "blue", "green"),pch=1)
numDays = 10
ymin = min(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_mult[1:(numDays*6*24),1])
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_mult[1:(numDays*6*24),1])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of VPMax - First 10 Days', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_vpmax_mult[1:(numDays*6*24),1],col = "blue")
points(y_gen_vpmax_adj[1:(numDays*6*24),1],col = "green")
legend('top',legend = c("Real Data", "Simulated Data - Original", "Simulated Data - Nudged"), 
       col = c("red", "blue", "green"),pch=1)
startDay = 366/2
numDays = 10
ymin = min(c(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax_mult[(1+6*24*startDay):(6*24*startDay+numDays*6*24),1]))
ymax = max(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax_mult[(1+6*24*startDay):(6*24*startDay+numDays*6*24),1])
plot(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of VPMax - Middle 10 Days', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)
axis(1, at = ticks, labels = round(ticks / 144,0) + startDay)
points(y_gen_vpmax_mult[(1+6*24*startDay):(6*24*startDay+numDays*6*24),1],col = "blue")
points(y_gen_vpmax_adj[(1+6*24*startDay):(6*24*startDay+numDays*6*24),1],col = "green")
legend('topleft',legend = c("Real Data", "Simulated Data - Original", "Simulated Data - Nudged"), 
       col = c("red", "blue", "green"),pch=1)

# Now try with ACF and spectral density
numSims = 20
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
ACF2 = acf(y_gen_vpmax_true, plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6,
     main='ACF Comparison', xlab = 'Lag',ylab = 'ACF')
lines(ACF2$lag, ACF2$acf, lwd = 6, col = "blue")
for (sim in 1:numSims){
  ACF3 = acf(y_gen_vpmax_adj[,sim], plot = "FALSE")
  lines(ACF3$lag, ACF3$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
numDays = 20
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE")
SpecDen2 = parzen.wge(y_gen_vpmax_true[1:(6*24*numDays)], plot = "FALSE")
ymax = max(SpecDen$pzgram,SpecDen2$pzgram)
ymin = min(SpecDen$pzgram,SpecDen2$pzgram)
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - First 20 Days', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_adj[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen3$freq,SpecDen3$pzgram, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
startDay = 180
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
SpecDen2 = parzen.wge(y_gen_vpmax_true[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
ymax = max(SpecDen$pzgram,SpecDen2$pzgram)
ymin = min(SpecDen$pzgram,SpecDen2$pzgram)
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - Middle 20 Days', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_adj[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE")
  lines(SpecDen3$freq,SpecDen3$pzgram, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
trunc = 1000
numDays = 20
lenPlot = 35
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE", trunc=trunc)
SpecDen2 = parzen.wge(y_gen_vpmax_true[1:(6*24*numDays)], plot = "FALSE", trunc=trunc)
ymax = max(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
ymin = min(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
plot(SpecDen$freq[1:lenPlot],SpecDen$pzgram[1:lenPlot], type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - First 40 Days - First Peak', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq[1:lenPlot],SpecDen2$pzgram[1:lenPlot], lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_adj[1:(6*24*numDays),sim], plot = "FALSE", trunc=trunc)
  lines(SpecDen3$freq[1:lenPlot],SpecDen3$pzgram[1:lenPlot], lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE", trunc=trunc)
SpecDen2 = parzen.wge(y_gen_vpmax_true[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE", trunc=trunc)
ymax = max(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
ymin = min(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
plot(SpecDen$freq[1:lenPlot],SpecDen$pzgram[1:lenPlot], type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - First 40 Days - First Peak', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq[1:lenPlot],SpecDen2$pzgram[1:lenPlot], lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_adj[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE", trunc=trunc)
  lines(SpecDen3$freq[1:lenPlot],SpecDen3$pzgram[1:lenPlot], lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
# That actually worked out quite well

# Let's try with some CCF plots
numSims = 20
CCF = ccf(weather$VPmax..mbar.,weather$Tdew..degC.,plot=FALSE)
plot(CCF$lag, CCF$acf , type = "l", lwd = 6, 
     main='CCF Comparison - VPMax and TDew', xlab = 'Lag',ylab = 'ACF')
for (sim in 1:numSims){
  CCF2 = ccf(y_gen_vpmax_adj[,sim], y_gen_tdew_adj[,sim],plot=FALSE)
  lines(CCF2$lag, CCF2$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1)
# That doesn't look great.  Only the lag 0 values are close

# What if I try some ahead
y_gen_vpmax_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_adj) <- paste0("sim", 1:numSims)
y_gen_tdew_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_adj) <- paste0("sim", 1:numSims)
numAhead = 40
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  
  # No lag
  diff = corVD - cor(y_gen_vpmax_mult[sim],y_gen_tdew_mult[sim])
  # alpha <- diff * (1-abs(diff))
  alpha <- diff * (1-abs(diff))^2
  y_gen_tdew_adj[,sim] <- (1 - alpha) * y_gen_tdew_mult[sim] + alpha * y_gen_vpmax_mult[sim]
  y_gen_vpmax_adj[,sim] <- (1 - alpha) * y_gen_vpmax_mult[sim] + alpha * y_gen_tdew_adj[,sim]
  # y_gen_vpmax_adj[y_gen_vpmax_adj<0] = 0 # remove for now, but I'll need to deal with negatives for VPMax
  
  # lag
  for (ahead in 1:numAhead){
    corr = cor(weather$VPmax..mbar.[(ahead+1):(6*24*366)],weather$Tdew..degC.[1:(6*24*366-ahead)])
    diff = corr - cor(y_gen_vpmax_mult[(ahead+1):(6*24*366),sim],y_gen_tdew_mult[1:(6*24*366-ahead),sim])
    # alpha <- diff * (1-abs(diff))
    alpha <- diff * (1-abs(diff))^2
    y_gen_tdew_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_tdew_mult[1:(6*24*366-ahead),sim] + alpha * y_gen_vpmax_mult[(ahead+1):(6*24*366),sim]
    y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_vpmax_mult[(ahead+1):(6*24*366),sim] + alpha * y_gen_tdew_adj[1:(6*24*366-ahead),sim]
    # y_gen_vpmax_adj[y_gen_vpmax_adj<0] = 0 # remove for now, but I'll need to deal with negatives for VPMax
  }
}
# This seems to shift the max to the last number

# What if I try shifting each way?
y_gen_vpmax_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_adj) <- paste0("sim", 1:numSims)
y_gen_tdew_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_adj) <- paste0("sim", 1:numSims)
numAhead = 1
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  
  # Initialize
  y_gen_vpmax_adj[sim] = y_gen_vpmax_mult[sim]
  y_gen_tdew_adj[sim] = y_gen_tdew_mult[sim]
  
  # lag
  for (ahead in numAhead:1){
    corr = cor(weather$VPmax..mbar.[(ahead+1):(6*24*366)],weather$Tdew..degC.[1:(6*24*366-ahead)])
    diff = corr - cor(y_gen_vpmax_adj[(ahead+1):(6*24*366),sim],y_gen_tdew_adj[1:(6*24*366-ahead),sim])
    alpha <- diff * (1-abs(diff))^2
    y_gen_tdew_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_tdew_mult[1:(6*24*366-ahead),sim] + alpha * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim]
    y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_tdew_adj[1:(6*24*366-ahead),sim]
    
    corr = cor(weather$VPmax..mbar.[1:(6*24*366-ahead)],weather$Tdew..degC.[(ahead+1):(6*24*366)])
    diff = corr - cor(y_gen_vpmax_adj[1:(6*24*366-ahead),sim],y_gen_tdew_adj[(ahead+1):(6*24*366),sim])
    alpha <- diff * (1-abs(diff))^2
    y_gen_tdew_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_tdew_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_vpmax_adj[1:(6*24*366-ahead),sim]
    y_gen_vpmax_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_vpmax_adj[1:(6*24*366-ahead),sim] + alpha * y_gen_tdew_adj[(ahead+1):(6*24*366),sim]
  }
  
  # No lag
  diff = corVD - cor(y_gen_vpmax_mult[sim],y_gen_tdew_mult[sim])
  alpha <- diff * (1-abs(diff))^2
  y_gen_tdew_adj[,sim] <- (1 - alpha) * y_gen_tdew_mult[sim] + alpha * y_gen_vpmax_adj[sim]
  y_gen_vpmax_adj[,sim] <- (1 - alpha) * y_gen_vpmax_adj[sim] + alpha * y_gen_tdew_adj[,sim]
  
}
numSims = 20
CCF = ccf(weather$VPmax..mbar.,weather$Tdew..degC.,plot=FALSE)
plot(CCF$lag, CCF$acf , type = "l", lwd = 6, ylim=c(0.5,0.73),
     main='CCF Comparison - VPMax and TDew', xlab = 'Lag',ylab = 'ACF')
for (sim in 1:numSims){
  CCF2 = ccf(y_gen_vpmax_adj[,sim], y_gen_tdew_adj[,sim],plot=FALSE)
  lines(CCF2$lag, CCF2$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1) # That did something different, but it's still not good
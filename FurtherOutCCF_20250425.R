# Demo

# 20 Simulations, original CCF
library(tswge)
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
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
numSims = 20
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
y_gen_vpmax_adj = y_gen_vpmax_mult
y_gen_tdew_adj = y_gen_tdew_mult
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
ACF2 = acf(y_gen_vpmax_true, plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6,
     main='ACF Comparison - Original simulations', xlab = 'Lag',ylab = 'ACF')
lines(ACF2$lag, ACF2$acf, lwd = 6, col = "blue")
for (sim in 1:numSims){
  ACF3 = acf(y_gen_vpmax_adj[,sim], plot = "FALSE")
  lines(ACF3$lag, ACF3$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
CCF = ccf(weather$VPmax..mbar.,weather$Tdew..degC.,plot=FALSE)
plot(CCF$lag, CCF$acf , type = "l", lwd = 6, ylim=c(0.55,0.73),
     main='CCF Comparison - VPMax and TDew - Original simulations', xlab = 'Lag',ylab = 'ACF')
for (sim in 1:numSims){
  CCF2 = ccf(y_gen_vpmax_adj[,sim], y_gen_tdew_adj[,sim],plot=FALSE)
  lines(CCF2$lag, CCF2$acf, lwd = 2, col = "red")
}
legend('topright',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1) 

# What I had originally, good ACF plot but bad CCF plot
y_gen_vpmax_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_adj) <- paste0("sim", 1:numSims)
y_gen_tdew_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_adj) <- paste0("sim", 1:numSims)
numAhead = 1
gapAhead = 0
# gapAhead = 40 / numAhead
numIter = 1
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  
  # Initialize
  y_gen_vpmax_adj[sim] = y_gen_vpmax_mult[sim]
  y_gen_tdew_adj[sim] = y_gen_tdew_mult[sim]
  
  # Try doing the nudging multiple times
  for (iter in 1:numIter){
    
    # No lag
    corr = cor(weather$VPmax..mbar.,weather$Tdew..degC.)
    diff = corr - cor(y_gen_vpmax_adj[sim],y_gen_tdew_adj[sim])
    alpha <- diff * (1-abs(diff))^2
    y_gen_tdew_adj[,sim] <- (1 - alpha) * y_gen_tdew_adj[sim] + alpha * y_gen_vpmax_adj[sim]
    y_gen_vpmax_adj[,sim] <- (1 - alpha) * y_gen_vpmax_adj[sim] + alpha * y_gen_tdew_adj[,sim]
    
    # lag
    for (aheadCounter in numAhead:1){
      ahead = aheadCounter * gapAhead
      
      # Positive lag
      corr = cor(weather$VPmax..mbar.[(ahead+1):(6*24*366)],weather$Tdew..degC.[1:(6*24*366-ahead)])
      diff = corr - cor(y_gen_vpmax_adj[(ahead+1):(6*24*366),sim],y_gen_tdew_adj[1:(6*24*366-ahead),sim])
      alpha <- diff * (1-abs(diff))^2
      y_gen_tdew_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_tdew_mult[1:(6*24*366-ahead),sim] + alpha * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim]
      y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_tdew_adj[1:(6*24*366-ahead),sim]
      
      # Negative lag
      corr = cor(weather$VPmax..mbar.[1:(6*24*366-ahead)],weather$Tdew..degC.[(ahead+1):(6*24*366)])
      diff = corr - cor(y_gen_vpmax_adj[1:(6*24*366-ahead),sim],y_gen_tdew_adj[(ahead+1):(6*24*366),sim])
      alpha <- diff * (1-abs(diff))^2
      y_gen_tdew_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_tdew_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_vpmax_adj[1:(6*24*366-ahead),sim]
      y_gen_vpmax_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_vpmax_adj[1:(6*24*366-ahead),sim] + alpha * y_gen_tdew_adj[(ahead+1):(6*24*366),sim]
    }
    
    # Make sure it isn't too low
    y_gen_vpmax_adj[y_gen_vpmax_adj < (min(weather$VPmax..mbar.)*0.5)] = min(weather$VPmax..mbar.)*0.5
    
  }
}
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
ACF2 = acf(y_gen_vpmax_true, plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6,
     main='ACF Comparison - 3 nudges at 0 lag', xlab = 'Lag',ylab = 'ACF')
lines(ACF2$lag, ACF2$acf, lwd = 6, col = "blue")
for (sim in 1:numSims){
  ACF3 = acf(y_gen_vpmax_adj[,sim], plot = "FALSE")
  lines(ACF3$lag, ACF3$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
CCF = ccf(weather$VPmax..mbar.,weather$Tdew..degC.,plot=FALSE)
plot(CCF$lag, CCF$acf , type = "l", lwd = 6, ylim=c(0.69,0.73),
     main='CCF Comparison - VPMax and TDew - 3 nudges at 0 lag', xlab = 'Lag',ylab = 'ACF')
for (sim in 1:numSims){
  CCF2 = ccf(y_gen_vpmax_adj[,sim], y_gen_tdew_adj[,sim],plot=FALSE)
  lines(CCF2$lag, CCF2$acf, lwd = 2, col = "red")
}
legend('topright',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1) 

# 20 nudges at every 5 lags, Bad ACF plot but good CCF plot
y_gen_vpmax_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_adj) <- paste0("sim", 1:numSims)
y_gen_tdew_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_adj) <- paste0("sim", 1:numSims)
numAhead = 8
# gapAhead = 0
gapAhead = 40 / numAhead
numIter = 20
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  
  # Initialize
  y_gen_vpmax_adj[sim] = y_gen_vpmax_mult[sim]
  y_gen_tdew_adj[sim] = y_gen_tdew_mult[sim]
  
  # Try doing the nudging multiple times
  for (iter in 1:numIter){
    
    # No lag
    corr = cor(weather$VPmax..mbar.,weather$Tdew..degC.)
    diff = corr - cor(y_gen_vpmax_adj[sim],y_gen_tdew_adj[sim])
    alpha <- diff * (1-abs(diff))^2
    y_gen_tdew_adj[,sim] <- (1 - alpha) * y_gen_tdew_adj[sim] + alpha * y_gen_vpmax_adj[sim]
    y_gen_vpmax_adj[,sim] <- (1 - alpha) * y_gen_vpmax_adj[sim] + alpha * y_gen_tdew_adj[,sim]
    
    # lag
    for (aheadCounter in numAhead:1){
      ahead = aheadCounter * gapAhead
      
      # Positive lag
      corr = cor(weather$VPmax..mbar.[(ahead+1):(6*24*366)],weather$Tdew..degC.[1:(6*24*366-ahead)])
      diff = corr - cor(y_gen_vpmax_adj[(ahead+1):(6*24*366),sim],y_gen_tdew_adj[1:(6*24*366-ahead),sim])
      alpha <- diff * (1-abs(diff))^2
      y_gen_tdew_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_tdew_mult[1:(6*24*366-ahead),sim] + alpha * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim]
      y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_tdew_adj[1:(6*24*366-ahead),sim]
      
      # Negative lag
      corr = cor(weather$VPmax..mbar.[1:(6*24*366-ahead)],weather$Tdew..degC.[(ahead+1):(6*24*366)])
      diff = corr - cor(y_gen_vpmax_adj[1:(6*24*366-ahead),sim],y_gen_tdew_adj[(ahead+1):(6*24*366),sim])
      alpha <- diff * (1-abs(diff))^2
      y_gen_tdew_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_tdew_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_vpmax_adj[1:(6*24*366-ahead),sim]
      y_gen_vpmax_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_vpmax_adj[1:(6*24*366-ahead),sim] + alpha * y_gen_tdew_adj[(ahead+1):(6*24*366),sim]
    }
    
    # Make sure it isn't too low
    y_gen_vpmax_adj[y_gen_vpmax_adj < (min(weather$VPmax..mbar.)*0.5)] = min(weather$VPmax..mbar.)*0.5
    
  }
}
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
ACF2 = acf(y_gen_vpmax_true, plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6,
     main='ACF Comparison - 20 nudges at every 5 lags', xlab = 'Lag',ylab = 'ACF')
lines(ACF2$lag, ACF2$acf, lwd = 6, col = "blue")
for (sim in 1:numSims){
  ACF3 = acf(y_gen_vpmax_adj[,sim], plot = "FALSE")
  lines(ACF3$lag, ACF3$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
CCF = ccf(weather$VPmax..mbar.,weather$Tdew..degC.,plot=FALSE)
plot(CCF$lag, CCF$acf , type = "l", lwd = 6, ylim=c(0.69,0.73),
     main='CCF Comparison - VPMax and TDew - 20 nudges at every 5 lags', xlab = 'Lag',ylab = 'ACF')
for (sim in 1:numSims){
  CCF2 = ccf(y_gen_vpmax_adj[,sim], y_gen_tdew_adj[,sim],plot=FALSE)
  lines(CCF2$lag, CCF2$acf, lwd = 2, col = "red")
}
legend('topright',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1) 

# 3 nudges at every 40 lags, OK ACF plot and OK CCF plot
y_gen_vpmax_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_adj) <- paste0("sim", 1:numSims)
y_gen_tdew_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_adj) <- paste0("sim", 1:numSims)
numAhead = 1
# gapAhead = 0
gapAhead = 40 / numAhead
numIter = 3
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  
  # Initialize
  y_gen_vpmax_adj[sim] = y_gen_vpmax_mult[sim]
  y_gen_tdew_adj[sim] = y_gen_tdew_mult[sim]
  
  # Try doing the nudging multiple times
  for (iter in 1:numIter){
    
    # No lag
    corr = cor(weather$VPmax..mbar.,weather$Tdew..degC.)
    diff = corr - cor(y_gen_vpmax_adj[sim],y_gen_tdew_adj[sim])
    alpha <- diff * (1-abs(diff))^2
    y_gen_tdew_adj[,sim] <- (1 - alpha) * y_gen_tdew_adj[sim] + alpha * y_gen_vpmax_adj[sim]
    y_gen_vpmax_adj[,sim] <- (1 - alpha) * y_gen_vpmax_adj[sim] + alpha * y_gen_tdew_adj[,sim]
    
    # lag
    for (aheadCounter in numAhead:1){
      ahead = aheadCounter * gapAhead
      
      # Positive lag
      corr = cor(weather$VPmax..mbar.[(ahead+1):(6*24*366)],weather$Tdew..degC.[1:(6*24*366-ahead)])
      diff = corr - cor(y_gen_vpmax_adj[(ahead+1):(6*24*366),sim],y_gen_tdew_adj[1:(6*24*366-ahead),sim])
      alpha <- diff * (1-abs(diff))^2
      y_gen_tdew_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_tdew_mult[1:(6*24*366-ahead),sim] + alpha * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim]
      y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_tdew_adj[1:(6*24*366-ahead),sim]
      
      # Negative lag
      corr = cor(weather$VPmax..mbar.[1:(6*24*366-ahead)],weather$Tdew..degC.[(ahead+1):(6*24*366)])
      diff = corr - cor(y_gen_vpmax_adj[1:(6*24*366-ahead),sim],y_gen_tdew_adj[(ahead+1):(6*24*366),sim])
      alpha <- diff * (1-abs(diff))^2
      y_gen_tdew_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_tdew_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_vpmax_adj[1:(6*24*366-ahead),sim]
      y_gen_vpmax_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_vpmax_adj[1:(6*24*366-ahead),sim] + alpha * y_gen_tdew_adj[(ahead+1):(6*24*366),sim]
    }
    
    # Make sure it isn't too low
    y_gen_vpmax_adj[y_gen_vpmax_adj < (min(weather$VPmax..mbar.)*0.5)] = min(weather$VPmax..mbar.)*0.5
    
  }
}
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
ACF2 = acf(y_gen_vpmax_true, plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6,
     main='ACF Comparison - 3 nudges at every 40 lags', xlab = 'Lag',ylab = 'ACF')
lines(ACF2$lag, ACF2$acf, lwd = 6, col = "blue")
for (sim in 1:numSims){
  ACF3 = acf(y_gen_vpmax_adj[,sim], plot = "FALSE")
  lines(ACF3$lag, ACF3$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
CCF = ccf(weather$VPmax..mbar.,weather$Tdew..degC.,plot=FALSE)
plot(CCF$lag, CCF$acf , type = "l", lwd = 6, ylim=c(0.69,0.73),
     main='CCF Comparison - VPMax and TDew - 3 nudges at every 40 lags', xlab = 'Lag',ylab = 'ACF')
for (sim in 1:numSims){
  CCF2 = ccf(y_gen_vpmax_adj[,sim], y_gen_tdew_adj[,sim],plot=FALSE)
  lines(CCF2$lag, CCF2$acf, lwd = 2, col = "red")
}
legend('topright',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1) 

################################################################################

# Last time it seemed like nudging +1/-1 for the lag caused lag 0 to go too high.
# I'm wondering if we could do something like +5/-5, or +10/-10, etc. to get around that problem.

# Get 20 simulations for testing CCFs
library(tswge)
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
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
numSims = 20
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

# Try messing with the CCF with more that steps of 1 for the lag
# What if I try shifting each way?
y_gen_vpmax_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_adj) <- paste0("sim", 1:numSims)
y_gen_tdew_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_adj) <- paste0("sim", 1:numSims)
numAhead = 1
gapAhead = 40
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  
  # Initialize
  y_gen_vpmax_adj[sim] = y_gen_vpmax_mult[sim]
  y_gen_tdew_adj[sim] = y_gen_tdew_mult[sim]
  
  # lag
  for (aheadCounter in numAhead:1){
    ahead = aheadCounter * gapAhead
    
    # Positive lag
    corr = cor(weather$VPmax..mbar.[(ahead+1):(6*24*366)],weather$Tdew..degC.[1:(6*24*366-ahead)])
    diff = corr - cor(y_gen_vpmax_adj[(ahead+1):(6*24*366),sim],y_gen_tdew_adj[1:(6*24*366-ahead),sim])
    alpha <- diff * (1-abs(diff))^2
    y_gen_tdew_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_tdew_mult[1:(6*24*366-ahead),sim] + alpha * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim]
    y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_tdew_adj[1:(6*24*366-ahead),sim]
    
    # Negative lag
    corr = cor(weather$VPmax..mbar.[1:(6*24*366-ahead)],weather$Tdew..degC.[(ahead+1):(6*24*366)])
    diff = corr - cor(y_gen_vpmax_adj[1:(6*24*366-ahead),sim],y_gen_tdew_adj[(ahead+1):(6*24*366),sim])
    alpha <- diff * (1-abs(diff))^2
    y_gen_tdew_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_tdew_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_vpmax_adj[1:(6*24*366-ahead),sim]
    y_gen_vpmax_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_vpmax_adj[1:(6*24*366-ahead),sim] + alpha * y_gen_tdew_adj[(ahead+1):(6*24*366),sim]
  }
  
  # No lag
  corr = cor(weather$VPmax..mbar.,weather$Tdew..degC.)
  diff = corr - cor(y_gen_vpmax_adj[sim],y_gen_tdew_adj[sim])
  alpha <- diff * (1-abs(diff))^2
  y_gen_tdew_adj[,sim] <- (1 - alpha) * y_gen_tdew_adj[sim] + alpha * y_gen_vpmax_adj[sim]
  y_gen_vpmax_adj[,sim] <- (1 - alpha) * y_gen_vpmax_adj[sim] + alpha * y_gen_tdew_adj[,sim]
  
}
CCF = ccf(weather$VPmax..mbar.,weather$Tdew..degC.,plot=FALSE)
plot(CCF$lag, CCF$acf , type = "l", lwd = 6, ylim=c(0.69,0.74),
     main='CCF Comparison - VPMax and TDew', xlab = 'Lag',ylab = 'ACF')
for (sim in 1:numSims){
  CCF2 = ccf(y_gen_vpmax_adj[,sim], y_gen_tdew_adj[,sim],plot=FALSE)
  lines(CCF2$lag, CCF2$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1) # This does well around -40 and 0, but 40 it jumps too high

# What does one round of this look like for the positive lag
y_gen_vpmax_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_adj) <- paste0("sim", 1:numSims)
y_gen_tdew_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_adj) <- paste0("sim", 1:numSims)
numAhead = 1
gapAhead = 40
sim = 1
y_gen_vpmax_adj[sim] = y_gen_vpmax_mult[sim]
y_gen_tdew_adj[sim] = y_gen_tdew_mult[sim]
corr = cor(weather$VPmax..mbar.[(ahead+1):(6*24*366)],weather$Tdew..degC.[1:(6*24*366-ahead)]) # 0.7
corrCurr = cor(y_gen_vpmax_adj[(ahead+1):(6*24*366),sim],y_gen_tdew_adj[1:(6*24*366-ahead),sim]) # 0.58
diff = corr - cor(y_gen_vpmax_adj[(ahead+1):(6*24*366),sim],y_gen_tdew_adj[1:(6*24*366-ahead),sim])
alpha <- diff * (1-abs(diff))^2
y_gen_tdew_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_tdew_mult[1:(6*24*366-ahead),sim] + alpha * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim]
y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_tdew_adj[1:(6*24*366-ahead),sim]
corrCurr = cor(y_gen_vpmax_adj[(ahead+1):(6*24*366),sim],y_gen_tdew_adj[1:(6*24*366-ahead),sim]) # 0.7
corr = cor(weather$VPmax..mbar.[1:(6*24*366-ahead)],weather$Tdew..degC.[(ahead+1):(6*24*366)])
diff = corr - cor(y_gen_vpmax_adj[1:(6*24*366-ahead),sim],y_gen_tdew_adj[(ahead+1):(6*24*366),sim])
alpha <- diff * (1-abs(diff))^2
y_gen_tdew_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_tdew_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_vpmax_adj[1:(6*24*366-ahead),sim]
y_gen_vpmax_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_vpmax_adj[1:(6*24*366-ahead),sim] + alpha * y_gen_tdew_adj[(ahead+1):(6*24*366),sim]
corrCurr = cor(y_gen_vpmax_adj[(ahead+1):(6*24*366),sim],y_gen_tdew_adj[1:(6*24*366-ahead),sim]) 
# Hm, now the correlation is higher, even though the nudging happened 80 points in the past

# Maybe an iterative approach to make sure they all are close
y_gen_vpmax_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_adj) <- paste0("sim", 1:numSims)
y_gen_tdew_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_adj) <- paste0("sim", 1:numSims)
numAhead = 8
gapAhead = 5
numIter = 20
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  
  # Initialize
  y_gen_vpmax_adj[sim] = y_gen_vpmax_mult[sim]
  y_gen_tdew_adj[sim] = y_gen_tdew_mult[sim]
  
  # Try doing the nudging multiple times
  for (iter in 1:numIter){
    
    # lag
    for (aheadCounter in numAhead:1){
      ahead = aheadCounter * gapAhead
      
      # Positive lag
      corr = cor(weather$VPmax..mbar.[(ahead+1):(6*24*366)],weather$Tdew..degC.[1:(6*24*366-ahead)])
      diff = corr - cor(y_gen_vpmax_adj[(ahead+1):(6*24*366),sim],y_gen_tdew_adj[1:(6*24*366-ahead),sim])
      alpha <- diff * (1-abs(diff))^2
      y_gen_tdew_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_tdew_mult[1:(6*24*366-ahead),sim] + alpha * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim]
      y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_tdew_adj[1:(6*24*366-ahead),sim]
      
      # Negative lag
      corr = cor(weather$VPmax..mbar.[1:(6*24*366-ahead)],weather$Tdew..degC.[(ahead+1):(6*24*366)])
      diff = corr - cor(y_gen_vpmax_adj[1:(6*24*366-ahead),sim],y_gen_tdew_adj[(ahead+1):(6*24*366),sim])
      alpha <- diff * (1-abs(diff))^2
      y_gen_tdew_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_tdew_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_vpmax_adj[1:(6*24*366-ahead),sim]
      y_gen_vpmax_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_vpmax_adj[1:(6*24*366-ahead),sim] + alpha * y_gen_tdew_adj[(ahead+1):(6*24*366),sim]
    }
    
    # No lag
    corr = cor(weather$VPmax..mbar.,weather$Tdew..degC.)
    diff = corr - cor(y_gen_vpmax_adj[sim],y_gen_tdew_adj[sim])
    alpha <- diff * (1-abs(diff))^2
    y_gen_tdew_adj[,sim] <- (1 - alpha) * y_gen_tdew_adj[sim] + alpha * y_gen_vpmax_adj[sim]
    y_gen_vpmax_adj[,sim] <- (1 - alpha) * y_gen_vpmax_adj[sim] + alpha * y_gen_tdew_adj[,sim]
    
  }
  
}
CCF = ccf(weather$VPmax..mbar.,weather$Tdew..degC.,plot=FALSE)
plot(CCF$lag, CCF$acf , type = "l", lwd = 6, ylim=c(0.69,0.73),
     main='CCF Comparison - VPMax and TDew', xlab = 'Lag',ylab = 'ACF')
for (sim in 1:numSims){
  CCF2 = ccf(y_gen_vpmax_adj[,sim], y_gen_tdew_adj[,sim],plot=FALSE)
  lines(CCF2$lag, CCF2$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1) # It looks like doing it with jumps of 5, for 20 iterations works pretty well

# Now let's see how much this messed up a single realization
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
# It noticeably dampens the maxes, and may actually have it drop below the minimum, but it doesn't look terrible

# Deal with minimum
y_gen_vpmax_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_adj) <- paste0("sim", 1:numSims)
y_gen_tdew_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_adj) <- paste0("sim", 1:numSims)
numAhead = 8
gapAhead = 5
numIter = 20
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  
  # Initialize
  y_gen_vpmax_adj[sim] = y_gen_vpmax_mult[sim]
  y_gen_tdew_adj[sim] = y_gen_tdew_mult[sim]
  
  # Try doing the nudging multiple times
  for (iter in 1:numIter){
    
    # lag
    for (aheadCounter in numAhead:1){
      ahead = aheadCounter * gapAhead
      
      # Positive lag
      corr = cor(weather$VPmax..mbar.[(ahead+1):(6*24*366)],weather$Tdew..degC.[1:(6*24*366-ahead)])
      diff = corr - cor(y_gen_vpmax_adj[(ahead+1):(6*24*366),sim],y_gen_tdew_adj[1:(6*24*366-ahead),sim])
      alpha <- diff * (1-abs(diff))^2
      y_gen_tdew_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_tdew_mult[1:(6*24*366-ahead),sim] + alpha * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim]
      y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_tdew_adj[1:(6*24*366-ahead),sim]
      
      # Negative lag
      corr = cor(weather$VPmax..mbar.[1:(6*24*366-ahead)],weather$Tdew..degC.[(ahead+1):(6*24*366)])
      diff = corr - cor(y_gen_vpmax_adj[1:(6*24*366-ahead),sim],y_gen_tdew_adj[(ahead+1):(6*24*366),sim])
      alpha <- diff * (1-abs(diff))^2
      y_gen_tdew_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_tdew_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_vpmax_adj[1:(6*24*366-ahead),sim]
      y_gen_vpmax_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_vpmax_adj[1:(6*24*366-ahead),sim] + alpha * y_gen_tdew_adj[(ahead+1):(6*24*366),sim]
    }
    
    # No lag
    corr = cor(weather$VPmax..mbar.,weather$Tdew..degC.)
    diff = corr - cor(y_gen_vpmax_adj[sim],y_gen_tdew_adj[sim])
    alpha <- diff * (1-abs(diff))^2
    y_gen_tdew_adj[,sim] <- (1 - alpha) * y_gen_tdew_adj[sim] + alpha * y_gen_vpmax_adj[sim]
    y_gen_vpmax_adj[,sim] <- (1 - alpha) * y_gen_vpmax_adj[sim] + alpha * y_gen_tdew_adj[,sim]
    
    # Make sure it isn't too low
    y_gen_vpmax_adj[y_gen_vpmax_adj < (min(weather$VPmax..mbar.)*0.5)] = min(weather$VPmax..mbar.)*0.5
    
  }
  
}
CCF = ccf(weather$VPmax..mbar.,weather$Tdew..degC.,plot=FALSE)
plot(CCF$lag, CCF$acf , type = "l", lwd = 6, ylim=c(0.69,0.73),
     main='CCF Comparison - VPMax and TDew', xlab = 'Lag',ylab = 'ACF')
for (sim in 1:numSims){
  CCF2 = ccf(y_gen_vpmax_adj[,sim], y_gen_tdew_adj[,sim],plot=FALSE)
  lines(CCF2$lag, CCF2$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1) # It looks like doing it with jumps of 5, for 20 iterations works pretty well
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

# Now look at ACF and spectral density
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
lenPlot = 35
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE", trunc=trunc)
SpecDen2 = parzen.wge(y_gen_vpmax_true[1:(6*24*numDays)], plot = "FALSE", trunc=trunc)
ymax = max(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
ymin = min(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
plot(SpecDen$freq[1:lenPlot],SpecDen$pzgram[1:lenPlot], type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - First 20 Days - First Peak', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
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
     main='Spectral Density Comparison - First 20 Days - First Peak', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq[1:lenPlot],SpecDen2$pzgram[1:lenPlot], lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_adj[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE", trunc=trunc)
  lines(SpecDen3$freq[1:lenPlot],SpecDen3$pzgram[1:lenPlot], lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
# ACF isn't quite as good, like it's a better ACF than the data, spectral densities are good though

# maybe try with less gap
y_gen_vpmax_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_adj) <- paste0("sim", 1:numSims)
y_gen_tdew_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_adj) <- paste0("sim", 1:numSims)
numAhead = 1
gapAhead = 20
# gapAhead = 40 / numAhead
numIter = 5
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  
  # Initialize
  y_gen_vpmax_adj[sim] = y_gen_vpmax_mult[sim]
  y_gen_tdew_adj[sim] = y_gen_tdew_mult[sim]
  
  # Try doing the nudging multiple times
  for (iter in 1:numIter){
    
    # lag
    for (aheadCounter in numAhead:1){
      ahead = aheadCounter * gapAhead
      
      # Positive lag
      corr = cor(weather$VPmax..mbar.[(ahead+1):(6*24*366)],weather$Tdew..degC.[1:(6*24*366-ahead)])
      diff = corr - cor(y_gen_vpmax_adj[(ahead+1):(6*24*366),sim],y_gen_tdew_adj[1:(6*24*366-ahead),sim])
      alpha <- diff * (1-abs(diff))^2
      y_gen_tdew_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_tdew_mult[1:(6*24*366-ahead),sim] + alpha * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim]
      y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_tdew_adj[1:(6*24*366-ahead),sim]
      
      # Negative lag
      corr = cor(weather$VPmax..mbar.[1:(6*24*366-ahead)],weather$Tdew..degC.[(ahead+1):(6*24*366)])
      diff = corr - cor(y_gen_vpmax_adj[1:(6*24*366-ahead),sim],y_gen_tdew_adj[(ahead+1):(6*24*366),sim])
      alpha <- diff * (1-abs(diff))^2
      y_gen_tdew_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_tdew_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_vpmax_adj[1:(6*24*366-ahead),sim]
      y_gen_vpmax_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_vpmax_adj[1:(6*24*366-ahead),sim] + alpha * y_gen_tdew_adj[(ahead+1):(6*24*366),sim]
    }
    
    # No lag
    corr = cor(weather$VPmax..mbar.,weather$Tdew..degC.)
    diff = corr - cor(y_gen_vpmax_adj[sim],y_gen_tdew_adj[sim])
    alpha <- diff * (1-abs(diff))^2
    y_gen_tdew_adj[,sim] <- (1 - alpha) * y_gen_tdew_adj[sim] + alpha * y_gen_vpmax_adj[sim]
    y_gen_vpmax_adj[,sim] <- (1 - alpha) * y_gen_vpmax_adj[sim] + alpha * y_gen_tdew_adj[,sim]
    
    # Make sure it isn't too low
    y_gen_vpmax_adj[y_gen_vpmax_adj < (min(weather$VPmax..mbar.)*0.5)] = min(weather$VPmax..mbar.)*0.5
    
  }
  
}
CCF = ccf(weather$VPmax..mbar.,weather$Tdew..degC.,plot=FALSE)
plot(CCF$lag, CCF$acf , type = "l", lwd = 6, ylim=c(0.69,0.73),
     main='CCF Comparison - VPMax and TDew', xlab = 'Lag',ylab = 'ACF')
for (sim in 1:numSims){
  CCF2 = ccf(y_gen_vpmax_adj[,sim], y_gen_tdew_adj[,sim],plot=FALSE)
  lines(CCF2$lag, CCF2$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1) # It looks like doing it with jumps of 5, for 20 iterations works pretty well
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
# There does seem to be a trade-off between a good CCF and an ACF that's too good

# Maybe try doing the no lag part first
y_gen_vpmax_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_adj) <- paste0("sim", 1:numSims)
y_gen_tdew_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_adj) <- paste0("sim", 1:numSims)
numAhead = 1
# gapAhead = 0
gapAhead = 40 / numAhead
numIter = 3
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  
  # Initialize
  y_gen_vpmax_adj[sim] = y_gen_vpmax_mult[sim]
  y_gen_tdew_adj[sim] = y_gen_tdew_mult[sim]
  
  # Try doing the nudging multiple times
  for (iter in 1:numIter){
    
    # No lag
    corr = cor(weather$VPmax..mbar.,weather$Tdew..degC.)
    diff = corr - cor(y_gen_vpmax_adj[sim],y_gen_tdew_adj[sim])
    alpha <- diff * (1-abs(diff))^2
    y_gen_tdew_adj[,sim] <- (1 - alpha) * y_gen_tdew_adj[sim] + alpha * y_gen_vpmax_adj[sim]
    y_gen_vpmax_adj[,sim] <- (1 - alpha) * y_gen_vpmax_adj[sim] + alpha * y_gen_tdew_adj[,sim]
    
    # lag
    for (aheadCounter in numAhead:1){
      ahead = aheadCounter * gapAhead
      
      # Positive lag
      corr = cor(weather$VPmax..mbar.[(ahead+1):(6*24*366)],weather$Tdew..degC.[1:(6*24*366-ahead)])
      diff = corr - cor(y_gen_vpmax_adj[(ahead+1):(6*24*366),sim],y_gen_tdew_adj[1:(6*24*366-ahead),sim])
      alpha <- diff * (1-abs(diff))^2
      y_gen_tdew_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_tdew_mult[1:(6*24*366-ahead),sim] + alpha * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim]
      y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_tdew_adj[1:(6*24*366-ahead),sim]
      
      # Negative lag
      corr = cor(weather$VPmax..mbar.[1:(6*24*366-ahead)],weather$Tdew..degC.[(ahead+1):(6*24*366)])
      diff = corr - cor(y_gen_vpmax_adj[1:(6*24*366-ahead),sim],y_gen_tdew_adj[(ahead+1):(6*24*366),sim])
      alpha <- diff * (1-abs(diff))^2
      y_gen_tdew_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_tdew_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_vpmax_adj[1:(6*24*366-ahead),sim]
      y_gen_vpmax_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_vpmax_adj[1:(6*24*366-ahead),sim] + alpha * y_gen_tdew_adj[(ahead+1):(6*24*366),sim]
    }
    
    # Make sure it isn't too low
    y_gen_vpmax_adj[y_gen_vpmax_adj < (min(weather$VPmax..mbar.)*0.5)] = min(weather$VPmax..mbar.)*0.5
    
  }
  
}
CCF = ccf(weather$VPmax..mbar.,weather$Tdew..degC.,plot=FALSE)
plot(CCF$lag, CCF$acf , type = "l", lwd = 6, ylim=c(0.69,0.73),
     main='CCF Comparison - VPMax and TDew', xlab = 'Lag',ylab = 'ACF')
for (sim in 1:numSims){
  CCF2 = ccf(y_gen_vpmax_adj[,sim], y_gen_tdew_adj[,sim],plot=FALSE)
  lines(CCF2$lag, CCF2$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1) # It looks like doing it with jumps of 5, for 20 iterations works pretty well
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
# Hm, the trade-off is still there
# numAhead = 1, gapAhead = 40, and numIter = 3 seems like a good middle point
# numAhead = 1, gapAhead = 0, and numIter = 1 seems like it produces the original way
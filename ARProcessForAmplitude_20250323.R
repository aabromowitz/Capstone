# In this notebook, I'm going to try to use an AR process to generate the mean and amplitude of the signal part.
# Maybe the variance in the noise as well.

## Demo

# Plot amplitude
# I wasn't sure a good way to model amplitude
library(tswge)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/vectors4.csv"
vecs = read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
A1vec = vecs$A1
absA1vec = abs(A1vec)
plot(absA1vec) # variance seems smaller in winter and larger in summer, not sure it follows AR pattern
ma.smooth.wge(absA1vec,order=75) # it looks kinda like a sine curve if you smooth it out a lot
t = (1:366)
res = absA1vec - 6.06861*sin(2*pi*(1/(366*2))*t)
plot(res) # if you model with a sine curve though, the residuals still seem to have a sine pattern

# Plot mean
# I feel like this increases pretty linearly
meanVec = vecs$mean
t = (1:366)
plot(meanVec) # looks less random than A1
ma.smooth.wge(meanVec,order=75) # looks like an absolute value graph
y = 24.15488 - abs(((366/2)-t))*0.1050146
plot(y)
res = y - meanVec
plot(res) # that does seem to be a bit more random

# Plot without averaging and with averaging
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
freq1 = 1/144
freq2 = 1/(144*366)
dayStart = 1
dayStop = 5
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24+1)
  psi1 = psi1mean
  A1 = abs(A1vec[day])
  b0 = meanVec[day]
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax[t] = st
}
plot(y_gen_vpmax)
numInDay = 144
overlap = 20
lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
}
plot(y_gen_vpmax)

# Plot all 3
freq1 = 1/144
freq2 = 1/(144*366)
psi1 = psi1mean
dayStart = 1
dayStop = 366
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  t = ((day-1)*6*24+1):(day*6*24+1)
  A1 = abs(A1vec[day])
  b0 = meanVec[day]
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax[t] = st
}
numInDay = 144
overlap = 20
lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
}
set.seed(1)
vara = mean(vecs$var)
pvec = vecs$p
phicoeff1vec = vecs$phi1
phicoeff2vec = vecs$phi2
phi1mean = mean(phicoeff1vec[pvec==2]) 
phi2mean = mean(phicoeff2vec[pvec==2])
phi = c(phi1mean,phi2mean)
zt = gen.arma.wge(n = 144*366, phi = phi, theta = 0, vara = vara, plot=FALSE)
numDays = 366
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax[1:(numDays*6*24)]+zt[1:(numDays*6*24)],col = "green")
numDays = 10
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax[1:(numDays*6*24)]+zt[1:(numDays*6*24)],col = "green")
startDay = 366/2
numDays = 10
ymin = min(c(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)])
plot(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "blue")
points(y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]+zt[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "green")

# ACF and spectral density
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6)
ACF2 = acf(y_gen_vpmax+zt, plot = "FALSE")
lines(ACF2$lag, ACF2$acf, lwd = 2, col = "red") # That looks really similar
numDays = 20
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*20)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
SpecDen2 = parzen.wge(y_gen_vpmax[1:(6*24*numDays)]+zt[1:(6*24*numDays)], plot = "FALSE")
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red") # it looks a little higher
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*180):(6*24*numDays+6*24*180)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
SpecDen2 = parzen.wge(y_gen_vpmax[(1+6*24*180):(6*24*numDays+6*24*180)]+zt[(1+6*24*180):(6*24*numDays+6*24*180)], plot = "FALSE")
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red") # better, but looks a little lower

# Try with 10 using b0 from data
numSims = 10
y_gen_vpmax_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_mult) <- paste0("sim", 1:10)
for (seed in 1:numSims){
  print(seed)
  set.seed(seed)
  meanVec = vecs$mean
  t = (1:366)
  y = 24.15488 - abs(((366/2)-t))*0.1050146
  res = y - meanVec
  maxRes = 8
  meanVar = var(res[abs(res)<maxRes]) 
  y = 25.2 - abs(((366/2)-t))*0.096
  b0sim = y + rnorm(n=366,mean=0,sd=sqrt(meanVar))
  dayStart = 1
  dayStop = 366
  y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
  for (day in dayStart:dayStop){
    t = ((day-1)*6*24+1):(day*6*24)
    psi1 = psi1mean
    A1 = abs(A1vec[day])
    b0 = meanVec[day] # use data
    st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
    y_gen_vpmax[t] = st
  }
  numInDay = 144
  overlap = 20
  lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
  firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
  for (ii in 1:overlap){
    beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
    afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
    factor = overlap+1-ii
    y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
    y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
  }
  zt = gen.arma.wge(n = 144*366, phi = phi, theta = 0, vara = vara, sn = 0, plot=FALSE)
  y_gen_vpmax_mult[, seed] = y_gen_vpmax + zt
}
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6)
for (sim in 1:numSims){
  ACF2 = acf(y_gen_vpmax_mult[,sim], plot = "FALSE")
  lines(ACF2$lag, ACF2$acf, lwd = 2, col = "red")
}
numDays = 20
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}
startDay = 180
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}

# Try 10 with b0 from simulation
numSims = 10
y_gen_vpmax_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_mult) <- paste0("sim", 1:10)
for (seed in 1:numSims){
  print(seed)
  set.seed(seed)
  b0sim = y + rnorm(n=366,mean=0,sd=sqrt(meanVar))
  dayStart = 1
  dayStop = 366
  y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
  for (day in dayStart:dayStop){
    t = ((day-1)*6*24+1):(day*6*24)
    psi1 = psi1mean
    A1 = abs(A1vec[day])
    b0 = b0sim[day] # use simulated b0
    st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
    y_gen_vpmax[t] = st
  }
  numInDay = 144
  overlap = 20
  lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
  firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
  for (ii in 1:overlap){
    beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
    afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
    factor = overlap+1-ii
    y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
    y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
  }
  zt = gen.arma.wge(n = 144*366, phi = phi, theta = 0, vara = vara, sn = 0, plot=FALSE)
  y_gen_vpmax_mult[, seed] = y_gen_vpmax + zt
}
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6)
for (sim in 1:numSims){
  ACF2 = acf(y_gen_vpmax_mult[,sim], plot = "FALSE")
  lines(ACF2$lag, ACF2$acf, lwd = 2, col = "red")
}
numDays = 20
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}
startDay = 180
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}

# Things I don't like:
#   - Variance seems to be much larger in the summer months, so average variance makes it high in winter and low in summer
#   - it would be good to not always make the signal match the data so closely
#   - not sure how to generate realistic amplitudes from scratch
#   - there are negative values, where it doesn't even equal 0 ever in the real dataset

################################################################################

# Pull in data
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
library(tswge)

# Create a new vector
f1 = 1/(6*24)
rangeVec = numeric(366)
midpointVec = numeric(366)
meanVec = numeric(366)
A1vec = numeric(366)
phi1vec = numeric(366)
pvec = numeric(366)
phicoeff1vec = numeric(366)
phicoeff2vec = numeric(366)
phicoeff3vec = numeric(366)
phicoeff4vec = numeric(366)
phicoeff5vec = numeric(366)
phicoeff6vec = numeric(366)
phicoeff7vec = numeric(366)
phicoeff8vec = numeric(366)
phicoeff9vec = numeric(366)
phicoeff10vec = numeric(366)
phicoeff11vec = numeric(366)
phicoeff12vec = numeric(366)
ljungVec = numeric(366)
varVec = numeric(366)
dayStart = 1
dayStop = 366
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24)
  vpmax_day = weather$VPmax..mbar.[t]
  rangeVec[day] = (max(vpmax_day) - min(vpmax_day))
  midpoint = (max(vpmax_day) - min(vpmax_day))/2 + min(vpmax_day)
  midpointVec[day] = midpoint
  meanVec[day] = mean(vpmax_day)
  y = vpmax_day - midpoint
  fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1), 
             start = list(A1 = 1, phi1 = 0),
             data = data.frame(t,y))
  A1vec[day] = coef(fit)['A1']
  phi1vec[day] = coef(fit)['phi1']
  res = resid(fit)
  p = aic5.ar.wge(res,p=0:12,method='mle',type='bic')[1,1]
  pvec[day] = p
  est <- tryCatch(est.ar.wge(res,p=p,method='mle'),
                    error=function(e){
                      return(est.ar.wge(res,p=p,method='burg'))
                    }
    )
  phicoeff1vec[day] = est$phi[1]
  phicoeff2vec[day] = est$phi[2]
  phicoeff3vec[day] = est$phi[3]
  phicoeff4vec[day] = est$phi[4]
  phicoeff5vec[day] = est$phi[5]
  phicoeff6vec[day] = est$phi[6]
  phicoeff7vec[day] = est$phi[7]
  phicoeff8vec[day] = est$phi[8]
  phicoeff9vec[day] = est$phi[9]
  phicoeff10vec[day] = est$phi[10]
  phicoeff11vec[day] = est$phi[11]
  phicoeff12vec[day] = est$phi[12]
  varVec[day] = est$avar
  fore = fore.arma.wge(res,phi=est$phi,plot=FALSE)
  ljung <- ljung.wge(fore$resid,K=24)
  ljungVec[day]=ljung$pval
}

# Save off the excel
library(writexl)
vecs <- data.frame(range = rangeVec, midpoint = midpointVec, mean = meanVec, A1 = A1vec, psi1 = phi1vec, p = pvec, 
                   phi1 = phicoeff1vec, phi2 = phicoeff2vec, phi3 = phicoeff3vec, phi4 = phicoeff4vec, 
                   phi5 = phicoeff5vec, phi6 = phicoeff6vec, phi7 = phicoeff7vec, phi8 = phicoeff8vec, 
                   phi9 = phicoeff9vec, phi10 = phicoeff10vec, phi11 = phicoeff11vec, phi12 = phicoeff12vec, 
                   ljung = ljungVec, var = varVec)
write_xlsx(vecs, "C:\\Users\\aabro\\OneDrive\\Desktop\\SMU Program\\Capstone\\Datasets\\vectors3.xlsx")

# Try something similar, but using the mean as the center instead of the midpoint between max and min
f1 = 1/(6*24)
rangeVec = numeric(366)
midpointVec = numeric(366)
meanVec = numeric(366)
A1vec = numeric(366)
phi1vec = numeric(366)
pvec = numeric(366)
phicoeff1vec = numeric(366)
phicoeff2vec = numeric(366)
phicoeff3vec = numeric(366)
phicoeff4vec = numeric(366)
phicoeff5vec = numeric(366)
phicoeff6vec = numeric(366)
phicoeff7vec = numeric(366)
phicoeff8vec = numeric(366)
phicoeff9vec = numeric(366)
phicoeff10vec = numeric(366)
phicoeff11vec = numeric(366)
phicoeff12vec = numeric(366)
ljungVec = numeric(366)
varVec = numeric(366)
dayStart = 1
dayStop = 366
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24)
  vpmax_day = weather$VPmax..mbar.[t]
  rangeVec[day] = (max(vpmax_day) - min(vpmax_day))
  midpoint = (max(vpmax_day) - min(vpmax_day))/2 + min(vpmax_day)
  midpointVec[day] = midpoint
  meanPoint = mean(vpmax_day)
  meanVec[day] = meanPoint
  # y = vpmax_day - midpoint
  y = vpmax_day - meanPoint
  fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1), 
             start = list(A1 = 1, phi1 = 0),
             data = data.frame(t,y))
  A1vec[day] = coef(fit)['A1']
  phi1vec[day] = coef(fit)['phi1']
  res = resid(fit)
  p = aic5.ar.wge(res,p=0:12,method='mle',type='bic')[1,1]
  pvec[day] = p
  est <- tryCatch(est.ar.wge(res,p=p,method='mle'),
                  error=function(e){
                    return(est.ar.wge(res,p=p,method='burg'))
                  }
  )
  phicoeff1vec[day] = est$phi[1]
  phicoeff2vec[day] = est$phi[2]
  phicoeff3vec[day] = est$phi[3]
  phicoeff4vec[day] = est$phi[4]
  phicoeff5vec[day] = est$phi[5]
  phicoeff6vec[day] = est$phi[6]
  phicoeff7vec[day] = est$phi[7]
  phicoeff8vec[day] = est$phi[8]
  phicoeff9vec[day] = est$phi[9]
  phicoeff10vec[day] = est$phi[10]
  phicoeff11vec[day] = est$phi[11]
  phicoeff12vec[day] = est$phi[12]
  varVec[day] = est$avar
  fore = fore.arma.wge(res,phi=est$phi,plot=FALSE)
  ljung <- ljung.wge(fore$resid,K=24)
  ljungVec[day]=ljung$pval
}
vecs <- data.frame(range = rangeVec, midpoint = midpointVec, mean = meanVec, A1 = A1vec, psi1 = phi1vec, p = pvec, 
                   phi1 = phicoeff1vec, phi2 = phicoeff2vec, phi3 = phicoeff3vec, phi4 = phicoeff4vec, 
                   phi5 = phicoeff5vec, phi6 = phicoeff6vec, phi7 = phicoeff7vec, phi8 = phicoeff8vec, 
                   phi9 = phicoeff9vec, phi10 = phicoeff10vec, phi11 = phicoeff11vec, phi12 = phicoeff12vec, 
                   ljung = ljungVec, var = varVec)
write_xlsx(vecs, "C:\\Users\\aabro\\OneDrive\\Desktop\\SMU Program\\Capstone\\Datasets\\vectors4.xlsx")
# It looks like this did the same thing

# Now see if there's a good way to generate a 366 vector of A1 values
plot(A1vec)
absA1vec = abs(A1vec)
plot(absA1vec)
plot(absA1vec,type='l')

# Maybe look at a moving average to smooth it out
ma.smooth.wge(absA1vec,order=5)
ma.smooth.wge(absA1vec,order=15)
ma.smooth.wge(absA1vec,order=25)
ma.smooth.wge(absA1vec,order=75) # this one looks pretty close to half a sine curve
A1SineCurve = ma.smooth.wge(absA1vec,order=75)
A1SineCurve = na.omit(A1SineCurve$smooth)
y = A1SineCurve
subnum = (75-1)/2
t=((1+subnum):(366-subnum))
fit <- nls(y ~ A1 * cos(2 * pi * (1/(366*2)) * t + 0), 
           start = list(A1 = 1),
           data = data.frame(t,y))
summary(fit) 
# This is saying amplitude = 0.56, which isn't right at all
# This makes me question how I was getting amplitude before as well

# You know what, let me just play around with a toy cosine curve
subnum = (75-1)/2
t=((1+subnum):(366-subnum))
y = 6*sin(2 * pi * (1/(366*2)) * t)
plot(y)

# Ohhh, it's sine
y = A1SineCurve
subnum = (75-1)/2
t=((1+subnum):(366-subnum))
fit <- nls(y ~ A1 * sin(2 * pi * (1/(366*2)) * t + 0), 
           start = list(A1 = 1),
           data = data.frame(t,y))
summary(fit) # now amplitude is 6.06861, that makes way more sense

# What would have happened with cosine if I gave it a psi?
y = A1SineCurve
subnum = (75-1)/2
t=((1+subnum):(366-subnum))
fit <- nls(y ~ A1 * cos(2 * pi * (1/(366*2)) * t + psi1), 
           start = list(A1 = 1, psi1 = 0),
           data = data.frame(t,y))
summary(fit) # amplitude is still around 6, that's good, gives me more confidence in the original approach

# Now let's look at the residuals
t = (1:366)
res = absA1vec - 6.06861*sin(2*pi*(1/(366*2))*t)
plot(res) # shoot, these have a pattern but it's a weird pattern

# What happens when we look at mean values over time?
plot(meanVec) # that one looks a little more regular
ma.smooth.wge(meanVec,order=5)
ma.smooth.wge(meanVec,order=15)
ma.smooth.wge(meanVec,order=25)
meanAbsCurve = ma.smooth.wge(meanVec,order=75)
meanAbsCurve = na.omit(meanAbsCurve$smooth)
max(meanAbsCurve) # 24.15488
min(meanAbsCurve) # 8.927757
m = (24.15488-8.927757) / ((366/2)-(1+subnum)) # 0.1050146
b = 8.927757 - m*(1+subnum) # 4.937201
t = (1:366)
y = 24.15488 - abs(((366/2)-t))*0.1050146
plot(y)
res = y - meanVec
plot(res) # that does seem to be a bit more random

# look at variances
plot(varVec)

# For now, I'd still like to try just using the A1 and b0 values gathered from the data
# This as opposed to generating those using a random process
# I feel like I don't know a good way to generate them using a random process, and it's just adding another problem I don't know how to solve
# I'll try to find an overall value for psi, p, and phis, just to make things easier
# I can start off by trying to generate an AR vector using the same variance as well
# This may need to be modified though, but I'd like to at least get something that stitches together properly first
# Then I can plot ACF and spectral density to see how close it is

# let's look at the p values, to see how much p=2 dominates
hist(pvec)
table(pvec) # more than half the values of p are 2
plot(phicoeff1vec[pvec==2]) # range from like 1.2 to 1.8
plot(phicoeff2vec[pvec==2]) # range from like -0.7 to -0.2

# what do some factor tables look like
factor.wge(phi=c(phicoeff1vec[1],phicoeff2vec[1])) # 0.9787, 0.4119
factor.wge(phi=c(phicoeff1vec[2],phicoeff2vec[2])) # 0.9443, 0.6688
factor.wge(phi=c(phicoeff1vec[219],phicoeff2vec[219])) # 0.9402, 0.2296
# These do seem to just be two AR(1) type frequencies, not a periodic frequency

# What does the mean of all these look like?
phi1mean = mean(phicoeff1vec[pvec==2]) # 1.41336
phi2mean = mean(phicoeff2vec[pvec==2]) # -0.4527101
factor.wge(phi=c(phi1mean,phi2mean)) # 0.9228, 0.4906

# What about the mean value for psi?
psi1vec = phi1vec
plot(psi1vec) # some crazy values
psi1vec = psi1vec %% (2*pi)
plot(psi1vec) # still two zones, likely because of the positive and negative A1 values
plot(psi1vec-pi)
psi1vec2 = numeric(366)
for (day in 1:366){
  if (psi1vec[day] < pi){
    psi1vec2[day] = psi1vec[day]
  } else { # psi1vec[day] > pi
    psi1vec2[day] = psi1vec[day] - pi
  }
}
plot(psi1vec2)
psi1mean = mean(psi1vec2) # 1.995956

# How does it look?
plot(cos(2*pi*(1:144)*1/144 + psi1mean))
plot(weather$VPmax..mbar.[1:144])
plot(weather$VPmax..mbar.[(144+1):(144+144)])
plot(weather$VPmax..mbar.[(2*144+1):(2*144+144)]) # wow, this day doesn't really seem sinusoidal at all
plot(weather$VPmax..mbar.[1:(144*10)])

# How about plotting 10 days where it uses the mean and amplitude values
freq1 = 1/144
freq2 = 1/(144*366)
A2 = -7.932237
psi2 = -0.379505
dayStart = 1
dayStop = 5
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24+1)
  psi1 = psi1mean
  A1 = abs(A1vec[day])
  b0 = meanVec[day]
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1) + A2 * cos(2 * pi * freq2 * t + psi2)
  y_gen_vpmax[t] = st
}
plot(y_gen_vpmax)
plot(weather$VPmax..mbar.[1:(144*dayStop)])

# What about the A2 part
plot(A2 * cos(2 * pi * freq2 * (1:(dayStop*144)) + psi2)) # hm, very negative, maybe just remove this part for now

# Try with no A2
freq1 = 1/144
freq2 = 1/(144*366)
dayStart = 1
dayStop = 5
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24+1)
  psi1 = psi1mean
  A1 = abs(A1vec[day])
  b0 = meanVec[day]
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax[t] = st
}
plot(y_gen_vpmax) # at least they're postive now
plot(weather$VPmax..mbar.[1:(144*dayStop)])

# Now let's try to connect the end of one day to the beginning of another day
numInDay = 144
overlap = 5
for (ii in 0:overlap){
  lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
  firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
  beforeDay = seq(numInDay-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])/2*factor + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])/2*factor + y_gen_vpmax[afterDay]*ii)/(overlap+1)
}
plot(y_gen_vpmax)

# That looks better, but some of the further away ones don't look great
# Let's try with a larger overlap
freq1 = 1/144
freq2 = 1/(144*366)
dayStart = 1
dayStop = 5
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24+1)
  psi1 = psi1mean
  A1 = abs(A1vec[day])
  b0 = meanVec[day]
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax[t] = st
}
plot(y_gen_vpmax)
numInDay = 144
overlap = 20
lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
  if (ii < 5){
    print(ii)
    print(y_gen_vpmax)
  }
}
plot(y_gen_vpmax)
# 10 looks better
# 20 looks really good

# What's going on with the 289 point?
ii = 1
overlap = 20
factor = overlap+1-ii
(y_gen_vpmax[288]+y_gen_vpmax[289])*factor/2 # 137.3795
y_gen_vpmax[289]*ii # 7.28271
((y_gen_vpmax[288]+y_gen_vpmax[289])*factor/2 + y_gen_vpmax[289]*ii)/(overlap+1)
afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
afterDay
((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)

ii = 2
factor = overlap+1-ii
afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
afterDay
((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)

# This looks pretty good, what does the whole year look like?
freq1 = 1/144
freq2 = 1/(144*366)
dayStart = 1
dayStop = 366
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24+1)
  psi1 = psi1mean
  A1 = abs(A1vec[day])
  b0 = meanVec[day]
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax[t] = st
}
plot(y_gen_vpmax)
numInDay = 144
overlap = 20
lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
}
plot(y_gen_vpmax)

# Plots on same graph
numDays = 366
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
# These look very similar, not quite as big of a range though

# Let's see how to add the variance part
phi1mean = mean(phicoeff1vec[pvec==2]) 
phi2mean = mean(phicoeff2vec[pvec==2])
phi = c(phi1mean,phi2mean)
vara = mean(vecs$var)
set.seed(1)
zt = gen.arma.wge(n = 144*366, phi = phi, theta = 0, vara = vara, sn = 0, plot=FALSE)
plot(zt + y_gen_vpmax)

# Plot all 3
numDays = 366
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax[1:(numDays*6*24)]+zt[1:(numDays*6*24)],col = "green")

# Only the first few days
numDays = 10
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax[1:(numDays*6*24)]+zt[1:(numDays*6*24)],col = "green")

# Some days in the middle of the year 
startDay = 366/2
numDays = 10
ymin = min(c(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)])
plot(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "blue")
points(y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]+zt[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "green")

# Let's try comparing ACFs
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6)
ACF2 = acf(y_gen_vpmax+zt, plot = "FALSE")
lines(ACF2$lag, ACF2$acf, lwd = 2, col = "red")
# That looks really similar

# How about some spectral densities
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*10)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
SpecDen2 = parzen.wge(y_gen_vpmax[1:(6*24*10)]+zt[1:(6*24*10)], plot = "FALSE")
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
# It looks higher

SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*180):(6*24*10+6*24*180)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
SpecDen2 = parzen.wge(y_gen_vpmax[(1+6*24*180):(6*24*10+6*24*180)]+zt[(1+6*24*180):(6*24*10+6*24*180)], plot = "FALSE")
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
# And that one looks a little lower

# What happens if I try to generate b0 instead of using it from the data
meanVec = vecs$mean
t = (1:366)
y = 24.15488 - abs(((366/2)-t))*0.1050146
res = y - meanVec
meanVar = var(res) # 17.71848
b0sim = y + rnorm(n=366,mean=0,sd=sqrt(meanVar))
plot(b0sim) # hm, that goes negative

# what if I try to remove large residuals?
maxRes = 5
meanVar = var(res[abs(res)<maxRes]) # 5.601831
b0sim = y + rnorm(n=366,mean=0,sd=sqrt(meanVar))
plot(b0sim)

plot(b0sim,col='blue')
points(meanVec,col='red')
legend('bottom',legend = c("Sim", "Real"), col = c("blue", "red"),pch=1)
# It doesn't seem to get that low to 0, maybe make max higher and slope less

# Try with different slope
meanVec = vecs$mean
t = (1:366)
y = 24.15488 - abs(((366/2)-t))*0.1050146
res = y - meanVec
maxRes = 8
meanVar = var(res[abs(res)<maxRes]) 
y = 25.2 - abs(((366/2)-t))*0.096
b0sim = y + rnorm(n=366,mean=0,sd=sqrt(meanVar))
plot(b0sim)
ymin = min(c(b0sim,meanVec))
ymax = max(c(b0sim,meanVec))
plot(b0sim,col='blue',ylim=c(ymin,ymax))
points(meanVec,col='red')
legend('bottom',legend = c("Sim", "Real"), col = c("blue", "red"),pch=1)

# That seems alright, what happens when we use that instead of the real b0?
set.seed(1)
meanVec = vecs$mean
t = (1:366)
y = 24.15488 - abs(((366/2)-t))*0.1050146
res = y - meanVec
maxRes = 8
meanVar = var(res[abs(res)<maxRes]) 
y = 25.2 - abs(((366/2)-t))*0.096
b0sim = y + rnorm(n=366,mean=0,sd=sqrt(meanVar))
psi1mean = mean(psi1vec2) 
freq1 = 1/144
freq2 = 1/(144*366)
dayStart = 1
dayStop = 5
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24+1)
  psi1 = psi1mean
  A1 = abs(A1vec[day])
  b0 = b0sim[day]
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax[t] = st
}
plot(y_gen_vpmax)
numInDay = 144
overlap = 20
lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
}
plot(y_gen_vpmax)

# That seems fine, what about for the full 366 days?
dayStart = 1
dayStop = 366
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24+1)
  psi1 = psi1mean
  A1 = abs(A1vec[day])
  b0 = b0sim[day]
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax[t] = st
}
plot(y_gen_vpmax)
numInDay = 144
overlap = 20
lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
}
plot(y_gen_vpmax)
set.seed(1)
zt = gen.arma.wge(n = 144*366, phi = phi, theta = 0, vara = vara, sn = 0, plot=FALSE)
numDays = 366
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax[1:(numDays*6*24)]+zt[1:(numDays*6*24)],col = "green")
numDays = 10
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax[1:(numDays*6*24)]+zt[1:(numDays*6*24)],col = "green")
startDay = 366/2
numDays = 10
ymin = min(c(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)])
plot(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "blue")
points(y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]+zt[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "green")

# That dips into the negatives
meanVec = vecs$mean
t = (1:366)
y = 24.15488 - abs(((366/2)-t))*0.1050146
res = y - meanVec
maxRes = 6
meanVar = var(res[abs(res)<maxRes]) 
y = 26 - abs(((366/2)-t))*0.094
set.seed(1)
b0sim = y + rnorm(n=366,mean=0,sd=sqrt(meanVar))
ymin = min(c(b0sim,meanVec))
ymax = max(c(b0sim,meanVec))
plot(b0sim,col='blue',ylim=c(ymin,ymax))
points(meanVec,col='red')
legend('bottom',legend = c("Sim", "Real"), col = c("blue", "red"),pch=1)
dayStart = 1
dayStop = 366
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24+1)
  psi1 = psi1mean
  A1 = abs(A1vec[day])
  b0 = b0sim[day]
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax[t] = st
}
numInDay = 144
overlap = 20
lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
}
set.seed(1)
zt = gen.arma.wge(n = 144*366, phi = phi, theta = 0, vara = vara, sn = 0, plot=FALSE)
numDays = 366
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax[1:(numDays*6*24)]+zt[1:(numDays*6*24)],col = "green")
legend('topleft',legend = c("VPMax", "Sim - no noise", "Sim - with noise"), col = c("red", "blue", "green"),pch=1)
numDays = 10
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax[1:(numDays*6*24)]+zt[1:(numDays*6*24)],col = "green")
legend('topleft',legend = c("VPMax", "Sim - no noise", "Sim - with noise"), col = c("red", "blue", "green"),pch=1)
startDay = 366/2
numDays = 10
ymin = min(c(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)])
plot(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "blue")
points(y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]+zt[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "green")
legend('bottomleft',legend = c("VPMax", "Sim - no noise", "Sim - with noise"), col = c("red", "blue", "green"),pch=1)

# Hm, I wonder if we can compare b0 to a1
plot(b0sim - abs(A1vec))
plot(b0sim,col='red')
points(abs(A1vec),col='blue')
plot(b0sim/abs(A1vec))
min(b0sim/abs(A1vec))

meanVec = vecs$mean
t = (1:366)
y = 24.15488 - abs(((366/2)-t))*0.1050146
res = y - meanVec
maxRes = 6
meanVar = var(res[abs(res)<maxRes]) 
y = 27 - abs(((366/2)-t))*0.094
set.seed(1)
b0sim = y + rnorm(n=366,mean=0,sd=sqrt(meanVar))
ampFac = 1.2
b0sim = pmax(b0sim,ampFac*abs(A1vec))
ymin = min(c(b0sim,meanVec))
ymax = max(c(b0sim,meanVec))
plot(b0sim,col='blue',ylim=c(ymin,ymax))
points(meanVec,col='red')
legend('bottom',legend = c("Sim", "Real"), col = c("blue", "red"),pch=1)
dayStart = 1
dayStop = 366
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  t = ((day-1)*6*24+1):(day*6*24+1)
  psi1 = psi1mean
  A1 = abs(A1vec[day])
  b0 = b0sim[day]
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax[t] = st
}
numInDay = 144
overlap = 20
lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
}
set.seed(1)
zt = gen.arma.wge(n = 144*366, phi = phi, theta = 0, vara = vara, sn = 0, plot=FALSE)
numDays = 366
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax[1:(numDays*6*24)]+zt[1:(numDays*6*24)],col = "green")
legend('topleft',legend = c("VPMax", "Sim - no noise", "Sim - with noise"), col = c("red", "blue", "green"),pch=1)
numDays = 10
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax[1:(numDays*6*24)]+zt[1:(numDays*6*24)],col = "green")
legend('topleft',legend = c("VPMax", "Sim - no noise", "Sim - with noise"), col = c("red", "blue", "green"),pch=1)
startDay = 366/2
numDays = 10
ymin = min(c(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)])
plot(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "blue")
points(y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]+zt[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "green")
legend('bottomleft',legend = c("VPMax", "Sim - no noise", "Sim - with noise"), col = c("red", "blue", "green"),pch=1)
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6)
ACF2 = acf(y_gen_vpmax+zt, plot = "FALSE")
lines(ACF2$lag, ACF2$acf, lwd = 2, col = "red") # That looks really similar
numDays = 50
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*20)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
SpecDen2 = parzen.wge(y_gen_vpmax[1:(6*24*numDays)]+zt[1:(6*24*numDays)], plot = "FALSE")
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red") # it looks a little higher
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*180):(6*24*numDays+6*24*180)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
SpecDen2 = parzen.wge(y_gen_vpmax[(1+6*24*180):(6*24*numDays+6*24*180)]+zt[(1+6*24*180):(6*24*numDays+6*24*180)], plot = "FALSE")
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red") # better, but looks a little lower

# This doesn't look that great to be honest, let's see what happens when we do 10 of them to get an ACF or parzen window
numSims = 10
y_gen_vpmax_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_mult) <- paste0("sim", 1:10)
for (seed in 1:numSims){
  print(seed)
  set.seed(seed)
  b0sim = y + rnorm(n=366,mean=0,sd=sqrt(meanVar))
  dayStart = 1
  dayStop = 366
  y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
  for (day in dayStart:dayStop){
    t = ((day-1)*6*24+1):(day*6*24)
    psi1 = psi1mean
    A1 = abs(A1vec[day])
    b0 = b0sim[day] # use simulated b0
    st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
    y_gen_vpmax[t] = st
  }
  numInDay = 144
  overlap = 20
  lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
  firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
  for (ii in 1:overlap){
    beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
    afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
    factor = overlap+1-ii
    y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
    y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
  }
  zt = gen.arma.wge(n = 144*366, phi = phi, theta = 0, vara = vara, sn = 0, plot=FALSE)
  y_gen_vpmax_mult[, seed] = y_gen_vpmax + zt
}
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6)
for (sim in 1:numSims){
  ACF2 = acf(y_gen_vpmax_mult[,sim], plot = "FALSE")
  lines(ACF2$lag, ACF2$acf, lwd = 2, col = "red")
}
numDays = 20
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}
startDay = 180
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}

# How about with the b0 from the data?
numSims = 10
y_gen_vpmax_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_mult) <- paste0("sim", 1:10)
for (seed in 1:numSims){
  print(seed)
  set.seed(seed)
  b0sim = y + rnorm(n=366,mean=0,sd=sqrt(meanVar))
  dayStart = 1
  dayStop = 366
  y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
  for (day in dayStart:dayStop){
    t = ((day-1)*6*24+1):(day*6*24)
    psi1 = psi1mean
    A1 = abs(A1vec[day])
    b0 = meanVec[day] # use data
    st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
    y_gen_vpmax[t] = st
  }
  numInDay = 144
  overlap = 20
  lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
  firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
  for (ii in 1:overlap){
    beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
    afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
    factor = overlap+1-ii
    y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
    y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
  }
  zt = gen.arma.wge(n = 144*366, phi = phi, theta = 0, vara = vara, sn = 0, plot=FALSE)
  y_gen_vpmax_mult[, seed] = y_gen_vpmax + zt
}
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6)
for (sim in 1:numSims){
  ACF2 = acf(y_gen_vpmax_mult[,sim], plot = "FALSE")
  lines(ACF2$lag, ACF2$acf, lwd = 2, col = "red")
}
numDays = 20
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}
startDay = 180
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}

# Demo
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
library(tswge)
f1 = 1/(6*24)
rangeVec = numeric(366)
b0vec = numeric(366)
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
dayStop = 21
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24)
  vpmax_day = weather$VPmax..mbar.[t]
  rangeVec[day] = (max(vpmax_day) - min(vpmax_day))
  b0 = (max(vpmax_day) - min(vpmax_day))/2 + min(vpmax_day)
  b0vec[day] = b0
  y = vpmax_day - b0
  fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1), 
             start = list(A1 = 1, phi1 = 0),
             data = data.frame(t,y))
  A1vec[day] = coef(fit)['A1']
  phi1vec[day] = coef(fit)['phi1']
  res = resid(fit)
  p = aic5.ar.wge(res,p=0:12,method='mle',type='bic')[1,1]
  pvec[day] = p
  capture.output(
    est <- tryCatch(est.ar.wge(res,p=p,method='mle'),
                    error=function(e){
                      return(est.ar.wge(res,p=p,method='burg'))
                    }
    )
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
  capture.output(ljung <- ljung.wge(fore$resid,K=24))
  ljungVec[day]=ljung$pval
}

# Save off the excel
library(writexl)
vecs <- data.frame(range = rangeVec, b0 = b0vec, A1 = A1vec, psi1 = phi1vec, p = pvec, 
                   phi1 = phicoeff1vec, phi2 = phicoeff2vec, phi3 = phicoeff3vec, phi4 = phicoeff4vec, 
                   phi5 = phicoeff5vec, phi6 = phicoeff6vec, phi7 = phicoeff7vec, phi8 = phicoeff8vec, 
                   phi9 = phicoeff9vec, phi10 = phicoeff10vec, phi11 = phicoeff11vec, phi12 = phicoeff12vec, 
                   ljung = ljungVec, var = varVec)
write_xlsx(vecs, "C:\\Users\\aabro\\OneDrive\\Desktop\\SMU Program\\Capstone\\Datasets\\vectors2.xlsx")

# Plot A1s, ps, and vars
plot(abs(A1vec))
plot(A1vec)
plot(phi1vec) # psi
plot(b0vec)
plot(pvec)
plot(varVec)

# First plot
y_gen_vpmax = numeric(6*24*366)
n = 6*24
b0 = median(weather$VPmax..mbar.)
sn = 2
A2 = -7.932237
psi2 = -0.379505
freq1 = 1/(6*24)
freq2 = 1/(6*24*366)
diff = 0
dayStart = 1
dayStop = 4
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24)
  est = est.ar.wge(res,p=pvec[day],method='mle')
  phi = est$phi
  # phi = c(phicoeff1vec[day],phicoeff2vec[day])
  vara = varVec[day]
  psi1 = phi1vec[day]
  A1 = A1vec[day]
  b0 = b0vec[day]
  st = b0 + diff + A1 * cos(2 * pi * freq1 * t + psi1) + A2 * cos(2 * pi * freq2 * t + psi2)
  zt = gen.arma.wge(n = n, phi = phi, theta = 0, vara = vara, sn = sn, plot=FALSE)
  xt = st + zt
  y_gen_vpmax[t] = xt
  diff = xt[6*24] - xt[1] # reset the diff after each day, to make sure that there is some continuity
}

numDays = dayStop
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")

# Things attempted
# 1. Constant psi1 values
# 2. Trying to set a different diff value
# 3. Changing the mu value
# 4. gen.arma.ada where you can set the starting value


################################################################################

# Files
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/vectors.csv"
prev_vecs <- read.csv(file_path, header = TRUE)

# I kinda want to see what the different phi values are
library(tswge)
day = 1
t = ((day-1)*6*24+1):(day*6*24)
f1 = 1/(6*24)
vpmax_1 = weather$VPmax..mbar.[t]
b0 = (max(vpmax_1) - min(vpmax_1))/2 + min(vpmax_1)
y = vpmax_1 - b0
fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1), 
           start = list(A1 = 1, phi1 = 0),
           data = data.frame(t,y))
res = resid(fit)
p = aic5.ar.wge(res,p=0:12,method='mle',type='bic')[1,1]
est <- est.ar.wge(res,p=p,method='mle')
# It's interesting, this says p = 2, but before it said p = 3 always
# This makes me think I did something wrong

day = 2
t = ((day-1)*6*24+1):(day*6*24)
vpmax_1 = weather$VPmax..mbar.[t]
b0 = (max(vpmax_1) - min(vpmax_1))/2 + min(vpmax_1)
y = vpmax_1 - b0
fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1), 
           start = list(A1 = 1, phi1 = 0),
           data = data.frame(t,y))
res = resid(fit)
p = aic5.ar.wge(res,p=0:12,method='mle',type='bic')[1,1]
est <- est.ar.wge(res,p=p,method='mle')
# The abs recip from the factor tables are 0.9787 and 0.4119 vs 0.9443 and 0.6688
# These are pretty close, but not exact

# I'll try this again, but also save off the different phi coeffs
f1 = 1/(6*24)
rangeVec = numeric(366)
b0vec = numeric(366)
A1vec = numeric(366)
phi1vec = numeric(366)
pvec = numeric(366)
phicoeff1vec = numeric(366)
phicoeff2vec = numeric(366)
phicoeff3vec = numeric(366)
phicoeff4vec = numeric(366)
ljungVec = numeric(366)
varVec = numeric(366)
dayStart = 21
dayStop = 366
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24)
  vpmax_day = weather$VPmax..mbar.[t]
  rangeVec[day] = (max(vpmax_day) - min(vpmax_day))
  b0 = (max(vpmax_day) - min(vpmax_day))/2 + min(vpmax_day)
  b0vec[day] = b0
  y = vpmax_day - b0
  fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1), 
             start = list(A1 = 1, phi1 = 0),
             data = data.frame(t,y))
  A1vec[day] = coef(fit)['A1']
  phi1vec[day] = coef(fit)['phi1']
  res = resid(fit)
  capture.output(p = aic5.ar.wge(res,p=0:12,method='mle',type='bic')[1,1])
  pvec[day] = p
  capture.output(
    est <- tryCatch(est.ar.wge(res,p=p,method='mle'),
                    error=function(e){
                      return(est.ar.wge(res,p=p,method='burg'))
                    }
    )
  )
  phicoeff1vec[day] = est$phi[1]
  phicoeff2vec[day] = est$phi[2]
  phicoeff3vec[day] = est$phi[3]
  
  phicoeff4vec[day] = est$phi[4]
  varVec[day] = est$avar
  fore = fore.arma.wge(res,phi=est$phi,plot=FALSE)
  capture.output(ljung <- ljung.wge(fore$resid,K=24))
  ljungVec[day]=ljung$pval
}

# Try making a generation using different A1, phi1, var, and phis
n = 6*24
b0 = median(weather$VPmax..mbar.)
sn = 1
A2 = -7.932237
psi2 = -0.379505
freq1 = 1/(6*24)
freq2 = 1/(6*24*366)

day = 1
t = ((day-1)*6*24+1):(day*6*24)
phi = c(phicoeff1vec[day],phicoeff2vec[day])
vara = varVec[day]
psi1 = phi1vec[day]
A1 = A1vec[day]
diff = 0
st = b0 + diff + A1 * cos(2 * pi * freq1 * t + psi1) + A2 * cos(2 * pi * freq2 * t + psi2)
zt = gen.arma.wge(n = n, phi = phi, theta = 0, vara = vara, sn = sn)
xt = st + zt
diff = xt[day*6*24] - xt[(day-1)*6*24+1]

ymin = min(weather$VPmax..mbar.[t],xt)
ymax = max(weather$VPmax..mbar.[t],xt)
plot(t,weather$VPmax..mbar.[t],col = "red",ylim=c(ymin,ymax))
points(t,xt,col = "blue")

day = 2
t = ((day-1)*6*24+1):(day*6*24)
phi = c(phicoeff1vec[day],phicoeff2vec[day])
vara = varVec[day]
psi1 = phi1vec[day]
A1 = A1vec[day]
diff = 0
st = b0 + diff + A1 * cos(2 * pi * freq1 * t + psi1) + A2 * cos(2 * pi * freq2 * t + psi2)
zt = gen.arma.wge(n = n, phi = phi, theta = 0, vara = vara, sn = sn)
xt = st + zt
diff = xt[day*6*24] - xt[(day-1)*6*24+1]

ymin = min(weather$VPmax..mbar.[t],xt)
ymax = max(weather$VPmax..mbar.[t],xt)
plot(t,weather$VPmax..mbar.[t],col = "red",ylim=c(ymin,ymax))
points(t,xt,col = "blue")

# That seems alright, let's try the whole thing
y_gen_vpmax = numeric(6*24*366)
n = 6*24
b0 = median(weather$VPmax..mbar.)
sn = 2
A2 = -7.932237
psi2 = -0.379505
freq1 = 1/(6*24)
freq2 = 1/(6*24*366)
diff = 0
dayStart = 1
dayStop = 366
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24)
  phi = c(phicoeff1vec[day],phicoeff2vec[day])
  vara = varVec[day]
  psi1 = phi1vec[day]
  A1 = A1vec[day]
  st = b0 + diff + A1 * cos(2 * pi * freq1 * t + psi1) + A2 * cos(2 * pi * freq2 * t + psi2)
  capture.output(zt = gen.arma.wge(n = n, phi = phi, theta = 0, vara = vara, sn = sn))
  xt = st + zt
  y_gen_vpmax[t] = xt
  diff = xt[6*24] - xt[1] # reset the diff after each day, to make sure that there is some continuity
}

plot(y_gen_vpmax)

numDays = 6
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")

# Hm, there is kind of a weird gap from day to day.
# I wonder if this is because of me using different phi1 values?
# Maybe I should just choose an average value and stick with that.
test = phi1vec %% (2*pi)
test[test>pi] = test[test>pi]-pi
plot(test)
psi1 = mean(test)

y_gen_vpmax = numeric(6*24*366)
n = 6*24
b0 = median(weather$VPmax..mbar.)
sn = 2
A2 = -7.932237
psi1s = phi1vec %% (2*pi)
psi1s[psi1s>pi] = psi1s[psi1s>pi]-pi
psi1 = mean(psi1s)
psi2 = -0.379505
freq1 = 1/(6*24)
freq2 = 1/(6*24*366)
diff = 0
dayStart = 1
dayStop = 366
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24)
  phi = c(phicoeff1vec[day],phicoeff2vec[day])
  vara = varVec[day]
  A1 = abs(A1vec[day])
  st = b0 + diff + A1 * cos(2 * pi * freq1 * t + psi1) + A2 * cos(2 * pi * freq2 * t + psi2)
  capture.output(zt = gen.arma.wge(n = n, phi = phi, theta = 0, vara = vara, sn = sn))
  xt = st + zt
  y_gen_vpmax[t] = xt
  diff = xt[6*24] - xt[1] # reset the diff after each day, to make sure that there is some continuity
}

numDays = 10
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")

# It still isn't working.  It just doesn't look like it's stitching correctly together.
# Maybe I should take the difference of what it would be the next day.
y_gen_vpmax = numeric(6*24*366)
n = 6*24
b0 = median(weather$VPmax..mbar.)
sn = 1
A2 = -7.932237
psi2 = -0.379505
freq1 = 1/(6*24)
freq2 = 1/(6*24*366)
diff = 0
mu = 0
dayStart = 1
dayStop = 5
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24+1)
  phi = c(phicoeff1vec[day],phicoeff2vec[day])
  vara = varVec[day]
  psi1 = phi1vec[day]
  A1 = A1vec[day]
  st = b0 + diff + A1 * cos(2 * pi * freq1 * t + psi1) + A2 * cos(2 * pi * freq2 * t + psi2)
  zt = gen.arma.wge(n = n+1, phi = phi, theta = 0, vara = vara, mu=mu, sn = sn, plot=FALSE) # It looks like zt can start at a random value that isn't related to the previous value
  xt = st[1:n] + zt[1:n]
  y_gen_vpmax[t] = xt
  if (day < dayStop){
    nextt = (day*6*24+1):((day+1)*6*24+1)
    nexts = b0 + diff + A1vec[day+1] * cos(2 * pi * freq1 * nextt + phi1vec[day+1]) + A2 * cos(2 * pi * freq2 * nextt + psi2)
    mu = nexts[1] - st[n+1] # make sure the AR part doesn't stray too far from expect
    nextz = gen.arma.wge(n = n+1, phi = c(phicoeff1vec[day+1],phicoeff2vec[day+1]), theta = 0, vara = varVec[day+1], mu=mu, sn = sn, plot=FALSE) # include the AR value, so it doesn't jump to much
    nextVal = nexts[1] + nextz[1]
    nextValWant = st[n+1] + zt[n+1] # connect the next zt vec, to the previous one
    diff = nextValWant - nextVal # reset the diff after each day, to make sure that there is some continuity
  }
}

numDays = dayStop
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
# This can make it run away increasing, apparently

# Make a function that can take in a starting value
gen.arma.ada = function(n, phi, theta, vara, starting_value, mu, sn)
{
  set.seed(sn)
  Xt = numeric(n)
  et = numeric(n)
  Xt[1] = starting_value
  et[1] = rnorm(1,0,sqrt(vara))
  p = length(phi)
  q = length(theta)
  for (time in 2: n)
  {
    et[time] = rnorm(1,0,sqrt(vara))
    Xpart = 0
    Apart = 0
    p_to_use = max(min(p,time-p+1),1) # trying to make sure it doesn't go too far back
    q_to_use = max(min(q,time-q+1,1),1) # making sure that it doesn't go below 1 also
    for (ii in 1:p_to_use)
    {
      Xpart = Xpart + phi[ii]*Xt[time-ii] 
    }
    for (ii in 1:q_to_use)
    {
      Apart = Apart + theta[ii]*et[time-ii]
    }
    Xt[time] = Xpart + mu*(1-sum(phi)) +  Apart + et[time]
  }
  return(Xt)
}

# I'd like to stitch it together, but also not have it run away
y_gen_vpmax = numeric(6*24*366)
n = 6*24
b0 = median(weather$VPmax..mbar.)
sn = 1
A2 = -7.932237
psi2 = -0.379505
freq1 = 1/(6*24)
freq2 = 1/(6*24*366)
diff = 0
mu = 0
starting_value = gen.arma.wge(1,phi=c(phicoeff1vec[1],phicoeff2vec[1]),theta=0,vara=varVec[1],mu=0,sn=sn)[1]
dayStart = 1
dayStop = 366
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24+1)
  phi = c(phicoeff1vec[day],phicoeff2vec[day])
  vara = varVec[day]
  psi1 = phi1vec[day]
  A1 = A1vec[day]
  st = b0 + diff + A1 * cos(2 * pi * freq1 * t + psi1) + A2 * cos(2 * pi * freq2 * t + psi2)
  zt = gen.arma.ada(n = n+1, phi = phi, theta = 0, vara = vara, starting_value = starting_value, mu=mu, sn = sn) # It looks like zt can start at a random value that isn't related to the previous value
  xt = st[1:n] + zt[1:n]
  y_gen_vpmax[t] = xt
  if (day < dayStop){
    nextt = (day*6*24+1):((day+1)*6*24+1)
    nexts = b0 + A1vec[day+1] * cos(2 * pi * freq1 * nextt + phi1vec[day+1]) + A2 * cos(2 * pi * freq2 * nextt + psi2)
    mu = nexts[1] - st[n+1] # make sure the AR part doesn't stray too far from expect
    starting_value = zt[n+1]
    nextz = gen.arma.ada(n = n+1, phi = c(phicoeff1vec[day+1],phicoeff2vec[day+1]), theta = 0, vara = varVec[day+1], starting_value = starting_value, mu=mu, sn = sn)
    nextVal = nexts[1] + nextz[1]
    nextValWant = st[n+1] + zt[n+1] # connect the next zt vec, to the previous one
    nextVal = nexts[1] + nextz[1]
    nextValWant = st[n+1] + zt[n+1] # connect the next zt vec, to the previous one
    diff = nextValWant - nextVal
    }
}

numDays = dayStop
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
# Apparently that doesn't really work

# Maybe just do 1 phi and vara throughout
y_gen_vpmax = numeric(6*24*366)
n = 6*24
b0 = median(weather$VPmax..mbar.)
sn = 1
A2 = -7.932237
psi2 = -0.379505
freq1 = 1/(6*24)
freq2 = 1/(6*24*366)
diff = 0
vara = median(varVec)
phi = c(median(phicoeff1vec),median(phicoeff2vec))
zt = gen.arma.wge(n = n*366, phi = phi, theta = 0, vara = vara, sn = sn, plot=FALSE)
dayStart = 1
dayStop = 366
for (day in dayStart:dayStop){
  if ((day %% 20) == 0){
    print(day)
  }
  t = ((day-1)*6*24+1):(day*6*24+1)
  psi1 = phi1vec[day]
  A1 = A1vec[day]
  st = b0 + diff + A1 * cos(2 * pi * freq1 * t + psi1) + A2 * cos(2 * pi * freq2 * t + psi2)
  xt = st[1:n] + zt[t]
  y_gen_vpmax[t] = xt
  if (day < dayStop){
    nextt = (day*6*24+1):((day+1)*6*24+1)
    nexts = b0 + diff + A1vec[day+1] * cos(2 * pi * freq1 * nextt + phi1vec[day+1]) + A2 * cos(2 * pi * freq2 * nextt + psi2)
    diff = st[n+1] - nexts[1] # reset the diff after each day, to make sure that there is some continuity
  }
}

numDays = dayStop
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")

acf(y_gen_vpmax[1:(numDays*6*24)])
acf(weather$VPmax..mbar.)

dayStart = 150
dayStop = 160
parzen.wge(y_gen_vpmax[((dayStart-1)*6*24+1):(dayStop*6*24+1)],trunc=100)
parzen.wge(weather$VPmax..mbar.[((dayStart-1)*6*24+1):(dayStop*6*24+1)],trunc=100)
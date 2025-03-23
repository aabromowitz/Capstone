# Get data
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)

# Let's try getting the amplitude for the first day
vpmax_1 = weather$VPmax..mbar.[1:(6*24)]
plot(vpmax_1)

# Figure out freq and amplitude of the data
f1 = 1/(6*24)
y = vpmax_1 - mean(vpmax_1)
t = 1:(6*24)
fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1), 
           start = list(A1 = 1, phi1 = 0),
           data = data.frame(t,y))

# Plot comparison
plot(t, y, pch = 16, col = "blue", xlab = "t", ylab = "y", main = "Original Data and Fitted Curve")
y_fit <- predict(fit)
lines(t, y_fit, col = "red", lwd = 2)
legend("topright", legend = c("Original Data", "Fitted Curve"), 
       col = c("blue", "red"), pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2))

# plot residuals
resid = y-y_fit
plot(resid) # Hm, these are periodic

# Try with 2
f1 = 1/(6*24)
y = vpmax_1 - mean(vpmax_1)
t = 1:(6*24)
fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f1 * t + phi2), 
           start = list(A1 = 1, A2 = 1, phi1 = 0, phi2 = 0),
           data = data.frame(t,y))
plot(resid(fit))

# Look at the freqs
library(tswge)
p1 = parzen.wge(y, trunc = 100)
p = parzen.wge(resid, trunc = 100) # Maybe try f1 and f1*2?

# two diff freqs
f1 = 1/(6*24)
f2 = 2*f1
y = vpmax_1 - mean(vpmax_1)
t = 1:(6*24)
fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2), 
           start = list(A1 = 1, A2 = 1, phi1 = 0, phi2 = 0),
           data = data.frame(t,y))
plot(resid(fit)) # It still looks like it has a frequency, sigh

# What about 1 freq and AR
f1 = 1/(6*24)
t = 1:(6*24)
vpmax_1 = weather$VPmax..mbar.[1:(6*24)]
y = vpmax_1 - mean(vpmax_1)
fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1), 
           start = list(A1 = 1, phi1 = 0),
           data = data.frame(t,y))
resid = resid(fit)
plot(resid)
y_sin_pred = predict(fit)
plot(y_sin_pred + b0)
p = aic5.ar.wge(resid,p=0:12,method='mle',type='bic')[1,1]
est=est.ar.wge(resid,p=p,method='mle')
fore = fore.arma.wge(resid,phi=est$phi,lastn=FALSE)
plot(t, vpmax_1, pch = 16, cex = 1.5, col = "red", xlab = "t", ylab = "VPMax", main = "VPMax First Day")
points(y_sin_pred + mean(vpmax_1), col = "blue", pch = 16, cex = 1.5)
points(y_sin_pred + resid + mean(vpmax_1) - fore$resid, col = "green", pch = 16, cex = 1.5) 
plot(fore$resid) # it looks kinda like white noise to me
ljung.wge(fore$resid,K=24) # default is K=24, want p>0.05 to show that it is white noise
ljung.wge(fore$resid,K=48)
# Neither are considered white noise

# What if I tried messing with the A value and the b0 value, since it seems like using a larger value might help.
b0 = (max(vpmax_1) - min(vpmax_1))/2 + min(vpmax_1) # use the mid-point instead of the mean value
y = vpmax_1 - b0
fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1), 
           start = list(A1 = 1, phi1 = 0),
           data = data.frame(t,y))
resid = resid(fit)
plot(resid) # This looks the exact same somehow
y_sin_pred = predict(fit)
plot(y_sin_pred + b0)
p = aic5.ar.wge(resid,p=0:12,method='mle',type='bic')[1,1]
est=est.ar.wge(resid,p=p,method='mle')
fore = fore.arma.wge(resid,phi=est$phi,lastn=FALSE)
plot(t, vpmax_1, pch = 16, cex = 1.5, col = "red", xlab = "t", ylab = "VPMax", main = "VPMax First Day")
points(y_sin_pred + b0, col = "blue", pch = 16, cex = 1.5) # This is shifted up, but the others are the same
points(y_sin_pred + resid + b0 - fore$resid, col = "green", pch = 16, cex = 1.5) 
plot(fore$resid) # it looks kinda like white noise to me
ljung.wge(fore$resid,K=24) # default is K=24, want p>0.05 to show that it is white noise
ljung.wge(fore$resid,K=48)

# I'll try with 2 frequencies, but I'm not really expecting it to be any different
t = 1:(6*24)
f1 = 1/(6*24)
f2 = 2*f1
vpmax_1 = weather$VPmax..mbar.[1:(6*24)]
b0 = mean(vpmax_1)
y = vpmax_1 - b0
fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2), 
           start = list(A1 = 1, A2 = 1, phi1 = 0, phi2 = 0),
           data = data.frame(t,y))
resid = resid(fit)
plot(resid)
y_sin_pred = predict(fit)
plot(y_sin_pred + b0)
p = aic5.ar.wge(resid,p=0:12,method='mle',type='bic')[1,1]
est=est.ar.wge(resid,p=p,method='mle')
fore = fore.arma.wge(resid,phi=est$phi,lastn=FALSE)
plot(t, vpmax_1, pch = 16, cex = 1.5, col = "red", xlab = "t", ylab = "VPMax", main = "VPMax First Day")
points(y_sin_pred + b0, col = "blue", pch = 16, cex = 1.5)
points(y_sin_pred + resid + b0 - fore$resid, col = "green", pch = 16, cex = 1.5) 
plot(fore$resid) # this looks different
ljung.wge(fore$resid,K=24) # default is K=24, want p>0.05 to show that it is white noise
ljung.wge(fore$resid,K=48) # wow, this is actually higher than 0.05

# Let's try with either 3 or 4 frequencies, since there are 3.5 peaks in the resid plot
parz = parzen.wge(resid,trunc = 100)
t = 1:(6*24)
f1 = 1/(6*24)
f2 = 2*f1
f3 = 2*f2
vpmax_1 = weather$VPmax..mbar.[1:(6*24)]
b0 = mean(vpmax_1)
y = vpmax_1 - b0
fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2) + A3 * cos(2 * pi * f3 * t + phi3), 
           start = list(A1 = 1, A2 = 1, A3 = 1, phi1 = 0, phi2 = 0, phi3 = 0),
           data = data.frame(t,y))
summary(fit)
resid = resid(fit)
plot(resid)
y_sin_pred = predict(fit)
plot(y_sin_pred + b0)
p = aic5.ar.wge(resid,p=0:12,method='mle',type='bic')[1,1]
est=est.ar.wge(resid,p=p,method='mle')
fore = fore.arma.wge(resid,phi=est$phi,lastn=FALSE)
plot(t, vpmax_1, pch = 16, cex = 1.5, col = "red", xlab = "t", ylab = "VPMax", main = "VPMax First Day")
points(y_sin_pred + b0, col = "blue", pch = 16, cex = 1.5)
points(y_sin_pred + resid + b0 - fore$resid, col = "green", pch = 16, cex = 1.5) 
plot(fore$resid) 
ljung.wge(fore$resid,K=24) # default is K=24, want p>0.05 to show that it is white noise
ljung.wge(fore$resid,K=48) # 0.01 and 0.07

# Try with 4
t = 1:(6*24)
f1 = 1/(6*24)
f2 = 2*f1
f3 = 2*f2
f4 = 2*f3
vpmax_1 = weather$VPmax..mbar.[1:(6*24)]
b0 = mean(vpmax_1)
y = vpmax_1 - b0
fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2) + 
             A3 * cos(2 * pi * f3 * t + phi3) + A4 * cos(2 * pi * f4 * t + phi4), 
           start = list(A1 = 1, A2 = 1, A3 = 1, A4 = 1, phi1 = 0, phi2 = 0, phi3 = 0, phi4 = 0),
           data = data.frame(t,y))
summary(fit)
resid = resid(fit)
plot(resid)
y_sin_pred = predict(fit)
plot(y_sin_pred + b0)
p = aic5.ar.wge(resid,p=0:12,method='mle',type='bic')[1,1]
est=est.ar.wge(resid,p=p,method='mle')
fore = fore.arma.wge(resid,phi=est$phi,lastn=FALSE)
plot(t, vpmax_1, pch = 16, cex = 1.5, col = "red", xlab = "t", ylab = "VPMax", main = "VPMax First Day")
points(y_sin_pred + b0, col = "blue", pch = 16, cex = 1.5)
points(y_sin_pred + resid + b0 - fore$resid, col = "green", pch = 16, cex = 1.5) 
plot(fore$resid) 
ljung.wge(fore$resid,K=24) # default is K=24, want p>0.05 to show that it is white noise
ljung.wge(fore$resid,K=48) # now p is 0.04, so 1 or 2 is honestly probably the best

# Is there any evidence for a 2nd frequency?
y = weather$VPmax..mbar. - mean(weather$VPmax..mbar.)
fft_result <- fft(y)
power_spectrum <- Mod(fft_result)^2
freqs <- (0:(length(y)-1)) / length(y)  # Normalize frequencies
dominant_indices <- order(power_spectrum, decreasing = TRUE)[1:5]  # Get top 5 frequencies
freqs[dominant_indices[1]] # 1.897389e-05 / 52704.01 hours / 366 days
freqs[dominant_indices[2]] # 0.999981
freqs[dominant_indices[3]] # 0.006944444 / 144 hours / 1 day
freqs[dominant_indices[4]] # 0.9930556
freqs[dominant_indices[5]] # 3.794778e-05
# So none of the top 5 mention the 2nd frequency

# Let's make vectors with the different values
t = 1:(6*24)
f1 = 1/(6*24)
rangeVec = numeric(366)
b0vec = numeric(366)
A1vec = numeric(366)
phi1vec = numeric(366)
pvec = numeric(366)
ljungVec = numeric(366)
varVec = numeric(366)
dayStart = 1
for (day in dayStart:366){
  if ((day %% 20) == 0){
    print(day)
  }
  vpmax_day = weather$VPmax..mbar.[((day-1)*6*24+1):(day*6*24)]
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
  varVec[day] = est$avar
  fore = fore.arma.wge(res,phi=est$phi,plot=FALSE)
  capture.output(ljung <- ljung.wge(fore$resid,K=24))
  ljungVec[day]=ljung$pval
}

plot(A1vec)
plot(abs(A1vec))
plot(rangeVec)
plot(b0vec)
plot(phi1vec)
plot(phi1vec %% (2*pi))
plot(varVec)
which.max(rangeVec)

# Save off the excel
library(writexl)
vecs <- data.frame(range = rangeVec, b0 = b0vec, A1 = A1vec, phi1 = phi1vec, p = pvec, ljung = ljungVec, var = varVec)
write_xlsx(vecs, "C:\\Users\\aabro\\OneDrive\\Desktop\\SMU Program\\Capstone\\Datasets\\vectors.xlsx")

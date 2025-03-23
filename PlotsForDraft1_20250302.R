# Pull in data
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather.csv"
weather <- read.csv(file_path, header = TRUE)

# Plot temp for full year
t = 1:(6*24*366)
t = t/(6*24)
len = length(weather$date)
plot(x=t[1:len],y=weather$T..degC.[1:len], type = "p", col = "blue", pch = 16, cex = 0.5,
     xlab = "Days since 1/1/2020", ylab = "Air Temperatuare (C)", 
     main = "Air Temperatuare - Full Year")

# Plot a week
len = 6*24*7
plot(x=t[1:len],y=weather$T..degC.[1:len], type = "p", col = "blue", pch = 16, cex = 0.8,
     xlab = "Days since 1/1/2020", ylab = "Air Temperatuare (C)", 
     main = "Air Temperatuare - Week")

# Plot as a line
len = 6*24*7
plot(x=t[1:len],y=weather$T..degC.[1:len], type = "l", col = "blue", pch = 16, cex = 1,
     xlab = "Days since 1/1/2020", ylab = "Air Temperatuare (C)", 
     main = "Air Temperatuare - Week")

# plot the linear simulation for a full year
weather$simLin = (5.509857-1.689e+00) + 0.981397*weather$Tdew..degC. + 6.500e-03*weather$PAR...mol.m..s.
len = length(weather$date)
plot(x=t[1:len],y=weather$T..degC.[1:len], type = "p", col = "blue", pch = 16, cex = 0.5,
     xlab = "Days since 1/1/2020", ylab = "Air Temperatuare (C)", 
     main = "Air Temperatuare - Full Year")
points(x=t[1:len],y=weather$simLin[1:len], col = "red", pch = 16, cex = 0.5)
legend("topleft", legend = c("Air Temp", "Dew Temp + PAR"), 
       col = c("blue", "red"), lwd = 2)

# Plot a week
len = 6*24*7
plot(x=t[1:len],y=weather$T..degC.[1:len], type = "p", col = "blue", pch = 16, cex = 0.5,
     xlab = "Days since 1/1/2020", ylab = "Air Temperatuare (C)", 
     main = "Air Temperatuare - Week")
points(x=t[1:len],y=weather$simLin[1:len], col = "red", pch = 16, cex = 0.5)
legend("topleft", legend = c("Air Temp", "Dew Temp + PAR"), 
       col = c("blue", "red"), lwd = 2)
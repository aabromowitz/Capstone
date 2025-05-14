# Pull in data
library(tswge)
library(vars)
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)

# Try the VARSelect
weather_sub = weather[,c('T..degC.','Tdew..degC.','VPmax..mbar.')]
# VARselect(weather_sub,lag.max=18,type="const",season=144,exogen=NULL)
# AIC(n)  HQ(n)  SC(n) FPE(n) 
# 17     13     12     17

# Set the VAR model
fit = VAR(weather_sub,p=12,type='const',season=144)
# summary(fit) # bad idea, this wrecked my computer

# Make the forecast
l = length(weather_sub$T..degC.)
h=720
preds=predict(fit,n.ahead=h)
ase = mean((weather_sub$T..degC.[(l-h+1):l]-preds$fcst$T..degC.[,1])^2)
ase # 29.4991
ymax = max(weather_sub$T..degC.[(l-h+1):l],preds$fcst$T..degC.[1:h,1])
ymin = min(weather_sub$T..degC.[(l-h+1):l],preds$fcst$T..degC.[1:h,1])
plot(seq(l-h+1,l,1),weather_sub$T..degC.[(l-h+1):l],type="b",pch=15,col='blue',ylim=c(ymin,ymax))
points(seq(l-h+1,l,1),preds$fcst$T..degC.[1:h,1],type="b",pch=15,col='red',ylim=c(ymin,ymax))
legend('topleft',legend = c("Real Data", "Predicted"), 
       col = c("blue", "red"),pch=15)
# This might actually be worse
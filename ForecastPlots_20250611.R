# Try looking at plots again

# Common variables
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Sims/temp_1_to_50_varTemp001.csv"
sims <- read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Sims/Signal.csv"
sig_tot = read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Forecasts/ARIMA_Forecasts_With_Trimmed_1_to_50_varTemp001.csv"
arima_fores = read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Forecasts/1_to_50_varTemp001_TFT_combined.csv"
tft_fores = read.csv(file_path, header = TRUE)
len_total = length(sims[,1])
h = 720
vlines = c(-720,0, 96, 192, 336)
b = h*2
sig = sig_tot[(len_total-h-b+1):len_total,2]
# interesting_plots = c(1,2,3,6,9,10,11,12,16,19,20,23,24,32,38,43,44,45,47,49)
interesting_plots = seq(1,50,1)
num_plots = length(interesting_plots)

# Loop through plots
for (iPlot in 1:num_plots){
  sim_num = interesting_plots[iPlot]
  arima_fore = arima_fores[,sim_num]
  tft_fore = tft_fores[,sim_num]
  sim = sims[(len_total-h-b+1):len_total,sim_num]
  y_min = min(c(sim,arima_fore,tft_fore,sig))
  y_max = max(c(sim,arima_fore,tft_fore,sig))
  out.path <- file.path(dest.dir, paste0("forecast_sim_", sim_num, ".png"))
  png(filename = out.path, width = 800, height = 600)
  plot(seq(1,h+b,1),sim, col='blue',ylim=c(y_min,y_max),
       main=paste('Forecast Comparison - Simulation', sim_num), xlab = 'Days',ylab = 'Air Temperature (C)',xaxt='n')
  for (iLine in 1:length(vlines)){
    vline = vlines[iLine] + b
    if (vline > 0){
      abline(v=vline,col="black",lty=2,lwd=0.5)
    }
  }
  if(-720+b > 0){
    abline(h=sim[-720+b],col="black",lty=2,lwd=0.5)
  }
  abline(h=sim[b],col="black",lty=2,lwd=0.5)
  points(seq(1,h+b,1),sig, col='green')
  points(seq(1+b,h+b,1),arima_fore,col = "red")
  points(seq(1+b,h+b,1),tft_fore,col = "magenta")
  legend('top',legend = c("Simulation", "Signal", "ARIMA Forecast", "ARIMA Forecast"), 
         col = c("blue", "green", "red", "magenta"),pch=1)
  ticks <- axTicks(1)
  axis(1, at = ticks, labels = round((ticks + len_total-h-b) / 144,0))
  dev.off()
}

# This is supposed to save off all the plots
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
plots.png.paths <- list.files(plots.dir.path, pattern="\\.png$", full.names = TRUE)
dest.dir <- "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Plots/RunForJun7"
file.copy(from=plots.png.paths, to=dest.dir)

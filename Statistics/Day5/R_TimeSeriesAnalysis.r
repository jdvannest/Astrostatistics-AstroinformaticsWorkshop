# Read the dataset
download.file(url = "https://git.psu.edu/stat/astro/datasets/-/raw/master/GX.dat", destfile = "GX.dat")
GX.dat <- scan("GX.dat") 
str(GX.dat)
summary(GX.dat)

# Convert to a class 'ts' (time series) R object

GX.time <- seq(from=0, to=512, length.out=length(GX.dat))
GX.ts <-  ts(GX.dat, GX.time) ; GX.ts.offset <- ts(GX.dat-30, GX.time)
str(GX.ts)

plot.ts(GX.ts, ylab='GX 5-1 counts', xlab='Time (x 1/128 sec)', 
   cex.lab=1.3, cex.axis=1.3, lwd=0.5)

# Compare histogram of counts to normal distribution

options(jupyter.plot_scale=1)
options(repr.plot.width = 7, repr.plot.height = 5)

hist(GX.dat, breaks=100, xlim=c(40,100), ylim=c(0,3500), xlab='GX 5-1 counts',
   font=2, font.lab=2, main='')
curve(dnorm(x,mean=mean(GX.dat), sd=sqrt(mean(GX.dat)))*65536, lwd=3, add=T)
sd(GX.dat) / sqrt(mean(GX.dat))  # counts are 1.24 x overdispersed compared to Gaussian noise

plot.ts(GX.ts[1:6000], ylab='GX 5-1 counts', xlab='Time (x 1/128 sec)', 
   cex.lab=1.3, cex.axis=1.3)  # Close-up view of 10% of the data

plot(GX.time,GX.dat, ylim=c(-10,115), xlab='Time (sec)', ylab='GX 5-1 counts',
   cex.lab=1.3, cex.axis=1.3, type='n')  # set up plot window but don't show any data
lines(ksmooth(GX.time, GX.dat+30, 'normal', bandwidth=7), lwd=2) 
text(450, 110, 'Normal kernel')  # Gaussian kernel density estimator with 7 bin FWHM bandwidth
lines(stats::filter(GX.ts, sides=2, rep(1,7)/7), lwd=2) 
text(450, 85, 'Moving average') # Moving average smoother with 7 bin bandwidth
lines(kernapply(GX.ts.offset, kernel('modified.daniell', 7)), lwd=2) 
text(450, 50, 'Modified Daniell')  # Moving average smoother with 1/2-weight at the end values of the span
lines(supsmu(GX.time, GX.dat-60, span=0.01), lwd=2) 
text(400, 20, "Friedman's super-smoother") # A smoother with adaptive bandwidth from Friedman (1984)
lines(lowess(GX.time, GX.dat-80, 0.02), lwd=2) 
text(400, 0, 'LOWESS local regression') # Cleveland's (1979) robust local polynomial smoother

par(mfrow=c(1,2))
acf(GX.dat)
pacf(GX.dat)
par(mfrow=c(1,1))

# Spectral periodogram

par(mfrow=c(3,1))  ;  par(mar=c(5,4,1,2))
spec.pgram(GX.ts, log='no', main='')
spec.pgram(GX.ts, spans=50, log='no', main='', sub='')
spec.pgram(GX.ts, spans=200, taper=0.15, log='no', main='', sub='')

# Autoregressive modeling: AR

ARmod <- ar(GX.ts, method='ols') 
print(ARmod)

ARspec <- spec.ar(GX.ts, plot=F)
GXspec <- spec.pgram(GX.ts, span=101, main='', sub='', lwd=2)
lines(ARspec$freq, ARspec$spec, col='green', lwd=2)
text(0.4,450, cex=2, paste0('AR(', ARmod$order, ')'))

acf(na.omit(ARmod$resid), ylim=c(-0.03,0.03))
Box.test(ARmod$resid, type='Ljung')  # test for autocorrelation in AR residuals
install.packages("goftest")
library(goftest)
ad.test(na.omit(ARmod$resid))  # test for normality in AR residuals

# Autoregressive modeling: ARIMA and ARFIMA
install.packages("forecast")
library(forecast)

ARIMA_fit <- auto.arima(GX.ts)
summary(ARIMA_fit)

par(mfrow=c(2,1))  ;  par(mar=c(5,4,1,2))
arima_fit <- auto.arima(GX.ts, stepwise=FALSE, approximation=FALSE, max.p=3, max.q=3, max.d=1)
acf(arima_fit$residuals, ylim=c(-0.1,0.1))
ad.test(na.omit(ARmod$resid))  # test for normality in AR residuals

arfima_fit <- arfima(GX.ts)
summary(arfima_fit)
acf(arfima_fit$residuals)
Box.test(arfima_fit$residuals, type='Ljung-Box')
sd(arfima_fit$x) ; sd(arfima_fit$residuals)


# Ingest and plot light curves for three Kepler stars
download.file(url = "https://git.psu.edu/stat/astro/datasets/-/raw/master/Kepler1.dat", destfile = "Kepler1.dat")
download.file(url = "https://git.psu.edu/stat/astro/datasets/-/raw/master/Kepler2.dat", destfile = "Kepler2.dat")
Kepler1 <- read.table('Kepler1.dat')[[1]] # KIC 007596240
Kepler2 <- read.table('Kepler2.dat')[[1]] # KIC 007609553

length(which(is.na(Kepler1))) / length(Kepler1) # 16% NAs
length(which(is.na(Kepler2))) / length(Kepler2) # 29% NAs

par(mfrow=c(2,1)) ; par(mar=c(5,4,1,2))
plot(Kepler1, type='l', xlab='Time')
plot(Kepler2, type='l', xlab='Time')

# Properties of the Kepler 1 lightcurve (KIC 007596240)

cat(' Kepler 1: Median flux = ', median(Kepler1, na.rm=TRUE), 
	' with InterQuartile Range = ', IQR(Kepler1, na.rm=TRUE))

par(mfrow=c(2,1)) ; par(mar=c(5,4,1,2))
acf(Kepler1, na.action=na.pass, ylim=c(-0.05, 0.2),
xlab='Kepler 1 lag', main='', ci.col='black')
hist(Kepler1, freq=FALSE, breaks=200, main='', xlim=c(-20,20),
xlab='Kepler 1 values')

library(MASS) # Comparison with normal distribution [computationally intensive]
Kep1_mn <- fitdistr(Kepler1[!is.na(Kepler1)], 'normal')[[1]][1]
Kep1_sd <- fitdistr(Kepler1[!is.na(Kepler1)], 'normal')[[1]][2]
curve(dnorm(x, Kep1_mn, Kep1_sd), -20, 20, add=TRUE)
install.packages("nortest")
library(nortest) ; ad.test(Kepler1)

cor.test(1:length(Kepler1), Kepler1, method='kendall')$p.value # Test for trend

Kep1_arima <- arima(Kepler1, order=c(2,1,2)) # Autoregressive model
print(Kep1_arima)# look at a summary of the 'arima' function output
str(Kep1_arima)  # look in detail at the 'arima' function output

acf(Kep1_arima$residuals, na.action=na.pass, ylim=c(-0.05, 0.2), xlab='Kepler 1 lag', ci.col='black')  # shows excellent improvement in autocorrelation ... 
IQR(Kepler1, na.rm=TRUE) ; IQR(Kep1_arima$residuals, na.rm=TRUE)  # ... but little improvement in noise

# Diagnostics tests on Kepler 1 ARIMA residuals

# install.packages('nortest')
# library(nortest)
ad.test(Kep1_arima$residuals) # test for normality

install.packages('imputeTS')
library(imputeTS)
arima_resids <- na_kalman(Kep1_arima$residuals)
x <- 1:length(arima_resids) ; y <- rnorm(length(arima_resids))
lmobject <- lm(y ~ x)
lmobject$residuals <- arima_resids
Box.test(arima_resids, type='Ljung') # Ljung-Box test shows no autocorrelation in the ARIMA residuals

Kepler2_impute <- na_kalman(Kepler2) # impute NAs with Kalman smoother
par(mfrow=c(2,1)) ; par(mar=c(5,4,1,2))
plot(Kepler2, type='l', xlim=c(60000,70000), xlab='Time', ylab='Kepler 2')
plot(Kepler2_impute, type='l', xlim=c(60000,70000), xlab='Time',
ylab='Kepler 2 imputed')

# Three periodograms: Fourier, Lomb-Scargle, epoch folding

par(mfrow=c(3,1)) ; par(mar=c(5,4,1,2))
spec.pgram(Kepler2_impute, xlim=c(0,0.005), spans=5, taper=0.0,
main='', ylab='Fourier', sub='')

install.packages('lomb') 
library(lomb)
lsp(Kepler2, from=0.00, to=0.005, ylab='Lomb-Scargle', main='')

install.packages('RobPer')
library(RobPer)
Kepler2_temp <- cbind(1:length(Kepler2), Kepler2)
Kepler2_irreg <- Kepler2_temp[!is.na(Kepler2_temp[,2]),]
PDM_per <- seq(from=1, to=10001, length.out=1000)
Kepler_PDM <- RobPer(Kepler2_irreg, weighting=FALSE, regression='L2',
model='step', steps=10, periods=PDM_per)
PDM_freq <- 1/PDM_per
plot(PDM_freq, Kepler_PDM, type='l', xlab='frequency', ylab='PDM',
xlim=c(0,0.005))

# Discrete wavelet transform for Kepler 2

install.packages('waveslim') 
library(waveslim)
Kepler2_wavdat <- Kepler2_impute[1:2^16]
Kepler2_dwt <- dwt(Kepler2_wavdat,n.levels=10)

par(mfrow=c(3,1)) ; par(mar=c(5,4,1,2))
plot.ts(up.sample(Kepler2_dwt[[5]],2^5), type='h', ylab='') 
abline(h=0) 
plot.ts(up.sample(Kepler2_dwt[[7]],2^7), type='h', ylab='', lwd=2)
abline(h=0) 
plot.ts(up.sample(Kepler2_dwt[[9]],2^{9}), type='h', ylab='', lwd=2) 
abline(h=0)

# Wavelet denoising

install.packages('wavethresh') 
library(wavethresh)
Kepler2_wd <- wd(Kepler2_wavdat)
Kepler2_wdth <- threshold(Kepler2_wd, policy='universal')
Kepler2_thresh <- wr(Kepler2_wdth)

par(mfrow=c(2,1)) ; par(mar=c(4,4,1,1))
plot(Kepler2_wavdat, type='l', ylab='Kepler 2', xlim=c(5000,8000), ylim=c(-150,200))
plot(Kepler2_thresh, type='l', ylab='Kepler 2', xlim=c(5000,8000), ylim=c(-150,200))

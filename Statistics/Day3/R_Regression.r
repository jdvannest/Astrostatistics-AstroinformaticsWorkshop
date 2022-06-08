# I. Construct large and small samples of 77K SDSS quasars

install.packages('astrodatR')  # package with datasets from the MSMA textbook
library(astrodatR)
data(SDSS_QSO)     # brings tabular data into an R data.frame

# Basic information about a data.frame
dim(SDSS_QSO)      # number of rows and columns
names(SDSS_QSO)    # names of columns
summary(SDSS_QSO)  # 5-number summary of each column

# Remove some bad photometry

qso <- SDSS_QSO[-which(SDSS_QSO[,3] == 0 | SDSS_QSO[,9] == 0),]
qso <- qso[-which(qso[,4] == 9.999 | qso[,12] == 9.999),]
qso[(qso[,4]<0.02),4 ] <- 0.02        # set threshold on magnitude errors
dim(qso) ; summary(qso)
attach(qso)

# Plot dataset of SDSS quasar i vs. u magnitudes showing 
# some of the heteroscedastic measurement errors, with contours for dense regions

options(jupyter.plot_scale=1)   # omit this if running on local computer
options(repr.plot.width = 7, repr.plot.height = 5)

plot(i_mag, u_mag, pch=20, cex=0.1, col='#00000040', xlim=c(16,21), 
   ylim=c(16.5,23.5), xlab="SDSS i (mag)", ylab="SDSS u (mag)")
for(i in 50:150) {
   lines(c(i_mag[i],i_mag[i]),c((u_mag[i]+sig_u_mag[i]),
      (u_mag[i]-sig_u_mag[i])), lwd=2, col='purple2')
   lines(c((i_mag[i]+sig_i_mag[i]),(i_mag[i]-sig_i_mag[i])),
      c(u_mag[i],u_mag[i]), lwd=2, col='purple2')   }

library(KernSmooth)
smqso <- bkde2D(cbind(i_mag, u_mag), bandwidth=c(0.05, 0.05), gridsize=c(400,400))
contour(smqso$x1, smqso$x2, smqso$fhat, add=T, col='gold', nlevels=9)

# II. Ordinary least squares fit

fit_ols <- lm(u_mag~i_mag)
summary(fit_ols) 
confint(fit_ols, level=0.90)          # 3 sigma equivalent for Gaussian distribution

plot(i_mag, u_mag, pch=20, cex=0.1, col='#00000040', xlim=c(16,21), 
   ylim=c(16.5,23.5), xlab="SDSS i (mag)", ylab="SDSS u (mag)")
abline(fit_ols$coef, lty=1, lwd=2)     # solid black line

# III. Weighted least squares fit

fit_wt <- lm(u_mag~i_mag, x=T, weights=1/(sig_u_mag*sig_u_mag))
summary(fit_wt)

plot(i_mag, u_mag, pch=20, cex=0.1, col='#00000040', xlim=c(16,21), 
   ylim=c(16.5,23.5), xlab="SDSS i (mag)", ylab="SDSS u (mag)")
abline(fit_wt$coef,lty=2,lwd=2, col='green')    # dashed green line

# Diagnostic plots involving regression residuals help identify outliers

par(mfrow=c(2,2))		
plot(fit_wt, which=c(2:5), caption='', sub.caption='' ,pch=20, cex=0.3, 
   cex.lab=1.3, cex.axis=1.3)

# Goodness-of-fit for a 2D regression

install.packages('cramer')
library(cramer)
cramer.test(fit_wt$model$u_mag[1:1000], fit_wt$fitted.values[1:1000], replicates=100)

plot(ecdf(fit_wt$model$u_mag[1:1000]), cex=0, col='royalblue', lwd=2, xlab='SDSS u (mag)', ylab='c.d.f.', main='Observed and fitted u magnitude')
plot(ecdf(fit_wt$fitted.values[1:1000]), cex=0, col='magenta', lwd=2, add=TRUE)
text(22, 0.5, 'Observed', cex=1.5, col='royalblue')
text(19, 0.8, 'Model', cex=1.5, col='magenta')

# Robust M-estimator

library(MASS)
fit_M <- rlm(u_mag~i_mag, method='M')	# robust fit with Huber's psi functon
summary(fit_M)  
aM <- fit_M$coef[[1]] ; bM <- fit_M$coef[[2]]

plot(i_mag, u_mag, pch=20, cex=0.1, col='#00000040', xlim=c(16,21), 
   ylim=c(16.5,23.5), xlab="SDSS i (mag)", ylab="SDSS u (mag)")
lines(c(16,21), c(aM+bM*16, aM+bM*21), lty=3, lwd=3, col='royalblue3') # dotted royal blue line

fit_Mwt <- rlm(u_mag~i_mag, method='M', weights=1/(sig_u_mag*sig_u_mag), 
   wt.method='inv.var')   # robust fit with measurement error weighting 
summary(fit_Mwt)  
aMwt <- fit_Mwt$coef[[1]] ; bMwt <- fit_Mwt$coef[[2]]

lines(c(16,21), c(aMwt+bMwt*16, aMwt+bMwt*21), lty=3, lwd=5, col='orange')
text(19.5, 17, 'u = 0.13 + 1.02*i', cex=1.3, col='orange')

abline(fit_wt$coef,lty=2,lwd=2, col='green')    # dashed green line (weighted LS fit)

# Unpack galaxy profile data

data(ell_gal_profile)
summary(ell_gal_profile)
NGC4472 <- ell_gal_profile[ell_gal_profile[,1] == 'NGC.4472',2:3]
NGC4472
radius <- NGC4472[,1]
surf_mag <- NGC4472[,2]

NGC4472.fit <-  nls(surf_mag ~ -2.5*log10(I.e * 10^(-(0.868*n-0.142)*
   ((radius/r.e)^{1/n}-1))) + 26, data=list(NGC4472), start=list(I.e=20.,
   r.e=120.,n=4.), model=T, trace=T)
summary(NGC4472.fit)
logLik(NGC4472.fit)


# Plot NGC 4472 data and best-fit model

par(mai=c(1,1,0.8,0.44))   # improve left-hand margin
plot(NGC4472.fit$model$radius, NGC4472.fit$model$surf_mag, pch=20, 
   xlab="r  (arcsec)", ylab=expression(mu ~~ (mag/sq.arcsec)), ylim=c(16,28), 
   cex.lab=1.5, cex.axis=1.5)
lines(NGC4472.fit$model$radius, fitted(NGC4472.fit))


# Details information about the nls fit

formula(NGC4472.fit)    # formula used
coef(NGC4472.fit)       # best-fit parameters
vcov(NGC4472.fit)       # best-fit parameter covariance matrix
logLik(NGC4472.fit)     # log-likelihood of best fit

confint(NGC4472.fit)    # 95% confidence intervals
profile(NGC4472.fit)    # profiles (cuts) around the best fit (unfortunately the console output is verbose here)

fitted(NGC4472.fit)     # fitted values
residuals(NGC4472.fit)  # residuals from the fitted values

# Residual plot

plot(NGC4472.fit$model$radius,residuals(NGC4472.fit), xlab="r (arcsec)", 
   ylab="Residuals", pch=20, cex.lab=1.5, cex.axis=1.5)
lines(supsmu(NGC4472.fit$model$radius, residuals(NGC4472.fit), span=0.05), 
   lwd=2)

# Test for normality (OK) and autocorrelation (not OK) of residuals
# For linear models, also use the Durbin-Watson test in CRAN packages lmtest and car

qqnorm(residuals(NGC4472.fit) / summary(NGC4472.fit)$sigma) 
abline(a=0,b=1)
shapiro.test(residuals(NGC4472.fit) / summary(NGC4472.fit)$sigma) 
acf(residuals(NGC4472.fit))

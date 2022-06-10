# Setup

install.packages('fpc', repos='https://cloud.r-project.org')           # Flexible Procedures for Clustering
library(fpc)
install.packages('class', repos='https://cloud.r-project.org')         # Functions for Classification
library(class)
install.packages('randomForest', repos='https://cloud.r-project.org')  # Random Forests
library(randomForest)
install.packages('e1071', repos='https://cloud.r-project.org')         # Support Vector Machine
library(e1071)

download.file(url = "https://git.psu.edu/stat/astro/datasets/-/raw/master/COMBO17_lowz.dat", destfile = "COMBO17_lowz.dat")
COMBO_loz <- read.table("COMBO17_lowz.dat", header = TRUE) 
str(COMBO_loz)
dim(COMBO_loz)
names(COMBO_loz)
names(COMBO_loz) <- c('MB', 'M280-MB')  ; names(COMBO_loz)

plot(COMBO_loz, pch=20, cex=0.5, xlim=c(-22,-7), ylim=c(-2,2.5), 
   xlab=expression(M[B]~~(mag)), ylab=expression(M[280] - M[B]~~(mag)), 
   main='')

# Two-dimensional kernel-density estimator

library(MASS)
COMBO_loz_sm <- kde2d(COMBO_loz[,1], COMBO_loz[,2], h=c(1.6,0.4), 
   lims = c(-22,-7,-2,2.5), n=500)
image(COMBO_loz_sm, col=grey(13:0/15), xlab=expression(M[B]~~(mag)), 
   ylab=expression(M[280] - M[B]~~(mag)), xlim=c(-22,-7), ylim=c(-2,2.5), 
   xaxp=c(-20,-10,2))
text(-16.5, -1, "Blue cloud", col='darkblue', pos=4, cex=0.8)
text(-17,-0.7, 'Green valley', col='darkgreen', pos=4, cex=0.8)
text(-13, -0.2, 'Red sequence', col='red', pos=4, cex=0.8)
text(-18.5, 1.7, 'Bright cluster galaxies', col='deeppink3', pos=4, cex=0.8)
dev.copy2pdf(file='COMBO17_CMD.pdf')

# Standardize variables and compute a Euclidean distance matrix

Mag_std <- scale(COMBO_loz[,1]) 
Color_std <- scale(COMBO_loz[,2])
COMBO_std <- cbind(Mag_std,Color_std)
COMBO_dist <- dist(COMBO_std)

# Hierarchical clustering

COMBO_hc <- hclust(COMBO_dist, method='ward.D')
plot(COMBO_hc, label=F)

# Cutting the tree at k=5 clusters

COMBO_hc5a <- rect.hclust(COMBO_hc, k=5, border='red') 
COMBO_hc5b <- cutree(COMBO_hc, k=5) 

# Show this hierarchical clustering solution

plot(COMBO_loz, pch=(COMBO_hc5b+18), cex=0.7, xlab=expression(M[B]~~(mag)), 
   ylab=expression(M[280] - M[B]~~(mag)), main='')

# Density-based clustering algorithm

# install.packages('fpc')
# library(fpc)
COMBO_dbs <-  dbscan(COMBO_std, eps=0.1, MinPts=5, method='raw')
print.dbscan(COMBO_dbs) ; COMBO_dbs$cluster

plot(COMBO_loz[COMBO_dbs$cluster==0,], pch=20, cex=0.7, xlab='M_B (mag)',
   ylab='M_280 - M_B (mag)')
points(COMBO_loz[COMBO_dbs$cluster==2,], pch=2, cex=1.0, col='blue2')
points(COMBO_loz[COMBO_dbs$cluster==1 | COMBO_dbs$cluster==3,], pch=1, cex=1.0, col='orangered3')

install.packages("dbscan")                                           
library(dbscan)

COMBO_hdbs <-  hdbscan(COMBO_std, minPts=25)
print(COMBO_hdbs) 
print(COMBO_hdbs$cluster_scores)
cat('COMBO-17 galaxy ', which.max(COMBO_hdbs$outlier_scores), ' is the most extreme outlier')

plot(COMBO_hdbs)
plot(COMBO_hdbs$hc)

# SDSS point sources test dataset, N=17,000 (mag<21, point sources, hi-qual)

download.file(url = "https://git.psu.edu/stat/astro/datasets/-/raw/master/SDSS_test.csv", destfile = "SDSS_test.csv")
SDSS <- read.csv('SDSS_test.csv', header = TRUE)
dim(SDSS) ; summary(SDSS)
SDSS_test <- data.frame(  # create data frame of SDSS colors
    u_g = SDSS$u_mag - SDSS$g_mag,
    g_r = SDSS$g_mag - SDSS$r_mag,
    r_i = SDSS$r_mag - SDSS$i_mag,
    i_z = SDSS$i_mag - SDSS$z_mag
    )
names(SDSS_test) <- c('u_g', 'g_r', 'r_i', 'i_z')
str(SDSS_test) 

par(mfrow=c(1,3))  # plot test set in black because labels are not yet established
plot(SDSS_test[,1], SDSS_test[,2], xlim=c(-0.7,3), ylim=c(-0.6,1.6),pch=20, col='#11111120',
	cex=0.4, cex.lab=1.5, cex.axis=1.5, main='Test dataset', xlab='u-g (mag)', ylab='g-r (mag)') 
plot(SDSS_test[,2], SDSS_test[,3], xlim=c(-0.6,1.6), ylim=c(-0.6,1.6), pch=20, col='#11111120',
	cex=0.4, cex.lab=1.5, cex.axis=1.5, main='', xlab='g-r (mag)', ylab='r-i (mag)') 
plot(SDSS_test[,3], SDSS_test[,4], xlim=c(-0.6,1.6), ylim=c(-0.9,1.0), pch=20, col='#11111120',
	cex=0.4, cex.lab=1.5, cex.axis=1.5, main='', xlab='r-i (mag)', ylab='i-z (mag)') 
par(mfrow=c(1,1))

# Quasar training set, N=2000 (Class 1)

download.file(url = "https://git.psu.edu/stat/astro/datasets/-/raw/master/SDSS_QSO.dat", destfile = "SDSS_QSO.dat")
qso1 <- read.table('SDSS_QSO.dat', header = TRUE)  
dim(qso1) ; summary(qso1)
bad_phot_qso <- which(qso1[,c(3,5,7,9,11)] > 21.0 | qso1[,3]==0)  # identify bad photometry
qso2 <- qso1[1:2000,-bad_phot_qso,]  # remove bad photometry
qso3 <- cbind((qso2[,3]-qso2[,5]), (qso2[,5]-qso2[,7]), (qso2[,7]-qso2[,9]), (qso2[,9]-qso2[,11])) # cbind concatenates columns
qso_train <- data.frame(cbind(qso3, rep(1, length(qso3[,1]))))
names(qso_train) <- c('u_g', 'g_r', 'r_i', 'i_z', 'Class')
dim(qso_train) ; summary(qso_train) 

# Star training set, N=5000 (Class 2)

download.file(url = "https://git.psu.edu/stat/astro/datasets/-/raw/master/SDSS_stars.csv", destfile = "SDSS_stars.csv")
temp2 <- read.csv('SDSS_stars.csv', header = TRUE)
dim(temp2) ; summary(temp2) 
star <- cbind((temp2[,1]-temp2[,2]), (temp2[,2]-temp2[,3]), (temp2[,3]-temp2[,4]), 
	(temp2[,4]-temp2[,5]))
star_train <- data.frame(cbind(star, rep(2, length(star[,1]))))
names(star_train) <- c('u_g','g_r','r_i','i_z','Class')
dim(star_train)  

# White dwarf training set, N=2000 (Class 3)

download.file(url = "https://git.psu.edu/stat/astro/datasets/-/raw/master/SDSS_wd.csv", destfile = "SDSS_wd.csv")
temp3 <- read.csv('SDSS_wd.csv', header = TRUE)
dim(temp3) ; summary(temp3)
temp3 <- na.omit(temp3)  # remove objects with missing data
wd <- cbind((temp3[1:2000,2]-temp3[1:2000,3]), (temp3[1:2000,3]-temp3[1:2000,4]),
	(temp3[1:2000,4]-temp3[1:2000,5]), (temp3[1:2000,5]-temp3[1:2000,6]))
wd_train <- data.frame(cbind(wd, rep(3, length(wd[,1]))))
names(wd_train) <- c('u_g', 'g_r', 'r_i', 'i_z', 'Class')
dim(wd_train) 

# Combine and plot the training set (9000 objects)

SDSS_train <- data.frame(rbind(qso_train, star_train, wd_train))  # rbind concatenates rows
names(SDSS_train) <- c('u_g', 'g_r', 'r_i', 'i_z', 'Class')
str(SDSS_train)

par(mfrow=c(1,3))  # plot training set in colors representing labeled classes
plot(SDSS_train[,1], SDSS_train[,2], xlim=c(-0.7,3), ylim=c(-0.6,1.6), pch=20, 
   	col=SDSS_train[,5], cex=0.6, cex.lab=1.6, cex.axis=1.6, main='Training dataset', xlab='u-g (mag)',
   	ylab='g-r (mag)')
legend(-0.5, 1.7, c('QSO','MS + RG','WD'), pch=20, col=c('black','red','green'), 
	cex=0.8)
plot(SDSS_train[,2], SDSS_train[,3], xlim=c(-0.6,1.6), ylim=c(-0.6,1.6), pch=20, 
	col=SDSS_train[,5], cex=0.6, cex.lab=1.6, cex.axis=1.6, main='', xlab='g-r (mag)',
	ylab='r-i (mag)') 
plot(SDSS_train[,3], SDSS_train[,4], xlim=c(-0.6,1.6), ylim=c(-0.9,1.0), pch=20, 
	col=SDSS_train[,5], cex=0.6, cex.lab=1.6, cex.axis=1.6, main='', xlab='r-i (mag)',
	ylab='i-z (mag)') 
par(mfrow=c(1,1))

## Function to evaluate classification performance

class_eval <- function(pred, act, plot=TRUE, ...){
  iact <- as.integer(act)
  ipred <- as.integer(pred)
  acc <- sum(ipred==iact)/length(iact)  # accuracy
  if (isTRUE(plot)){
    plot(jitter(ipred), jitter(iact), pch=20, cex=0.5, xlab='Predicted Class', ylab='True class',lab=c(3,3,1), ...)
    mtext(paste("Accuracy =", round(acc, 3)))
  }
  return(list("Confusion Table"=table("True Class"=iact, "Predicted Class"=ipred), Accuracy=acc))
}

## Copy original full dataset before splitting

SDSS_train_full <- SDSS_train

set.seed(456)
## Save 20% of training set for validation

val_set <- sample(nrow(SDSS_train), round(nrow(SDSS_train)*0.2))
SDSS_val <- SDSS_train_full[val_set,]
SDSS_train <- SDSS_train_full[-val_set,]

## Test of classification evaluation:  random class assignment

set.seed(123)
SDSS_rand_train_pred <- sample(1:3, nrow(SDSS_train), replace=T)
par(mfcol=c(1, 1))
class_eval(SDSS_rand_train_pred, SDSS_train$Class, main="Random Assignment")

# Unsupervised k-means partitioning

SDSS.kmean <- kmeans(SDSS_test,6)
print(SDSS.kmean$centers)
plot(SDSS_test[,1], SDSS_test[,2], pch=20, cex=0.3, col=gray(SDSS.kmean$cluster/7), 
   xlab='u-g (mag)', ylab='g-r (mag)', xlim=c(-0.5,3), ylim=c(-0.6,1.5), main='k-means clustering') 

# Linear discriminant analysis

library(MASS)
SDSS_lda <- lda(SDSS_train[,1:4], as.factor(SDSS_train[,5]))
SDSS_lda_train_pred <- predict(SDSS_lda)$class
SDSS_lda_val_pred <- predict(SDSS_lda, SDSS_val[,1:4])$class
SDSS_lda_test_pred <- predict(SDSS_lda, SDSS_test[,1:4])$class

par(mfrow=c(1,2))
plot(SDSS_val[,1],SDSS_val[,2], xlim=c(-0.7,3), ylim=c(-0.7,1.8), pch=20,
     col=SDSS_lda_val_pred, cex=0.5, main='', xlab='u-g (mag)', ylab='g-r (mag)')

class_eval(SDSS_lda_val_pred, SDSS_val$Class, main="LDA Classification")

# k-nn classification

library(class)
SDSS_knn_test_pred <- knn(SDSS_train[,1:4], SDSS_test, as.factor(SDSS_train[,5]), k=5, prob=T)
SDSS_knn_val_pred <- knn(SDSS_train[,1:4], SDSS_val[,1:4], as.factor(SDSS_train[,5]), k=5, prob=T)

plot(SDSS_val[,1], SDSS_val[,2], xlim=c(-0.7,3), ylim=c(-0.7,1.8), pch=20,
     col=SDSS_knn_val_pred, cex=0.5, main='', xlab='u-g (mag)', ylab='g-r (mag)')

class_eval(SDSS_knn_val_pred, SDSS_val$Class, main="kNN Classification")

# Neural network 

library(nnet)
SDSS_nnet <- nnet(as.factor(Class) ~ u_g + g_r + r_i + i_z, SDSS_train,size=5)

SDSS_nnet_train_pred <- predict(SDSS_nnet, type="class")
SDSS_nnet_val_pred <- predict(SDSS_nnet, SDSS_val, type="class")
SDSS_nnet_test_pred <- predict(SDSS_nnet, SDSS_test, type="class")

plot(SDSS_val[,1], SDSS_val[,2], xlim=c(-0.7,3), ylim=c(-0.7,1.8), pch=20,
     col=SDSS_nnet_val_pred, cex=0.5, main='', xlab='u-g (mag)', ylab='g-r (mag)')

class_eval(SDSS_nnet_val_pred, SDSS_val$Class, main="Neural Net Classification")

# Classification And Regression Tree

library(rpart)
SDSS_rpart <- rpart(SDSS_train[,5] ~., data=SDSS_train[,1:4], method="class")
summary(SDSS_rpart)
str(SDSS_rpart)

SDSS_rpart_train_pred <- predict(SDSS_rpart, type="class")
SDSS_rpart_val_pred <- predict(SDSS_rpart, SDSS_val, type="class")
SDSS_rpart_test_pred <- predict(SDSS_rpart, SDSS_test, type="class")

plot(SDSS_val[,1], SDSS_val[,2], xlim=c(-0.7,3), ylim=c(-0.7,1.8), pch=20,
   col=SDSS_rpart_val_pred, cex=0.5,
     main='', xlab='u-g (mag)', ylab='g-r (mag)')

class_eval(SDSS_rpart_val_pred, SDSS_val$Class, main="Tree Classification")

# Additional plots for decision tree

plot(SDSS_rpart, branch=0.5, margin=0.05)
text(SDSS_rpart, digits=3, use.n=T, cex=0.8)
plotcp(SDSS_rpart, lwd=2, cex.axis=1.3, cex.lab=1.3)

# Random Forests

library(randomForest)
SDSS_rf <- randomForest(as.factor(Class) ~ u_g + g_r + r_i + i_z, data=SDSS_train, mtry=2, importance=TRUE, do.trace=TRUE, ntree=100)

print(SDSS_rf)

SDSS_rf_train_pred <- predict(SDSS_rf, SDSS_train)
SDSS_rf_oob_pred <- predict(SDSS_rf)
SDSS_rf_val_pred <- predict(SDSS_rf, SDSS_val)
SDSS_rf_test_pred <- predict(SDSS_rf, SDSS_test)

plot(SDSS_val[,1], SDSS_val[,2], xlim=c(-0.7,3), ylim=c(-0.7,1.8), pch=20,
   col=SDSS_rf_val_pred, cex=0.5,
     main='', xlab='u-g (mag)', ylab='g-r (mag)')

class_eval(SDSS_rf_train_pred, SDSS_train$Class, main="Random Forest Classification")  
class_eval(SDSS_rf_oob_pred, SDSS_train$Class, main="Random Forest Classification")
class_eval(SDSS_rf_val_pred, SDSS_val$Class, main="Random Forest Classification")

# Support Vector Machine model, prediction and validation

library(e1071)
SDSS_svm <- svm((SDSS_train[,5]) ~.,data=SDSS_train[,1:4],cost=100, gamma=1)
SDSS_svm <- svm(as.factor(SDSS_train[,5]) ~.,data=SDSS_train[,1:4],cost=100, gamma=1)

summary(SDSS_svm)

SDSS_svm_train_pred <- predict(SDSS_svm)
SDSS_svm_val_pred <- predict(SDSS_svm, SDSS_val)
SDSS_svm_test_pred <- predict(SDSS_svm, SDSS_test)

plot(SDSS_test[,1], SDSS_test[,2], xlim=c(-0.7,3), ylim=c(-0.7,1.8), pch=20,
   col=round(as.numeric(SDSS_svm_test_pred)), cex=0.5, main='',
   xlab='u-g (mag)', ylab='g-r (mag)')

class_eval(SDSS_svm_val_pred, SDSS_val$Class, main="Random Forest Classification")

# Final plot of test set results

par(mfrow=c(1,3))
plot(SDSS_test[,1], SDSS_test[,2], xlim=c(-0.7,3), col=round(as.numeric(SDSS_svm_test_pred)), 
   ylim=c(-0.7,1.8), pch=20, cex=0.5, main='Support Vector Machine', xlab='u-g (mag)',ylab='g-r (mag)') 
plot(SDSS_test[,2], SDSS_test[,3], xlim=c(-0.7,1.8), col=round(as.numeric(SDSS_svm_test_pred)),
   ylim=c(-0.7,1.8), pch=20, cex=0.5, main='Classification of', xlab='g-r (mag)',ylab='r-i (mag)') 
plot(SDSS_test[,3], SDSS_test[,4], xlim=c(-0.7,1.8), col=round(as.numeric(SDSS_svm_test_pred)),
   ylim=c(-1.1,1.3), pch=20, cex=0.5, main='SDSS Test Dataset', xlab='r-i (mag)',ylab='i-z (mag)') 
par(mfrow=c(1,1))


# Write table of input data with output classifications onto disk

SDSS_test_svm_out <- cbind(SDSS[,6], SDSS[,7], SDSS_test, round(as.numeric(SDSS_svm_test_pred)))
names(SDSS_test_svm_out)[c(1,2,7)] <- c('R.A.', 'Dec', 'SVM Class')
write.table(format(SDSS_test_svm_out), file='SDSS_test_svm.out',sep='\t',quote=F)

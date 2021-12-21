# clear workspace:  
rm(list=ls()) 
if (!require("rjags")) install.packages("rjags")
library("rjags")


### Simulate data for testing
I = 36  # Number of participants
J = 6   # Number of items (ie, stimulus intensities / durations)

alpha = rnorm(I, 0, .5) # ability or latent true score
theta = abs(rnorm(I,0,1)) # improvement parameter/slope
beta = c(2,1,seq(0,-3, len=J-2)) # item/stimulus difficulty, center stimuli at 0 indicating duration/intensity at or slightly above chance performance

q = pmax(outer(alpha, beta, '-') * theta, 0) # this is going to be simulated true/latent dprime

# simulate the criterion
criterion = rnorm(I,0.5,1)

# number of signal trials for each subject and stimulation intensity
s = matrix(rep(c(420,420,100,100,100,100),I),nrow=I,ncol=J,byrow=T)
n = rep( 100, I) # and the noise trials

# initialize response matrices for detected stimuli and false alarms
y_f<- rep(0,I)
y_h<- matrix(rep(y_f,J),nrow=I, ncol=J, byrow=T)

# convert latent dprime and criterion to probabilities and simulate experiment
p_h <- pnorm(q/2 - matrix(rep(criterion,J),nrow=I,ncol=J)) # hit probability
p_f <- pnorm(-criterion) # false alarm probability
y_f <- rbinom(p_f,n,p_f) # random number of false alarms according to prior probabilities and sampling number
for (i in 1:J) {y_h[,i] <- rbinom(p_h[,i],s[1,i],p_h[,i])} # same for hits

# calculate empirical performance values
hr <- as.matrix(y_h/s)
fA <- as.matrix(y_f/n)

# log-linear transformation for extreme values (see Hautus et al.)
y_Yes <- cbind(y_f,y_h)
y_No <- cbind(n,s) - y_Yes
yY_hauted <-  y_Yes + 0.5
yN_hauted <- y_No + 0.5
s_hauted <- yY_hauted[,seq(2,ncol(y_Yes),1)] + yN_hauted[,seq(2,ncol(y_No),1)]
n_hauted <- yY_hauted[,1] + yN_hauted[,1]
hr_hauted <- as.matrix(yY_hauted[,seq(2,ncol(y_Yes),1)]/s_hauted)
fA_hauted <- as.matrix(yY_hauted[,1]/n_hauted)
dP <- qnorm(hr_hauted) - qnorm(matrix(rep(fA_hauted,ncol(y_h)),nrow=nrow(y_h),ncol=ncol(y_h)))

####### End simulation of data (see bottom of script for an alternative but equivalent simulation approach)

#### Begin analysis
# Needed variables:
# I: Number of participants
# J: Number of intensities/durations
# y_h: I by J matrix of correct judgments
# y_f: I by J matrix of false positive judgments
# s: I by J matrix of total signal judgments
# n: I by J matrix of total noise judgments
########

# parameter list to hand over to JAGS
forJags <- list(
  y_h = y_h,
  y_f = y_f,
  n = n,
  s = s,
  I = nrow(y_h),
  J = ncol(y_h)
)

# Number of MCMC iterations
## 4000 is only for testing
M = 40000

## compile JAGS
foo <- jags.model(file=paste0(getwd(),"/mac3dP.BUG"),
                  data=forJags)

# fit the model to the simulated data
out <- coda.samples(foo,M,
                    variable.names=c("alpha","beta","theta",
                                     "p_h","p_f","dprime0","c","prec_alpha","prec_theta"))

# extract parameters
cnames = colnames(out[[1]])
alpha_cols = out[[1]][,grepl("alpha[", cnames, fixed = TRUE)]
beta_cols = out[[1]][,grepl("beta[", cnames, fixed = TRUE)]
theta_cols = out[[1]][,grepl("theta[", cnames, fixed = TRUE)]
x_cols = out[[1]][,grepl("dprime0[", cnames, fixed = TRUE)]
c_cols = out[[1]][,grepl("c[",cnames,fixed = TRUE)]
p_h_cols = out[[1]][,grepl("p_h[",cnames,fixed = TRUE)]
p_f_cols = out[[1]][,grepl("p_f[",cnames,fixed = TRUE)]

# Compare with true values
plot(colMeans(alpha_cols), alpha)
abline(0,1)
# 
plot(colMeans(beta_cols), beta)
abline(0,1)
# 
plot(colMeans(theta_cols), theta)
abline(0,1)
#
plot(colMeans(x_cols), q)
abline(0,1)
#
plot(colMeans(c_cols), criterion)
abline(0,1)
#
plot(colMeans(p_h_cols), p_h)
abline(0,1)
#
plot(colMeans(p_f_cols), p_f)
abline(0,1)

# Chain diagnostic plots
pdf("diagnostics.pdf", ver="1.4")
plot(out[[1]])
dev.off()

## Posterior probability at chance
postPrAtChance <- sapply(1:M, function(i) x_cols[i,] == 0)
dim(postPrAtChance) <- c(I, J, M)
postPrAtChance <- apply(postPrAtChance, c(1,2), mean) 

## Plot of posterior probabilities by performance
plot(dP,postPrAtChance, col=rep(1:J,each=I),pch=19, ylab="Post. prob. at chance", xlab="Observed sensitivity")

for(i in 1:I)
  lines((dP)[i,],postPrAtChance[i,],lty=3, col="gray")
abline(h=.9,col="gray", lty=2)
abline(v=0, col="gray",lty=2)
legend(0.6,1,1:J,1:J,1:J,ncol = 6)

## Residuals
x_est = colMeans(x_cols)
dim(x_est) = c(I,J)
crit_est = colMeans(c_cols)
p_est = pnorm(x_est/2 - matrix(rep(crit_est,J),nrow=I,ncol=J))
std_resid = (y_h - s*p_est) / sqrt(s*p_est * ( 1 - p_est) )
#std_resid = (y_Yes - cbind(n,s)*p_est) / sqrt(cbind(n,s)*p_est * ( 1 - p_est) )

plot(p_est, std_resid, pch=19, col="black", ylab="Standardized residual", xlab="Transformed prediction")
for(i in 1:I)
  lines(p_est[i,],std_resid[i,],lty=1, col="gray")
abline(h=0, lty=2)
#abline(v=mean(p_est[,1]))
abline(h=c(-1.96,1.96), col="red", lty=2)


## posterior probability above chance
postPrAbChance <- sapply(1:M, function(i) x_cols[i,] > 0)
dim(postPrAbChance) <- c(I,J, M)
postPrAbChance <- apply(postPrAbChance, c(1,2), mean) 
# plot posterior probability above chance
plot(dP,postPrAbChance, col=rep(1:J,each=I),pch=19, ylab="Post. prob. above chance", xlab="Observed sensitivity")


for(i in 1:I)
  lines((dP)[i,],postPrAbChance[i,],lty=3, col="gray")
abline(h=.9,col="gray", lty=2)
abline(v=0, col="gray",lty=2)
legend(0.6,0.89,1:J,1:J,1:J,ncol = 6)


#######################################
# alternative data simulation approach#
#######################################
nTrials = 1340
noiseMean = 0
y_f<- rep(0,I)
y_h<- matrix(rep(y_f,J),nrow=I, ncol=J, byrow=T)

variance = 1;

stim = c(rep(0,100),rep(1,420),rep(2,420),rep(3,100),rep(4,100),rep(5,100),rep(6,100))
stim = sample(stim)
resp = rep(0,nTrials);
signalMean = q;


for (i_subj in 1:I)
{
  for (i in 1:nTrials) 
  {
    if (stim[i] == 0)  #noise trial
    {internalResponse = rnorm(1)*sqrt(variance)+noiseMean}
    if (stim[i] == 1) #sub1 trial
    {internalResponse = rnorm(1)*sqrt(variance)+signalMean[i_subj,1]}
    if (stim[i] == 2) #sub2 trial
    {internalResponse = rnorm(1)*sqrt(variance)+signalMean[i_subj,2]}
    if (stim[i] == 3) #ADTH trial
    {internalResponse = rnorm(1)*sqrt(variance)+signalMean[i_subj,3]}
    if (stim[i] == 4) #NTH1 trial
    {internalResponse = rnorm(1)*sqrt(variance)+signalMean[i_subj,4]}
    if (stim[i] == 5) #NTH2 trial
    {internalResponse = rnorm(1)*sqrt(variance)+signalMean[i_subj,5]}
    if (stim[i] == 6) #STH trial
    {internalResponse = rnorm(1)*sqrt(variance)+signalMean[i_subj,6]}
    
    #Decision:
    resp[i] = internalResponse>criterion[i_subj]
  }
  
  y_f[i_subj] = sum(resp==1 & stim ==0)
  
  y_h[i_subj,1] = sum(resp==1 & stim==1)
  y_h[i_subj,2] = sum(resp==1 & stim==2)
  y_h[i_subj,3] = sum(resp==1 & stim==3)
  y_h[i_subj,4] = sum(resp==1 & stim==4)
  y_h[i_subj,5] = sum(resp==1 & stim==5)
  y_h[i_subj,6] = sum(resp==1 & stim==6)
  
}

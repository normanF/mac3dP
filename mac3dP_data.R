# fitting behavioral 1-AFC data with the mac3dP model
# clear workspace:  
rm(list=ls()) 
library(rjags)

### load response data for all participants, requires detected yes|no counts
respY <- read.csv('YesResponses.csv',header=T,sep=",")
respN <- read.csv('NoResponses.csv',header=T,sep=",")

subjSelec <- c(1:36) # subject selection if necessary

#### Begin analysis
# Needed variables:
# I: Number of participants
# J: Number of intensities/durations
# y_h: I by J matrix of correct judgments
# y_f: I by J matrix of false positive judgments
# s: I by J matrix of total signal judgments
# n: I by J matrix of total noise judgments
########

y_h <- respY[subjSelec,seq(2,7,1)] # adapt indices pointing to the stimulation conditions in the dataframe, here 2 to 7 are signal conditions
y_f <- respY[subjSelec,1] # 1 is the noise condition (no stimulus)
n <- respY[subjSelec,1] + respN[subjSelec,1]
s <- respY[subjSelec,seq(2,7,1)] + respN[subjSelec,seq(2,7,1)] # adapt indices pointing to the stimulation conditions

# prepare input parameters for JAGS
forJags <- list(
  y_h = y_h,
  y_f = y_f,
  n = n,
  s = s,
  I = nrow(y_h),
  J = ncol(y_h)
)

# Number of MCMC iterations
## 4000 is only for testing, better 100k
M = 4000

## compile JAGS
foo <- jags.model(file=paste0(getwd(),"/mac3dP.BUG"),
                  data=forJags)

out <- coda.samples(foo,M,
                    variable.names=c("alpha","beta","theta",
                                     "p_h","p_f","dprime0","c","prec_alpha","prec_theta"))

# extract the important information: ability, difficulty, slope, dprime, criterion, probability for hits and fA
cnames = colnames(out[[1]])
alpha_cols = out[[1]][,grepl("alpha[", cnames, fixed = TRUE)]
beta_cols = out[[1]][,grepl("beta[", cnames, fixed = TRUE)]
theta_cols = out[[1]][,grepl("theta[", cnames, fixed = TRUE)]
x_cols = out[[1]][,grepl("dprime0[", cnames, fixed = TRUE)]
c_cols = out[[1]][,grepl("c[",cnames,fixed = TRUE)]
p_h_cols = out[[1]][,grepl("p_h[",cnames,fixed = TRUE)]
p_f_cols = out[[1]][,grepl("p_f[",cnames,fixed = TRUE)]

# Chain diagnostic plots
pdf("diagnostics.pdf", ver="1.4")
plot(out[[1]])
dev.off()

## Posterior probability at chance
postPrAtChance <- sapply(1:M, function(i) x_cols[i,] == 0)
dim(postPrAtChance) <- c(nrow(y_h), ncol(y_h), M)
postPrAtChance <- apply(postPrAtChance, c(1,2), mean) 

## Plot of posterior probabilities by log-linear transformed dprime
respY_Hauted = respY + 0.5
respN_Hauted = respN + 0.5
n_Hauted = respY_Hauted[subjSelec,1] + respN_Hauted[subjSelec,1]
s_Hauted = respY_Hauted[subjSelec,seq(2,7,1)] + respN_Hauted[subjSelec,seq(2,7,1)]
fA_Hauted = as.matrix(respY_Hauted[subjSelec,1]/n_Hauted)
hr_Hauted = as.matrix(respY_Hauted[subjSelec,seq(2,7,1)]/s_Hauted)
dP <- qnorm(hr_Hauted) - qnorm(matrix(rep(fA_Hauted,ncol(y_h)),nrow=nrow(y_h),ncol=ncol(y_h)))
criterion <- -qnorm(fA_Hauted)

# plot posterior probability against observed sensitivity
plot(dP,postPrAtChance, col=rep(1:ncol(y_h),each=nrow(y_h)),pch=19, ylab="Post. prob. at chance", xlab="Observed sensitivity")

for(i in 1:nrow(y_h))
  lines((dP)[i,],postPrAtChance[i,],lty=3, col="gray")
abline(h=.9,col="gray", lty=2) # credibility threshold
abline(v=0, col="gray",lty=2)
legend(0.6,1,1:ncol(y_h),1:ncol(y_h),1:ncol(y_h),ncol = 6)

# posterior probability above chance
postPrAbChance <- sapply(1:M, function(i) x_cols[i,] > 0)
dim(postPrAbChance) <- c(nrow(y_h),ncol(y_h), M)
postPrAbChance <- apply(postPrAbChance, c(1,2), mean) 

# plot posterior probability above chance
plot(dP,postPrAbChance, col=rep(1:ncol(y_h),each=nrow(y_h)),pch=19, ylab="Post. prob. above chance", xlab="Observed sensitivity")

for(i in 1:nrow(y_h))
  lines((dP)[i,],postPrAbChance[i,],lty=3, col="gray")
abline(h=.9,col="gray", lty=2)
abline(v=0, col="gray",lty=2)
legend(0.6,0.89,1:ncol(y_h),1:ncol(y_h),1:ncol(y_h),ncol = 6)

## Residuals
I = nrow(y_h)
J = ncol(y_h)
x_est = colMeans(x_cols)
dim(x_est) = c(I,J)
crit_est = colMeans(c_cols)
p_est = as.matrix(pnorm(x_est/2 - matrix(rep(crit_est,J),nrow=I,ncol=J)))
std_resid = (as.matrix(y_h) - as.matrix(s)*p_est) / sqrt(as.matrix(s)*p_est * ( 1 - p_est) )

plot(p_est, std_resid, pch=19, col="black", ylab="Standardized residual", xlab="Transformed prediction")
for(i in 1:I)
  lines(p_est[i,],std_resid[i,],lty=1, col="gray")
abline(h=0, lty=2)
abline(h=c(-1.96,1.96), col="red", lty=2)


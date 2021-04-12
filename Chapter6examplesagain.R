library(astsa)
set.seed(999); num=100
x=arima.sim(n=num+1,list(ar=.8),sd=1)
y=ts(x[-1]+rnorm(num,0,1))
u=ts.intersect(y,lag(y,-1),lag(y,-2))
varu=var(u)
coru=cor(u)
varu;coru           
phi=coru[1,3]/coru[1,2]; phi
(phi0=varu[1,3]/varu[1,2])
(q=(1-phi^2)*varu[1,2]/phi)
(r=varu[1,1]-q/(1-phi^2))
(r=varu[1,1]-varu[1,2]/phi)

dat=data.frame(x=c(1,2,3,4,5,6),  y=c(1,3,5,6,8,12))
min.RSS <- function(data, par) {
  with(data, sum((par[1] + par[2] * x - y)^2))
}
(result <- optim (par = c(0, 1), fn = min.RSS, data = dat))




library(astsa)
data(package="astsa")
data(globtemp)
globtemp=astsa::globtemp
globtemp1=astsa::globtemplglobtemp1=data(globtempl)
# Setup
y = cbind(globtemp, globtempl); num = nrow(y); input = rep(1,num)
A = array(rep(1,2), dim=c(2,1,num))
mu0 = -.35; Sigma0 = 1; Phi = 1
# Function to Calculate Likelihood
Linn = function(para){
  cQ = para[1] # sigma_w
  cR1 = para[2] # 11 element of chol(R)
  cR2 = para[3] # 22 element of chol(R)
  cR12 = para[4] # 12 element of chol(R)
  cR = matrix(c(cR1,0,cR12,cR2),2) # put the matrix together
  drift = para[5]
  kf = Kfilter1(num,y,A,mu0,Sigma0,Phi,drift,0,cQ,cR,input)
  return(kf$like) }
# Estimation
init.par = c(.1,.1,.1,0,.05) # initial values of parameters

est = optim(init.par, Linn, gr=NULL, method="BFGS", hessian=TRUE)
SE = sqrt(diag(solve(est$hessian)))
cbind(estimate=c(phi=est$par[1],sigw=est$par[2],sigv=est$par[3]), SE) 
SE = sqrt(diag(solve(est$hessian)))
# Display estimates
u = cbind(estimate=est$par, SE)
rownames(u)=c('sigw','cR11', 'cR22', 'cR12', 'drift'); u
# Smooth (first set parameters to their final estimates)
cQ = est$par[1]
cR1 = est$par[2]
cR2 = est$par[3]
cR12 = est$par[4]
cR = matrix(c(cR1,0,cR12,cR2), 2)
(R = t(cR)%*%cR) # to view the estimated R matrix
drift = est$par[5]
ks = Ksmooth1(num,y,A,mu0,Sigma0,Phi,drift,0,cQ,cR,input)
# Plot
xsm = ts(as.vector(ks$xs), start=1880)
rmse = ts(sqrt(as.vector(ks$Ps)), start=1880)
plot(xsm, ylim=c(-.6, 1), ylab='Temperature Deviations')
xx = c(time(xsm), rev(time(xsm)))
yy = c(xsm-2*rmse, rev(xsm+2*rmse))
polygon(xx, yy, border=NA, col=gray(.6, alpha=.25))
lines(globtemp, type='o', pch=2, col=4, lty=6)
lines(globtempl, type='o', pch=3, col=3, lty=6)


library(nlme) # loads package nlme
# Generate data (same as Example 6.6)
set.seed(999); num = 100
x = arima.sim(n=num+1, list(ar = .8), sd=1)
y = ts(x[-1] + rnorm(num,0,1))
# Initial Estimates (same as Example 6.6)
u = ts.intersect(y, lag(y,-1), lag(y,-2))
varu = var(u); coru = cor(u)
phi = coru[1,3]/coru[1,2]
q = (1-phi^2)*varu[1,2]/phi
r = varu[1,1] - q/(1-phi^2)
# EM procedure - output not shown

(em = EM0(num, y, A=1, mu0=0, Sigma0=2.8, Phi=phi, cQ=sqrt(q), cR=sqrt(r),
          max.iter=75, tol=.00001))
# Standard Errors (this uses nlme)
phi = em$Phi; cq = sqrt(em$Q); cr = sqrt(em$R)
mu0 = em$mu0; Sigma0 = em$Sigma0
para = c(phi, cq, cr)
Linn = function(para){ # to evaluate likelihood at estimates
  kf = Kfilter0(num, y, 1, mu0, Sigma0, para[1], para[2], para[3])
  return(kf$like) }
emhess = fdHess(para, function(para) Linn(para))
SE = sqrt(diag(solve(emhess$Hessian)))
# Display Summary of Estimation
estimate = c(para, em$mu0, em$Sigma0); SE = c(SE, NA, NA)
u = cbind(estimate, SE)
rownames(u) = c('phi','sigw','sigv','mu0','Sigma0'); u






library(tsm)
data(sarb_quarter)

dat <- sarb_quarter$KBP6006L/ sarb_quarter$KBP6006D
dat.tmp <- diff(log(na.omit(dat)) * 100, lag = 1)

inf <- ts(dat.tmp, start = c(1960, 2), frequency = 4)
plot.ts(inf)
###To remove some data from the time series for inflation, we can set a number of observations equal to NA, which in this case relates to observations 70 to 82.
inf.mis <- inf
inf.mis[70:82] <- NA
plot(inf.mis)

###We can use a local level model structure again, which may be constructed as before, where the model will now be applied to the variable inf.mis.
fn <- function(parm) {
  dlmModPoly(order = 1, dV = exp(parm[1]), dW = exp(parm[2]))
}

## Estimate parameters & generate statistics
fit <- dlmMLE(inf.mis, rep(0, 2), build = fn, hessian = TRUE)
conv <- fit$convergence  # zero for converged

loglik <- dlmLL(inf.mis, dlmModPoly(1))
n.coef <- 2
r.aic <- (2 * (loglik)) + 2 * (sum(n.coef))  #dlmLL caculates the neg. LL
r.bic <- (2 * (loglik)) + (log(length(inf))) * (n.coef)

mod <- fn(fit$par)
obs.error.var <- V(mod)
state.error.var <- W(mod)

filtered <- dlmFilter(inf.mis, mod = mod)
smoothed <- dlmSmooth(filtered)
resids <- residuals(filtered, sd = FALSE)
mu <- dropFirst(smoothed$s)

conf.tmp <- dlmSvd2var(smoothed$U.S, smoothed$D.S)
conf <- ts(as.numeric(conf.tmp)[-1], start = c(1960, 2), 
           frequency = 4)
wid <- qnorm(0.05, lower = FALSE) * sqrt(conf)

conf.pos <- mu + wid
conf.neg <- mu - wid

comb.state <- cbind(mu, conf.pos, conf.neg)

cat("AIC", r.aic)
## AIC 401.036
cat("BIC", r.bic)
## BIC 407.8947
cat("V.variance", obs.error.var)
## V.variance 1.365551
cat("W.variance", state.error.var)
## W.variance 0.03032414
###To forecast forward we use the dlmForecast command, where in this case we are interested in forecasting ten steps ahead.
forecast <- dlmForecast(filtered, nAhead = 12)
var.2 <- unlist(forecast$Q)
wid.2 <- qnorm(0.05, lower = FALSE) * sqrt(var.2)
comb.fore <- cbind(forecast$f, forecast$f + wid.2, forecast$f - 
                     wid.2)

result <- ts(rbind(comb.state, comb.fore), start = c(1960, 
                                                     2), frequency = 4)

###To draw a graph that contains the results, we would utilise the following commands:
  par(mfrow = c(1, 1), mar = c(2.2, 2.2, 1, 1), cex = 0.8)
plot.ts(comb.state, col = c("black", "red", "red"), plot.type = "single", 
        xlab = "", ylab = "", lty = c(1, 2, 2), ylim = c(-2, 
                                                         8))
lines(inf.mis, col = "darkgrey", lwd = 1.5)
legend("topright", legend = c("Observed Deflator", "Stochastic level", 
                              "Confidence"), lwd = c(1.5, 1), col = c("darkgrey", 
                                                          
                                                                    
                                                              
library(astsa)                                                                      
data(WBC)    ; data(HCT); data(PLT)                               
y    = cbind(WBC, PLT, HCT)
num  = nrow(y)       
A    = array(0, dim=c(3,3,num))  # creates num 3x3 zero matrices
for(k in 1:num) if (y[k,1] > 0) A[,,k]= diag(1,3) 

# Initial values 
mu0    = matrix(0,3,1) 
Sigma0 = diag(c(.1,.1,1) ,3)
Phi    = diag(1,3)
cQ     = diag(c(.1,.1,1), 3)
cR     = diag(c(.1,.1,1), 3)  
(em = EM1(num, y, A, mu0, Sigma0, Phi, cQ, cR, 100, .001))    

# Graph smoother
ks  = Ksmooth1(num, y, A, em$mu0, em$Sigma0, em$Phi, 0, 0, chol(em$Q), chol(em$R), 0)
ks2  = Ksmooth2(num, y, A, em$mu0, em$Sigma0, em$Phi, 0,0,0, chol(em$Q), chol(em$R),0,0)
y1s = ks$xs[1,,] 
y2s = ks$xs[2,,] 
y3s = ks$xs[3,,]
p1  = 2*sqrt(ks$Ps[1,1,]) 
p2  = 2*sqrt(ks$Ps[2,2,]) 
p3  = 2*sqrt(ks$Ps[3,3,])
par(mfrow=c(3,1))
tsplot(WBC, type='p', pch=19, ylim=c(1,5), xlab='day')
lines(y1s) 
lines(y1s+p1, lty=2, col=4) 
lines(y1s-p1, lty=2, col=4)
tsplot(PLT, type='p', ylim=c(3,6), pch=19, xlab='day')
lines(y2s)
lines(y2s+p2, lty=2, col=4)
lines(y2s-p2, lty=2, col=4)
tsplot(HCT, type='p', pch=19, ylim=c(20,40), xlab='day')
lines(y3s)
lines(y3s+p3, lty=2, col=4) 
lines(y3s-p3, lty=2, col=4) 





##### Shumway, Stoffer, 6-13        Shumway, Stoffer, 6-13  ##### 
##### Shumway, Stoffer, 6-13        Shumway, Stoffer, 6-13  ##### 

##### Shumway, Stoffer, 6-13        Shumway, Stoffer, 6-13  ##### 
##### Shumway, Stoffer, 6-13        Shumway, Stoffer, 6-13  ##### 

library(plyr) # used for displaying progress
tol = sqrt(.Machine$double.eps) # determines convergence of optimizer

nboot = 500 # number of bootstrap replicates
y = window(qinfl, c(1953,1), c(1965,2)) # inflation
z = window(qintr, c(1953,1), c(1965,2)) # interest
num = length(y)
A = array(z, dim=c(1,1,num))
input = matrix(1,num,1)
# Function to Calculate Likelihood
Linn = function(para, y.data){ # pass data also
  phi = para[1]; alpha = para[2]
  b = para[3]; Ups = (1-phi)*b
  cQ = para[4]; cR = para[5]
  kf = Kfilter2(num,y.data,A,mu0,Sigma0,phi,Ups,alpha,1,cQ,cR,0,input)
  return(kf$like) }
# Parameter Estimation
mu0 = 1; Sigma0 = .01
init.par = c(phi=.84, alpha=-.77, b=.85, cQ=.12, cR=1.1) # initial values
est = optim(init.par, Linn, NULL, y.data=y, method="BFGS", hessian=TRUE,
            control=list(trace=1, REPORT=1, reltol=tol))
SE = sqrt(diag(solve(est$hessian)))
phi = est$par[1]; alpha = est$par[2]
b = est$par[3]; Ups = (1-phi)*b
cQ = est$par[4]; cR = est$par[5]
round(cbind(estimate=est$par, SE), 3)
###estimate SE
####  phi 0.865 0.223
##### b 0.788 0.226
##### cQ 0.115 0.107
##### cR 1.135 0.147
# BEGIN BOOTSTRAP
# Run the filter at the estimates
kf = Kfilter2(num,y,A,mu0,Sigma0,phi,Ups,alpha,1,cQ,cR,0,input)
# Pull out necessary values from the filter and initialize
xp = kf$xp
innov = kf$innov
sig = kf$sig
K = kf$K
e = innov/sqrt(sig)
e.star = e # initialize values
y.star = y
xp.star = xp
k = 4:50 # hold first 3 observations fixed
para.star = matrix(0, nboot, 5) # to store estimates
init.par = c(.84, -.77, .85, .12, 1.1)
pr <- progress_text() # displays progress
pr$init(nboot)
for (i in 1:nboot){
  pr$step()
  e.star[k] = sample(e[k], replace=TRUE)
  for (j in k){ xp.star[j] = phi*xp.star[j-1] +
    Ups+K[j]*sqrt(sig[j])*e.star[j] }
  y.star[k] = z[k]*xp.star[k] + alpha + sqrt(sig[k])*e.star[k]
  est.star = optim(init.par, Linn, NULL, y.data=y.star, method="BFGS",
                   control=list(reltol=tol))
  para.star[i,] = cbind(est.star$par[1], est.star$par[2], est.star$par[3],
                        abs(est.star$par[4]), abs(est.star$par[5])) }

# Some summary statistics
rmse = rep(NA,5) # SEs from the bootstrap
for(i in 1:5){rmse[i]=sqrt(sum((para.star[,i]-est$par[i])^2)/nboot)
cat(i, rmse[i],"\n") }
# Plot phi and sigw
phi = para.star[,1]
sigw = abs(para.star[,4])
phi = ifelse(phi<0, NA, phi) # any phi < 0 not plotted
library(psych) # load psych package for scatter.hist
scatter.hist(sigw, phi, ylab=expression(phi), xlab=expression(sigma[~w]),
             smooth=FALSE, correl=FALSE, density=FALSE, ellipse=FALSE,
             title='', pch=19, col=gray(.1,alpha=.33),
             panel.first=grid(lty=2), cex.lab=1.2)







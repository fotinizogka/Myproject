# R CODE:Goodness-of-fit tests for the Poisson distribution (bootstrap method)

library(digest)
library(energy)#for SR
library(stats)#for uniform 
library(distr6)#for negative binomial
library(extraDistr)#for discrete uniform
library(HMMpa)#for generalized poisson
library(VGAM)
library(remotes)

#Usquare
U.stat.test <- function(x)
{
  var(x) #sample variance
  mean(x)
  FI <- var(x)/mean(x)
  FI
  Dn <- (n-1)*FI
  Dn
  u <- 1/(2*n)*(Dn-n)^2
  return(u)
}
#Kolmogorov Smirnov
KS.stat.test <- function(x)
{
  options(warn = -1)
  if (sum(x) == 0) {0}
  else {ecdfx = ecdf(x)
  f1 = function(y) {(ecdfx(y)-ppois(y, mean(x)))}
  resf11 = f1(min(x):max(x))
  ks = max(resf11)
  return(ks)
  }
}
#Cramer von Mises
Cn.stat.test <- function(x)
{
  if (sum(x) == 0) {0}
  else {ecdfx = ecdf(x)
  f1 = function(y) {(ecdfx(y)-ppois(y,mean(x)))^2*dpois(y,mean(x))}
  resf11 = sum(f1(min(x):max(x)))
  cn = length(x)*resf11
  return(cn)
  }
}
#Cramer von Mises Cn*
Cnstar.stat.test <- function(x)
{
  options(warn = -1)
  if (sum(x) == 0) {0}
  else {ecdfx = ecdf(x)
  f1 = function(y) {(ecdfx(y)-ppois(y,mean(x)))^2*(ecdfx(y)-ecdfx(y-1))}
  resf11 = sum(f1(min(x):max(x)))
  Cnstar = length(x)*resf11
  }
}
#Klar from Cramer von Mises
Klar.stat.test <- function(x)
{
  options(warn = -1)
  if (sum(x) == 0) {0}
  else {ecdfx = ecdf(x)
  f1 = function(y) {sqrt(length(x))*abs((ecdfx(y) - ppois(y, mean(x))))}
  resf11 = sum(f1(min(x):max(x)))
  tn = return(resf11)
  }
}
#Klar
In.stat.test <- function(x)
{
  options(warn = -1)
  if (sum(x) == 0) {0}
  else {ecdfx = ecdf(x)
  f1 = function(y) {abs (ecdfx(y)-ppois(y, mean(x)))}
  resf11 = sum(f1(min(x):max(x-1)))
  i. = sqrt(length(x))*max(resf11)
  return(i.)
  }
}
#Kocherlakota
K.stat.test <- function(t,lambda)
{
  t <- runif(1,min=0,max=1) 
  t
  g <- function(t)
  {
    lambda = mean(x)
    exp(lambda*(t-1))
  }
  gn <- function(t)
  {
    1/n*sum(t^x)
  }
  gn(t)
  sqrt(n)*(gn(t)-g(t))/exp(lambda*(t^2-1))-exp(2*lambda*(t-1))*(1+lambda*(t-1)^2)
}
#Rueda et.al.
R.stat.test <- function(x)
{
  Gn = function(t,z){
    (sqrt(length(z))*((1/length(z))*sum(t^z) -exp(mean(z)*(t-1))))^2
  }
  i = integrate(f = Vectorize(Gn,vectorize.args = 't'), lower = 0, upper = 1,z = x)
  i$value
}
#Baringhaus Gurtler Henze
Ra.stat.test <- function(x,a)
{
  Gna = function(t,y,b){
    (sqrt(length(y))*((1/length(y))*sum(t^y) -exp(mean(y)*(t-1))))^2*t^b
  }
  i = integrate(f = Vectorize(Gna,vectorize.args = 't'), lower = 0, upper = 1, y = x, b = a)
  i$value
}
#Baringhaus Henze
T.stat.test <- function(x)
{
  suma = outer(x,x,"+")
  mult = outer(x,x,"*")
  t. = table(x)
  dd = as.numeric(names(t.))+1
  ni = rep(0,(max(x)+1))
  ni[dd] = as.numeric(t.)
  msum = suma[mult>0]
  mmult = mult[mult>0]
  t1 = (1/length(x))*(sum(mean(x)^2/(suma+1))+sum(mmult/ (msum-1))) mean(x)*(length(x)-(1/length(x))*(ni[1])^2)
  return(t1)
}
#Treutler
Ta.stat.test <- function(x,a)
{
  f2nder = function(t) {mean(x) * mean(sum(t^x))}
  f2der = function(t) {mean(sum(t^x*x/t))}
  f3_a = function(t,b) {((f2nder(t) - f2der(t)) ^2) * t^b}
  b_a = integrate(Vectorize(f3_a, vectorize.args = 't'), lower = 0, upper = 1, b = a)
  b_as = as.numeric(b_a[1])
  return(b_as/length(x))
}
#Meintanis Nikitin
MNast.stat.test <- function(x,a)
{
  barx <- mean(x)
  gx<-function(x){(x/(x+a)) -(barx/(x+a+1)) 
    -(x-barx)*(1/barx)*mean(x/(x+a))}
  
  var.=mean((gx(x))^2)
  ÌÍ= ((sqrt(1/n))*sum(x/(x+a)-barx/(x+a+1)))/sqrt(var.)
  
  
  return(ÌÍ)
}
#W^2
W2.stat.test = function(x)
{
  options(warn = -1)
  if (sum(x) == 0) {0}
  else {ecdfx = ecdf(x)
  f1 = function(y) {(ecdfx(y)-ppois(y, mean(x)))^2*
      (ppois(y,mean(x))-ppois(y-1,mean(x)))}
  resf11 = sum(f1(min(x):max(x)))
  w2 = length(x)*resf11
  }
}
#A^2
A2.stat.test <- function(x)
{
  options(warn = -1)
  if (sum(x) == 0) {0}
  else {ecdfx=ecdf(x)
  f1 = function(y) {((ecdfx(y) - ppois(y, mean(x)))^2)/( ppois(y, mean(x))*(1-ppois(y, mean(x))))*(ppois(y, mean(x)) - ppois(y-1, mean(x)))}
  resf11 = sum(f1(min(x):max(x)))
  a2 = length(x)*resf11
  }
}
#Szekely Rizzo

#SR = poisson.m(x)
#Nakamura Perez Abreu
NPA.stat.test<-function(x)
{
  np<-0
  for(i in 1:n){
    for(j in 1:n){
      for(k. in 1:n){
        for(l. in 1:n){
          
          t<-runif(1,0,1)
          
          I=sum(t^(x[i]+x[j]))
          
          
          np=(np+(x[i]*(x[i]-x[j]-1)*x[k.]*(x[k.]-x[l.]-1))*I)
          
        }
      }
      
    }
  }
  npa=(length(x)*(1/length(x)^4)*np)/mean(x)^1.45
  return(npa)
}

n=50 #or n=100, n=200
#For Poisson distribution
lambda=0.5 #or lambda=1,lambda=5,lambda=10,lambda=30
s=0;z.=0;m=0;b.=0;c=0;d=0;e=0;f=0;h=0;k=0;l=0;o=0;
r=0;w=0;u=0;q
n.boot=500 #replications
a=5 #parametre 
#simulation runs
simulations=10000
for(i in 1:simulations)
{
  x <- rpois(n,lambda) #poisson distribution
  #x<-rbinom(n,2,0.5) #binomial distribution
  #x<-rnbinom(n,1,0.5)#negative binomial distribution
  
  #x<-rdunif(n,0,2)   #discrete uniform distribution
  #x<-rgenpois(n,4,0.1)#generalized poisson distribution
  #x<-rmixpois(n,c(1,5),c(0.25,0.75)) #poisson mixtures
  lambda.hat = mean(x)
  
  Usq = U.stat.test(x) #usquare
  KS = KS.stat.test(x) #ks
  Cn = Cn.stat.test(x) #cm
  Cnstar = Cnstar.stat.test(x) #cn*
  Tn = Klar.stat.test(x) #Klar Tn (cm)
  K = K.stat.test(t,lambda.hat) #kocherlakota
  In = In.stat.test(x) #klar In
  R = R.stat.test(x) #rueda
  Ra = Ra.stat.test(x,a) #barin gurtler henze
  T = T.stat.test(x) #barin henze
  Ta = Ta.stat.test(x,a) #treutler
  MNa =abs(MNast.stat.test(x,a)) #MNa
  W2 = W2.stat.test(x) #w2 cm
  A2 = A2.stat.test(x) #A2 cm
  SR = poisson.m(x) #szekely rizzo
  NPA=NPA.stat.test(x)#nakamura 
  #bootstrap 
  T.boot1 = rep(n.boot,0)
  T.boot2 = rep(n.boot,0)
  T.boot3 = rep(n.boot,0)
  T.boot4 = rep(n.boot,0)
  T.boot5 = rep(n.boot,0)
  T.boot6 = rep(n.boot,0)
  T.boot7 = rep(n.boot,0)
  T.boot8 = rep(n.boot,0)
  T.boot9 = rep(n.boot,0)
  T.boot10= rep(n.boot,0)
  T.boot11= rep(n.boot,0)
  T.boot12= rep(n.boot,0)
  T.boot13= rep(n.boot,0)
  T.boot14= rep(n.boot,0)
  T.boot15= rep(n.boot,0)
  T.boot16= rep(n.boot,0)
  
  #bootstrap replications  
  for (j in 1:n.boot){
    x.boot = rpois(n,lambda.hat)
    T.boot1[j] = U.stat.test(x.boot)
    T.boot2[j] = KS.stat.test(x.boot)
    T.boot3[j] = Cn.stat.test(x.boot)
    T.boot4[j] = Cnstar.stat.test(x.boot)
    T.boot5[j] = Klar.stat.test(x.boot)
    T.boot6[j] = K.stat.test(t,lambda.hat)
    T.boot7[j] = In.stat.test(x.boot)
    T.boot8[j] = R.stat.test(x.boot)
    T.boot9[j]= Ra.stat.test(x.boot,a)
    T.boot10[j]= T.stat.test(x.boot)
    T.boot11[j]= Ta.stat.test(x.boot,a)
    T.boot12[j]= MNast.stat.test(x.boot,a)
    T.boot13[j]= W2.stat.test(x.boot)
    T.boot14[j]= A2.stat.test(x.boot)
    T.boot15[j]= poisson.m(x.boot)
    T.boot16[j]= NPA.stat.test(x.boot)
  }
  
  if (Usq>quantile(na.omit(T.boot1),0.95)) s=s+1
  if (KS>quantile(T.boot2,0.95)) z.=z.+1
  if (Cn>quantile(T.boot3,0.95)) m=m+1
  if (Cnstar>quantile(T.boot4,0.95)) b.=b.+1
  if (Tn>quantile(T.boot5,0.95)) c=c+1
  if (K>quantile(T.boot6,0.95)) d=d+1
  if (In>quantile(T.boot7,0.95)) e=e+1
  if (R>quantile(T.boot8,0.95)) f=f+1
  if (Ra>quantile(T.boot9,0.95)) h=h+1
  if (T>quantile(T.boot10,0.95)) k=k+1
  if (Ta>quantile(T.boot11,0.95)) l=l+1
  if (MNa>quantile(T.boot12,0.95)) o=o+1
  if (W2>quantile(T.boot13,0.95)) r=r+1
  if (A2>quantile(T.boot14,0.95)) w=w+1
  if (SR>quantile(T.boot15,0.95)) u=u+1
  if (NPA>quantile(T.boot18,0.95)) q=q+1
}
s/simulations #usquare
z./simulations #ks
m/simulations #cm
b./simulations #cn*
c/simulations #tn klar
d/simulations #kocherlakota
e/simulations #klar In
f/simulations #rueda
h/simulations #barin gurtler henze
k/simulations #barin henze
l/simulations #treutler
o/simulations #MNa
r/simulations #W2 from (cn)
w/simulations #A2 from (cn)
u/simulations #Szekely Rizzo
q/simulations #Nakamura Perez Abreu

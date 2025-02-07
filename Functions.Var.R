# Initial Functions
library(MASS)
library(MASS)

########################################  Functions for the code ######################################
## Calculates c_hat(k1,k2) =  (1/(2M+1))* sum_s(J_(k1+s)J*_(k2+s))
c_hatM = function(J,k1,k2,M){
  s = seq(-M,M)
  n = length(J)
  ind10 = k1+s
  ind20 = k2+s 
  ind1 = ((ind10) %% n)+1
  ind2 = ((ind20) %% n)+1
  Jk = J[ind1]
  Jkl = Conj(J[ind2])
  ch = mean(Jk*Jkl)
  return(ch)
}


## Calculates Yhat(l1,l2) = J_(k+l1+s,n)J*_(k+l2+s,n) - c_hat(k,l1,l2)
Yhat = function(J,k,l1,l2,M){
  s = seq(-M,M)
  n = length(J)
  ind1 = k+l1+s
  ind2 = k+l2+s
  ind1 = ((ind1) %% n)+1
  ind2 = ((ind2) %% n)+1
  J1 = J[ind1]
  J2 = Conj(J[ind2])
  Y = J1*J2 - c_hatM(J,k+l1,k+l2,M)
  return(Y)
}

## Calsulates S_1M and S_2M
## Inputs:
##        J = fft
##        k = leading freq
##        l = Lag
##        L, M = Tuning parameters
##        nP = n/P 

S_var_c_estP = function(J,k,l,L,M,nP){
  Y01 = Yhat(J,k,0,l,M)
  Y23 = Y01
  s = seq(-M,M)
  
  ##Tmat[s1,s2] = Yhat_s1  x Yhat^*_s2
  Tmat = outer(Y01, Conj(Y23),"*")
  
  
  P = length(J)/nP
  
  
  RL0 = c(seq(-L,L)) 
  
  RL = RL0
  
  for(j in 1:P)
  {
    RL = c(RL,seq(j*nP-L,j*nP+L),seq(-j*nP-L,-j*nP+L))
  } #RL should pick jn_P +/- L
  
  O1 = outer(s,s,"-")
  O2 = outer(s,s,"+")+(2*k+l)
  
  I1 <- `dim<-`(O1 %in% RL0, dim(O1)) #I1[s1,s2]=TRUE is |s_1 - s_2| <= L 
  I2 <- `dim<-`(O2 %in% RL, dim(O2)) #I2[s1,s2]=TRUE if |s_1 + s+2 + 2k +l|_nP <=L
  
  
  That1 = sum(I1*Tmat)/((2*M+1)^2) #S_1,M
  
  That2 = sum(I2*Tmat)/((2*M+1)^2) #S_2,M
  
  return(c(That1,That2))
}

T_var_c_conj_estP = function(J,k,l,L,M,nP){
  Y01 = Yhat(J,k,0,l,M)
  Y23 = Yhat(J,k,l,0,M)
  s = seq(-M,M)
  Tmat = outer(Y01, Conj(Y23),"*")
  
  #That1 = 0
  #That2 = 0
  
  P = length(J)/nP
  
  
  RL0 = c(seq(-L,L))
  RL = RL0
  
  for(j in 1:P)
  {
    RL = c(RL,seq(j*nP-L,j*nP+L),seq(-j*nP-L,-j*nP+L))
  }
  
  
  O1_1 = outer(s,s,"-")+l
  O1_2 = outer(s,s,"-")-l
  O2_1 = outer(s,s,"+")+(2*k)
  O2_2 = outer(s,s,"+")+(2*k + 2*l)
  
  I1_1 = `dim<-`(O1_1 %in% RL, dim(O1_1))
  I1_2 = `dim<-`(O1_2 %in% RL, dim(O1_2))
  I2_1 = `dim<-`(O2_1 %in% RL, dim(O2_1))
  I2_2 = `dim<-`(O2_2 %in% RL, dim(O2_2))
  
  I1 = I1_1*I1_2
  I2 = I2_1*I2_2
  
  That1 = sum(I1*Tmat)/((2*M+1)^2) #T_1,M
  That2 = sum(I2*Tmat)/((2*M+1)^2) #T_2,M
  
  return(c(That1,That2))
}


### Functions to simulate time series (all taken from SSR's rmd file)

ar2complex = function(lambda,omega,n1){
  n = (n1+200)
  gen = rnorm(n)
  x = rep(0,n)
  for(j in c(3:n)){
    x[j] = (2*lambda*cos(omega)*x[j-1] - (lambda**2)*x[j-2] + gen[j])
  }
  x1 = x[-c(1:200)]  
  return(x1)
}

tvar1 = function(beta0,beta1,n1){
  gen  = rnorm(n1)    
  x = gen
  phi = 0.9-0.7/(1+exp(beta0*(x-beta1)))
  for(j in c(2:n1)){
    x[j] = (phi[j]*x[j-1] + gen[j])
  }
  return(x)
}

tvar2 = function(beta0,beta1,omega,n1){
  gen = rnorm(n1)
  x = c(1:n1)
  # Since |phi|<1 this ensures the transition matrix    
  # has largest eigenvalue less than one; hence it is stable.
  phi = 0.9-0.7/(1+exp(beta0*(x-beta1)))
  for(j in c(3:n1)){
    x[j] = (2*phi[j]*cos(omega)*x[j-1] - (phi[j]**2)*x[j-2] + gen[j])
  }
  return(x)
}

VARP = function(n1){
  n2 = n1+200
  B1<-matrix(c(0.9, 0.4, 0.6, 0, - 0.9, 0.2, 0,0,0.7), 3 , byrow = TRUE)
  # B1 is a VAR transition matrix, its eigenvalues should be less than one
  # to ensure it is causal.   
  mu1 = rep(0,3)
  Sigma1 = matrix(c(1,0,0,0,1,0,0,0,1),3,byrow = TRUE)
  gen = mvrnorm(n=n2,mu = mu1, Sigma  = Sigma1)
  x = gen
  for(j in c(2:n2)){
    x[j,] = B1%*%x[j-1,] + gen[j,]
  }
  x1 = x[-c(1:200),]  
  x2 = c(t(x1))
  return(x2)
}

tvVARP = function(beta0,beta1,n1){
  B1<-matrix(c(0.9, 0.4, 0.6, 0, -0.9, 0.2, 0,0,0.7), 3 , byrow = TRUE)
  mu1 = rep(0,3)
  Sigma1 = matrix(c(1,0,0,0,1,0,0,0,1),3,byrow = TRUE)
  gen = mvrnorm(n1,mu = mu1, Sigma  = Sigma1)
  x = gen
  t = c(1:n1)
  phi = 0.9-0.7/(1+exp(beta0*(t-beta1)))
  for(j in c(2:n1)){
    x[j,] = phi[j]*B1%*%x[j-1,] + gen[j,]
  }
  x1 = c(t(x))
  return(x1)
}



R=1000

k1 = 300
k2 = 500
n = 600
M = 30
L = 0

ch = rep(0,R)
S11=rep(0,R)
S21=rep(0,R)
T11=rep(0,R)
T21=rep(0,R)
var_re = rep(0,R)
var_im = rep(0,R)

set.seed(9)


var = function(J,k1,k2,L,M,n) {
  # Input is our time series of N length , M dimensional
  
  xt = VARP(n/3)
  J = fft(xt)/sqrt(n)
  
  ch[r] = c_hatM(J,k1,k2,M)  
  
  v = S_var_c_estP(J,k1,k2-k1,L,M,n)
  vc = T_var_c_conj_estP(J,k1,k2-k1,L,M,n)
  
  S11[r] = v[1]
  S21[r] = v[2]
  T11[r] = vc[1]
  T21[r] = vc[2]
  
  var_re[r] = Re((v[1]+v[2]+vc[1]+vc[2])/2)
  var_im[r] = Re((v[1]+v[2]-vc[1]-vc[2])/2)

  
}


# M= 30, sd=0.0058
# M= 60, sd=0.0029

#L=0  sd = 0.001269747  
#bias = 0.006042356 - 0.004399761 = 0.0016
# MSE = 0.0018

#L = 1  sd = 0.00178791
#bias = 0.006284646 - 0.004935146 = 0.0013
#MSE = 0.0019

#L = 2   sd = 0.00228
# bias =  - 0.00009
# MSE = 0.00176

var(Re(ch))
mean((var_re))
mean(abs(var_re-var(Re(ch))))
sd(var_re)
hist(var_re, main="var(real part)")
abline(v=var(Re(ch)), col="red",lty=2)

var(Im(ch))
mean((var_im))
mean(abs(var_im-var(Im(ch))))
hist(var_im, main="var(Imaginary part)")
abline(v=var(Im(ch)), col="red",lty=2)

z_r = Re(ch)/sqrt(var_re)
sum(abs(z_r)>1.96)


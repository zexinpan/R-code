# Test code for "Dimension-independent convergence rates of randomized nets using median-of-means"
# Paper available at https://arxiv.org/abs/2505.13815

source('auxiliary functions.R')
col = read.table('sobol_Cs.col') # load Sobol' sequence from https://people.cs.kuleuven.be/~dirk.nuyens/qmc-generators/

# Set up the test function used in the paper
gamma=2 # decay rate of relative norm
c_alpha=2/sqrt(5*exp(2)-1) # normalizing constant for ||f||_0,1=1
# c_alpha=2/sqrt(13*exp(2)-5) # normalizing constant for ||f||_1,1=1
f<-function(x){ # input x is n by d
  mean(apply(x,1,function(x){
    prod=1
    for(i in 1:length(x)){
      prod=prod*(1+i^(-gamma)*c_alpha*(x[i]*exp(x[i])-1))
    }
    prod
  }))
}
function_mean=1

# Set up experiment parameters
d=100 # dimension
M=16 # test n=2^1,2^2,...,2^M
r=2*10-1 # number of repeats to take median or mean
precision=64 # number of digits to generate

# Order 1 digital net, mean vs median
errormean=rep(0,M)
errormedian=rep(0,M)
C = array(0,c(M,M,d)) 
for( j in 1:d ){
  for( k in 1:M )
    C[,k,j] = .rsobol.int2bits(col[j,k],M) # load generating matrices
}
for(m in 1:M){
  print(m)
  n=2^m
  I=matrix(0,nrow=m,ncol= n) 
  for(i in 0:(n-1)) I[,i+1]=.rsobol.int2bits(i,m) # store binary expansion of 0,...,2^m-1
  u=matrix(nrow=n,ncol=d)
  temp=rep(0,r)
  L=matrix(0,nrow=precision,ncol=m)
  for(i in 1:m) L[i,i]=1
  for(trial in 1:r){
    for(k in 1:d){
      for(j in 1:m) L[(j+1):precision,j]=sample(c(0,1),precision-j,replace=TRUE) # randomize scrambling matrix
      LC=(L%*%C[1:m,1:m,k])%%2 # randomize generating matrix
      u[,k]=digital_shift((LC%*%I)%%2) # scrambled Sobol' sequence
    }
    temp[trial]=f(u)
  }
  errormean[m]=mean(temp)-function_mean # mean of scrambled nets
  errormedian[m]=median(temp)-function_mean # median of scrambled nets
}
error=errormean
save(error,file="STDalpha0gamma2d100.RData")
error=errormedian
save(error,file="RLSalpha0gamma2d100.RData")

# Median QMC with complete random designs
error=rep(0,M)
for(m in 1:M){
  print(m)
  I=matrix(0,nrow=m,ncol=2^m)
  for(i in 0:(2^m-1)){
    I[,i+1]=.rsobol.int2bits(i,m)  # store binary expansion of 0,...,2^m-1
  }
  u=matrix(nrow=2^m,ncol=d)
  temp=rep(0,r)
  for(trial in 1:r){
    for(k in 1:d) u[,k]=digital_shift((matrix(sample(c(0,1),m*precision,replace=TRUE),nrow = precision,ncol=m)%*%I)%%2) # scrambling under completely random designs
    temp[trial]=f(u)
  }
  error[m]=median(temp)-function_mean
}
save(error,file="CRDalpha0gamma2d100.RData")

# Plotting
y=log(abs(error),2)
x=1:length(y)
lm(y~x)
plot(x,y,ylim=c(-50,0),type="b",xlab="m",ylab="Log2 ABS")
abline(a=0,b=-1)
abline(a=0,b=-1.5)
abline(a=0,b=-2)



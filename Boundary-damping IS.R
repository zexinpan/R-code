# Test code for "Quasi-Monte Carlo integration over R^s with boundary-damping importance sampling"
# Paper available at https://arxiv.org/abs/2509.07509

source('auxiliary functions.R')
col = read.table('sobol_Cs.col') # load Sobol' sequence from https://people.cs.kuleuven.be/~dirk.nuyens/qmc-generators/

# Set up the test function used in the paper
gamma=2 # decay rate of relative norm
M_f=0.3 # growth parameter M in the paper
f<-function(x){
  prod=1
  for(i in 1:length(x)){
    prod=prod*(1+i^(-gamma)*(sqrt(1-2*M_f)*exp(M_f*x[i]^2)-1))
  }
  prod
}
function_mean=1

# A lazy solution to f(x) evaluates to Inf for outliers
# Setting Inf to 0 introduces a small bias negligible for this experiment
# Better replace Inf by the asymptotics if integration error is small
robust_fval<-function(x){
  apply(x,1,function(x){
    ans=f(x)
    if(is.infinite(ans)) return(0)
    ans
  })
}

# Set up experiment parameters
d=5 # dimension
M=18 # test n=2^1,2^2,...,2^M
nexpr=30 # repeat 30 times to estimate the variance
precision=32 # number of digits to generate

# Set up QMC integrator
C = array(0,c(M,M,d)) 
for( j in 1:d ){
  for( k in 1:M )
    C[,k,j] = .rsobol.int2bits(col[j,k],M) # load generating matrices
}
QMCerror<-function(expr="which_experiment"){
  error=matrix(nrow=nexpr,ncol=M)
  for(m in 1:M){
    print(m)
    n=2^m
    if(expr=="truncation") a=sqrt(log(n)/(0.5-M_f)) # Truncation parameter
    I=matrix(0,nrow=m,ncol=n)
    for(i in 0:(n-1)){
      I[,i+1]=.rsobol.int2bits(i,m) # store binary expansion of 0,...,2^m-1
    } 
    u=matrix(nrow=n,ncol=d)
    L=matrix(0,nrow=precision,ncol=m)
    for(i in 1:m) L[i,i]=1
    for(rep in 1:nexpr){
      w=rep(1,n)
      for(k in 1:d){
        for(j in 1:m) L[(j+1):precision,j]=sample(c(0,1),precision-j,replace=TRUE)  # randomize scrambling matrix
        LC=(L%*%C[1:m,1:m,k])%%2 # randomize generating matrix
        u[,k]=digital_shift((LC%*%I)%%2) # scrambled Sobol' sequence
        # Choose transport map and weight according to expr
        if(expr=="IS"){ # boundary-damping Improtance sampling
          w=w*sapply(u[,k],function(x){weightf(x,theta[k])}) # update weights
          u[,k]=sapply(u[,k],function(x){transf(x,theta[k])}) # compute T(u)
        }
        else if(expr=="qnorm"){ # inversion method
          u[,k]=qnorm(u[,k])
        }
        else if(expr=="cot"){ # Mobius transformation
          w=w*weightf(u[,k]) # update weights
          u[,k]=transf(u[,k]) # compute T(u)
        }
        else if(expr=="truncation"){ # Truncation method
          w=w*sapply(u[,k],function(x){weightf(x,a)}) # update weights
          u[,k]=sapply(u[,k],function(x){transf(x,a)}) # compute T(u)
        }
        else return(NA)
      }
      error[rep,m]=mean(robust_fval(u)*w)-function_mean
    }
  }
  return(error)
}

# Boundary-damping improtance sampling
theta=0.1/(1:d)^2 # set theta parameters
transf<-function(x,theta){ # transport map T(u)
  if(x>1/2) u=1-x else u=x
  if(u<theta){
    u=u/theta
    if(u<1/2) p=theta*exp(2-1/u)/8
    else{
      u=1-u
      p=theta*(1/2-u+exp(2-1/u)/8)
    }
  }else{
    p=u-theta/2
  }
  if(x>1/2) return(qnorm(1-p/(1-theta)))
  else return(qnorm(p/(1-theta)))
}
weightf<-function(x,theta){ # weights w(u)
  if(x>1/2) x=1-x
  if(x>theta) ans=1 
  else{
    u=x/theta
    if(u<1/2) ans=exp(2-1/u)/u^2/8
    else{
      u=1-u
      ans=1-exp(2-1/u)/u^2/8
    }
  }
  if(is.nan(ans)) return(0) # division by 0 occurs when evaluating exp(-Inf)
  else return(ans/(1-theta))
}
error=QMCerror(expr="IS")
save(error,file="ISM03gamma2d5.RData")

# Inversion method
error=QMCerror(expr="qnorm")
save(error,file="qnormM03gamma2d5.RData")

# Mobius transformation
transf<-function(x){ # transport map T(u)
  -1/tan(pi*x)
}
weightf<-function(x){ # weights w(u)
  pi/sin(pi*x)^2*exp(-1/2/tan(pi*x)^2)/sqrt(2*pi)
}
error=QMCerror(expr="cot")
save(error,file="MobiusM03gamma2d5.RData")

# Truncation method
transf<-function(x,a){ # transport map T(u)
  a*x-a*(1-x)
}
weightf<-function(x,a){ # weights w(u)
  exp(-(a*x-a*(1-x))^2/2)/sqrt(2*pi)*2*a
}
error=QMCerror(expr="truncation")
save(error,file="TruncationM03gamma2d5.RData")

# Plotting 
MSE=apply(error^2,2,mean)
y=log(MSE,2)/2
x=1:length(y)
lm(y~x)
plot(x,y,ylim=c(-50,0),type="b",xlab="m",ylab="Log2 RMSE")
abline(a=0,b=-0.5)
abline(a=0,b=-1)

# Modified from https://artowen.su.domains/code/rsobol.R

.rsobol.int2bits = function(x,M=32){
  # Convert an integer x into M bits
  # For large x (or small M) you get bits of x modulo 2^M
  # as.vector() is there to force answer to be a vector.
  # This does just one integer (not a vector of them).
  
  ans = rep(0,M)
  for( j in 1:M ){
    ans[j] = as.vector(x %% 2)
    x = (x-ans[j])/2
  }
  ans
}

.rsobol.bits2unif = function(bits){
  # Turn sequence of bits into a point in [0,1)
  # First bits are highest order
  if(is.vector(bits)){
    ans = 0
    for( j in length(bits):1){
      ans = (ans+bits[j])/2
    }
    return(ans)
  }
  
  ans=rep(0,dim(bits)[2])
  for( j in dim(bits)[1]:1){
    ans = (ans+bits[j,])/2
  }
  ans
}

.rsobol.bits2int = function(b){
  # Convert bits b into integers.
  # Inverse of int2bits
  # This is vectorized: each row of the matrix b is a vector of bits
  #
  if( is.vector(b) )
    b = matrix(b,nrow=1)
  
  n = nrow(b)
  p = ncol(b)
  ans = rep(0,n)
  
  for( j in p:1 )
    ans = ans*2+b[,j]
  ans
}


digital_shift<-function(y){
  # Randomization by digital shifts
  #
  E=dim(y)[1]
  for(j in 1:E) if(runif(1)>0.5) y[j,]=1-y[j,]
  .rsobol.bits2unif(y)+runif(1)*2^(-E)
}
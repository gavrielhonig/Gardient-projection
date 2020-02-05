#gavriel honig
#200864734

GRLS = function(f, a, b, tol, maxIter){
  GR = (5^0.5-1)/2
  i = 0
  while ( abs(b-a)>tol & i<maxIter) {
    d =GR *(b-a)
    i = i+1
    if (f(b-d)>f(a+d)) {
      a = b - d}
    else{
      b = a + d
    }
  }
  return((a+b)/2)
}

grd=function(x,y,beta){
  return(-2*t(x)%*%(y-x%*%beta))
}


GP = function(X, Y, b0, maxIter, tol, GRLStol){
  n=0
  loss <- sum((Y-X%*%b0)^2) # f(b0)
  while (n<maxIter) {
  n=n+1
  if(n==maxIter){return(list('b'=b0,'W'=w0,'N'=n,'V'=v))}
  loss <- c(loss,sum((Y-X%*%b0)^2)) # f(b0)
  if (n==1) {w0=which(abs(b0)<tol)}
  if (length(w0)>0) {
  A=matrix(0, nrow = length(w0), ncol = length(b0))
  for (i in 1:nrow(A)) {
    for (j in  1:ncol(A)){
      if (j==w0[i]) {
        A[i,j]=-1  
      }
    }
   }
  P=diag(1,nrow(b0),nrow(b0))-(t(A)%*%solve(A%*%t(A))%*%A)*(dim(A)[1]>0)
  }else{P=diag(1,nrow(b0),nrow(b0))}
  g=grd(X,Y,b0)
  d=-P%*%g
  if(sqrt(t(d)%*%d) > tol){
  a1_demo=c()
  for (i in 1:length(b0)) {
    if (d[i]<0) {
      a1_demo[i]=b0[i]/d[i]
    }else
      {a1_demo[i]=0}
  }
  if (sum(a1_demo==0)==length(a1_demo)){a1=100000000000} 
  else{
  a1=-1*max(a1_demo[a1_demo!=0])} #alpha 1
  idx=which(a1_demo==-a1)
  f<-function(a){return(sum((Y-X%*%(b0+a*d))^2))}
  a2=GRLS(f,0,a1,GRLStol,maxIter) #alpha 2
  if (abs(a1-a2)>tol) {
    b0=b0+a2*d
  }else{
    w0=sort(unique(c(idx,w0)))
    print(paste("Adding constraint", paste(as.character(idx), collapse = ",")))
    b0=b0+a2*d
    }
  }
    if(sqrt(t(d)%*%d) < tol){
    if (length(w0)>0) {
      Mu=-1*solve(A%*%t(A))%*%A%*%g  
    }else{
    Mu=-g}
    v=sum((Y-X%*%b0)^2) # f(b0)
    if (sum(Mu>=-tol)==length(Mu) |n==30000){
      loss=loss[-1]
      index=1:n
      plot(index,loss)
      print(ifelse(length(w0)==0, "No active constraints in the last iteration",
                   paste("Constraints",paste(as.character(c(w0)), collapse = ","),"are active")))
      print(paste("Minimal value", v, "is reached after", n, "updates"))
      
      return(list('b'=b0,'W'=w0,'N'=n,'V'=v))
   }
   min.Mu=min(Mu)
   idx2=which(Mu==min.Mu)
   print(paste("Removing constraint", paste(as.character(w0[idx2]), collapse = ",")))
   w0=sort(w0[-idx2])
  }
 }
} 

x=matrix(c(0.5377,0.3188,3.5784,0.7254,1.8339,-1.3077,2.7694,-0.0631,-2.2588,-0.4336,-1.3499,0.7147,
         0.8622,0.3426,3.0349,-0.2050),4,4)
true_x=matrix(c(1,0,0.5,0),4,1)
y=x%*%true_x
x0=matrix(c(1,1,0,0),4,1)
max_Iter=1000
tol=10^-8
GRLStol=10^-8
solution=GP(x,y,x0,max_Iter,tol,GRLStol)
solution



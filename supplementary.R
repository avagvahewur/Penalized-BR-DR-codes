### Supplementary R-functions
expit <- function(x) exp(x)/(1+exp(x))

# Double robust estimator
U <- function(R,Y,X,gamma,beta){
  n <- length(R)
  int.X<-cbind(rep(1,n),X)
  (R/expit(int.X%*%gamma)*(Y-int.X%*%beta)+(int.X%*%beta))
}

# Double robust estimator for Post-Lasso 
Upost <- function(R,Y,X1,X2,gamma,beta){
  n <- length(R)
  int.X1<-cbind(rep(1,n),X1)
  int.X2<-cbind(rep(1,n),X2)
  (R/expit(int.X1%*%gamma)*(Y-int.X2%*%beta)+(int.X2%*%beta))
}

# Confidence intervals accross replications
colCI <- function(average,se,m1)
{
  p <- dim(average)[2]
  a <- rep(0,dim(average)[2]) 
  for (i in seq(1,p, by=1))
  {
    a[i] <- mean(ifelse((m1<=average[,i]+1.96*se[,i])&(m1>=average[,i]-1.96*se[,i]),1,0),na.rm=TRUE)
  }
  a
}


# CV selection of the penalty term for BR-DR 
crossval<-function(r, x, nf=5, nstep=10, wstart=lam)
{
  region<-wstart+(-nstep/2+3):(nstep/2+2)*wstart/nstep
  N<-length(r)
  winsize<-(N-N%%nf)/nf
  windows<-c(rep(winsize,nf-1),winsize+N%%nf)
  CV<-rep(0,nstep)
  SER<-rep(0,nstep)
  t<-1
  nonc<-rep(0,nstep)
  for (lam in region)
  {
    
    loss<-rep(0,nf)
    for (f in seq(1,nf, by=1))
    {
      first<-1+(f-1)*winsize
      last<-first+windows[f]-1
      
      XX_out<-x[-(first:last),]
      XX_out<-scale(XX_out)
      
      XX_in<-x[(first:last),]
      
      RR_out<-r[-(first:last)]
      RR_in<-r[(first:last)]
      
      int.XX_out<-cbind(rep(1,N-windows[f]),XX_out)
      int.XX_in<-cbind(rep(1,windows[f]),XX_in)
      
      gamma <- Variable(p+1)
      objective <- Minimize( mean((RR_out*exp(-int.XX_out%*%gamma)) + ((1-RR_out)*(int.XX_out%*%gamma)) )+lam*sum(abs(gamma)))
      prob <- Problem(objective)
      sol<-try(result <- solve(prob), silent=TRUE)
      
      if (typeof(sol[1])=="character" | sol$status == "solver_error") {
        
        nonc[t] <- nonc[t] +1
        gamma<-as.vector(rep(0, dim(int.XX_in)[2]))
      } else { gamma<-result$getValue(gamma) }
      
      fit<-int.XX_in%*%gamma
      loss[f]<- mean((RR_in*exp(-fit)) + ((1-RR_in)*(fit)) )
      
    }
    CV[t]<-mean(loss)
    SER[t]<-sqrt(var(loss)/length(loss))
    t<-t+1
  }
  plot(CV, type="l")
  mCV<-min(CV[nonc==0])
  min_position<-(1:nstep)[CV==mCV]
  se_sel<-CV[min_position]+SER[min_position]
  #1se situation
  #se_sel<-CV-se_sel
  #lambda_se<-instep+lstep*(max((1:nstep)[se_sel<0])-1)
  lambda_min<-region[min_position]
  Ret<-matrix(rep(1,1*(nstep+1)), nrow=1)
  Ret[1,1]<-lambda_min
  #Ret[2,1]<-lambda_se
  Ret[1,2:(nstep+1)]<-CV
  #Ret[2,2:(nstep+1)]<-SER
  if(sum(nonc[5:10])!=0) {
    Ret[1,1]<-100
  }
  Ret
}

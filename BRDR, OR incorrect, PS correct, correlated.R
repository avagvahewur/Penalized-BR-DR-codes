install.packages("glmnet")
install.packages("CVXR")

library(glmnet)
library(CVXR)

#Number of replications
nu <- 1000
#Number of observations
n <- 200
#Number of covariates
p <- 100

#Covariance structure
sigma <- diag(p)

for (i in seq(1:p)) 
{
  for(j in seq(1:p))
  { 
    sigma[i,j] <- 0.3^(abs(i-j))
  }
}

e <- eigen(sigma)
sqrtsigma <- e$vectors %*% diag(e$values^0.5) %*% t(e$vectors)

#True target parameter
m1 <- 5

unconv<-rep(0,nu)
lambdas<-rep(0,nu)
exp_problem<-rep(0,nu)
out_problem<-rep(0,nu)
t<-1
result<-matrix(nrow=nu,ncol=8)
se<-matrix(nrow=nu,ncol=8)

for  (i in 1:nu)
{
  set.seed(i)
  
  x<-(matrix(rep(1,p*n),ncol=p))
  for (j in seq(1,p))
  {
    x[,j]<-rnorm(n,0,1)
  }
  ##Correlated vs uncorrelated variables
  x<-t(sqrtsigma%*%t(x))
  
  ##PS model, correct
  r <-rbinom(n,1,expit(t(-0.75*x[,1] + 0.5*x[,2] - 0.25*x[,3] - 0.1*x[,4])))
  ##PS model, incorrect
  #r <-rbinom(n,1,expit(t(-0.75*(x[,1])^2 + 0.5*abs(x[,2])^0.5 - 0.25*(x[,3]) - 0.1*x[,4])))
  
  ##OR model, correct
  #y <-rnorm(n, 5+ 1.5*x[,1] + 1.5*x[,2] + 1.5*x[,3] + 1.5*x[,4] , 2 )
  ##OR model, correct
  y <-rnorm(n, 3.5+ 1.5*x[,1]^2 + 1.5*x[,2]+ 1.5*x[,3] + 1.5*x[,4], 2 )
  
  int.x<-cbind(rep(1,n),x)
  
  #Lasso model for the outcome model
  eps<-0
  repeat{
    cv2 <- cv.glmnet(x[r==1,], y[r==1],alpha=1, nfolds = 5)
    #lam1<-cv2$lambda.min
    lasso <- glmnet(x[r==1,], y[r==1], lambda=cv2$lambda.min, alpha=1)
    lasso <-as.vector(coef(lasso))
    S_outcome<-(1:p)[lasso[-1]!=0]
    if(length(S_outcome)!=0)   {break}
    eps<-eps+cv$lambda.min/10
    out_problem[i]<-1
  }
  
  result[i,1]<-mean(int.x%*%(lasso))
  se[i,1]<-sd(int.x%*%(lasso))/sqrt(n)
  
  #Lasso model for PS model.
  eps<-0
  repeat{
    cv<-cv.glmnet(x,r,family="binomial", nfolds = 5)
    mlerlasso<-glmnet(x,r,family="binomial", lambda=cv$lambda.min-eps, alpha=1)
    mlerlasso<-as.vector(coef(mlerlasso))
    S_exposure<-(1:p)[mlerlasso[-1]!=0]
    
    if(length(S_exposure)!=0)   {break}
    eps<-eps+cv$lambda.min/10
    exp_problem[i]<-1
  }
  
  result[i,2]<-sum(r*y/expit(int.x%*%mlerlasso))/sum((r/expit(int.x%*%mlerlasso)))
  se[i,2]<-sqrt(n)*sd(r*(y-result[i,2])/expit(int.x%*%mlerlasso))/sum((r/expit(int.x%*%mlerlasso)))
  
  #DR estimator with Lasso models
  result[i,3]<-mean(U(r,y,x,mlerlasso,lasso))
  se[i,3]<-sd(U(r,y,x,mlerlasso,lasso))/sqrt(n)
  
  #Post-Selection Lasso
  #Active set of variables in PS model
  S_exposure <- (1:p)[mlerlasso[-1]!=0]
  X_exp <- as.matrix(x[,S_exposure])
  #Active set of variables in OR model
  S_outcome <- (1:p)[lasso[-1]!=0]
  X_out <- as.matrix(x[,S_outcome])
  
  postlasso <- lm(y[r==1]~ X_out[r==1,])
  postlasso <- as.vector(coef(postlasso))
  
  postmlerlasso <- glm(r~X_exp, family="binomial")
  postmlerlasso <- as.vector(coef(postmlerlasso))
  
  #DR estimator with Post-Lasso
  result[i,4] <- mean(Upost(r,y,X_exp, X_out,postmlerlasso, postlasso))
  se[i,4] <- sd(Upost(r,y,X_exp, X_out, postmlerlasso,postlasso))/sqrt(n)
  
  #Double-Selection lasso
  #Joined active set of variables
  S_ds <- sort(union(S_exposure, S_outcome))
  X_ds <- as.matrix(x[,S_ds])
  
  dslasso <- lm(y[r==1] ~ X_ds[r==1,])
  dslasso <- as.vector(coef(dslasso))
  
  dsmlerlasso <- glm(r~ X_ds,family="binomial")
  dsmlerlasso <- as.vector(coef(dsmlerlasso))
  
  result[i,5] <- mean(U(r,y,X_ds,dsmlerlasso, dslasso))
  se[i,5] <- sd(U(r,y,X_ds, dsmlerlasso,dslasso))/sqrt(n)
  
  #  Proposed BR-DR estimator
  # Prespecified penalty parameter, warm starting point
  lam <- (1.1/2)*qnorm(1-0.05/max(n,p*log(n)))/sqrt(n)
  #Selecting the penalty parameter
  lambda <- crossval(r,x,nf=5, nstep=10, wstart=lam)
  if(lambda[1]==100 | is.nan(lambda[1]))  {
    lambda<-lam
  } else {lambda<-lambda[1]}
  
  lambdas[i]<-lambda
  
  # BR-DR estimated PS nuisance parameter
  gamma <- Variable(p+1)
  objective <- Minimize( mean((r*exp(-int.x%*%gamma)) + ((1-r)*(int.x%*%gamma)) )+lambda*sum(abs(gamma)))
  prob <- Problem(objective)
  sol <- try(solve(prob), silent=TRUE)
  
  if (typeof(sol[1])=="character" | is.na(sol$num_iters)) {
    gamma <- Variable(p+1)
    objective <- Minimize( mean((r*exp(-int.x%*%gamma)) + ((1-r)*(int.x%*%gamma)) )+lam*sum(abs(gamma)))
    prob <- Problem(objective)
    sol <- try(solve(prob), silent=TRUE)
    unconv[i]<-1
    if (typeof(sol[1])=="character" | is.na(sol$num_iters)) {
      
      gammabr<-as.vector((mlerlasso))
      unconv[i]<-1
    }     else { gammabr<-as.vector(sol$getValue(gamma)) }
    
  } else { gammabr<-as.vector(sol$getValue(gamma)) }
  
  
  weight<-as.vector(1/exp(gammabr%*%t(int.x)))
  
  cv1 <- cv.glmnet(x[r==1,],y[r==1],weights=weight[r==1], alpha=1, nfolds = 5)
  lassobr <- glmnet(x[r==1,],y[r==1],weights=weight[r==1],lambda=cv1$lambda.min, alpha=1)
  lassobr <- coef(lassobr)
  
  result[i,6] <- mean(int.x%*%((lassobr)))
  se[i,6] <- sd(int.x%*%((lassobr)))/sqrt(n)
  
  result[i,7] <- (sum(r*y/expit(int.x%*%(gammabr))))/(sum(r/expit(int.x%*%(gammabr))))
  se[i,7] <- sqrt(n)*sd(r*(y-result[i,7])/expit(int.x%*%(gammabr)))/(sum(r/expit(int.x%*%(gammabr))))
  
  result[i,8] <- mean(U(r,y,x,gammabr,lassobr))
  se[i,8] <- sd((U(r,y,x,gammabr,lassobr)))/sqrt(n)
  
}


Method <- c("OR LASSO","SB-IPTW LASSO", "DR LASSO",
            "Post LASSO", "DS LASSO", "OR BRDR","SB-IPTW BRDR", "PBRDR LASSO")
Measures <- c("Bias", "RMSE", "MAE", "MCSD", "ASSE", "COV")

Res <- matrix(rep(1,6*8), ncol=6)
Res[,1] <- (colMeans(result)-m1)
Res[,2] <- sqrt(colMeans((result-m1)^2))
Res[,3] <- apply((abs(result-m1)), FUN = median, MARGIN = 2)
Res[,4] <- sqrt(apply(result, FUN=var, MARGIN = 2))
Res[,5] <- colMeans(se)
Res[,6] <- colCI(result,se,m1)
rownames(Res) <- Method
colnames(Res) <- Measures


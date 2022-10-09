method1<-function(xbar,s,n,mu=0,a=0.05){
  p<-length(xbar)
  df2<-n-p
  T2 <- n * t(xbar - mu) %*% solve(s) %*% (xbar - mu)
  Fstat <- T2 / (p * (n-1)) *df2
  pval <- 1 - pf(Fstat, df1=p, df2=df2)
  kritis<-qf(1-a,df1 = p,df2 = df2)
  data.frame(T2=as.numeric(T2), Fstat=as.numeric(Fstat), df1=p, df2=df2,
            T.Kritis=kritis, row.names="")
}

method2<-function(xbar,S,n,mu=0,a=0.05){
  Z2 <- n * t(xbar - mu) %*% solve(S) %*% (xbar - mu)
  p<-length(xbar)
  chisq<-qchisq(1-a,p)
  data.frame(Z2=as.numeric(Z2),df1=p,T.Kritis=chisq,row.names = "")
}


#' Inferences about Mean Vector
#' @author Muhammad Sirojuddin
#' @name One_Population
#' @param xbar xbare
#' @param s matriks cova
#' @param n jumlah sampel
#' @export
onepopulation<-function(xbar,s,n,mu=0,a=0.05,sigma=c("known","unknown")){
  if(sigma=="known") method2(xbar,s,n,mu,a)
  if(sigma=="unknown") method1(xbar,s,n,mu,a)
}


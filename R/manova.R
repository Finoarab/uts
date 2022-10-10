xbartot<-function(xbar,n){
  p<-nrow(xbar)
  g<-ncol(xbar)
  xbartot<-matrix(0,p,1)
  for(i in 1:g){
    xbartot<-xbartot+xbar[,i]*n[i]
  }
  return(xbartot/sum(n))
}

lamda<-function(xbar,s,n,b=0,w=0){

  xbartot<-xbartot(xbar,n)
  for(l in 1:ncol(xbar)){
    b <- b + (n[l] *(xbar[,l]-xbartot) %*% t(xbar[,l]))
    w <- w + ((n[l]-1) * s[[l]])
  }
  lamda<-det(w)/det(w+b)
  return(list(b=b,w=w,lamda=lamda))
}

onesim<-function(xbar,s,n,a=0.05,method=c("kecil","besar")){
  g<-ncol(xbar)
  p<-nrow(xbar)
  xbartot<-xbartot(xbar,n)
  wbl<-lamda(xbar,s,n)
  lamda<-wbl[['lamda']]

  fuji<-(sum(n)-1-(p+g)/2)*-log(lamda)
  kritis<-qchisq(1-a,p*(g-1))
  if(method=='kecil'){
    fuji<-(sum(n)-p-2)/p*(1-sqrt(lamda))/sqrt(lamda)
    kritis<-qf(1-a,2*p,2*sum(n)-p-2)
  }
  return(list(fuji=fuji,
              kritis=kritis))
}


oneci<-function(xbar,s,n,a=0.05,j){
  g<-ncol(xbar)
  p<-nrow(xbar)
  t<-qt((1-a)/(p*g*(g-1)),sum(n[j])-g)
  w<-lamda(xbar,s,n)[['w']]
  tmp<-data.frame(rbind(1,2))
  xbar1 <- xbar[,j[1]]
  xbar2 <- xbar[,j[2]]
  for(i in 1:p){
  lower<-xbar1[i]-xbar2[i]+t*sqrt((1/n[j[1]]+1/n[j[2]])*w[i,i]/(sum(n[j])-g))
  upper<-xbar1[i]-xbar2[i]-t*sqrt((1/n[j[1]]+1/n[j[2]])*w[i,i]/(sum(n[j])-g))
  tmp[,i]<-rbind(lower,upper)
  }
  return(tmp)
}

#' One way MANOVA
#' @author Muhammad Sirojuddin
#' @name OnewayManova
#' @param xbar dadekno siji px2
#' @param s dadekno list(s1,s2)
#' @param n c(n1,n2)
#' @param method c("kecil","besar") ukuran sample
#' ci[[1]] menunjukkan ci antara grup 1 dan 2
#' @export
oneway<-function(xbar,s,n,a=0.05,method=c("kecil",'besar')){
  g<-ncol(xbar)
  xbartot<-xbartot(xbar,n)
  wbl<-lamda(xbar,s,n)
  simtest<-onesim(xbar,s,n,a,method = method)
  ci<-list()
  for(j in 1:g-1)
    x<-oneci(xbar,s,n,a,j=c(j,j+1))
  ci[[j]]<-x
  return(list(
    xbarTotal=xbartot,
    wbl=wbl,
    simtest=simtest,
    ci=ci
  ))
}



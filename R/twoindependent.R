spooled<-function(s,n){
  spooled<-((n[1]-1)*s[[1]]+(n[2]-1)*s[[2]])/(sum(n)-2)
  return(spooled)
}

t2<-function(xbar,s,n,mu=0){
  si<-(sum(n)/(n[1]*n[2]))*spooled(s,n)
  xdiff<-xbar[,1]-xbar[,2]
  t2<-t(xdiff)%*%solve(si)%*%xdiff
  return(t2)
}

c2<-function(xbar,n,a=0.05){
  p<-nrow(xbar)
  f<-qf(1-a,p,sum(n)-p-1)
  c2<-(sum(n)-2)*p/(sum(n)-p-1)*f
  return(c2)
}

ci<-function(xbar,s,n,a=0.05){
  ci<-data.frame(rbind("lower","upper"))
  p<-nrow(xbar)
  xdiff<-xbar[,1]-xbar[,2]
  si<-spooled(s,n)
  t<-qt(1-a/(2*p),df = sum(n)-2)
  for(i in 1:p){
    lower<-xdiff[i]-t*sqrt((1/n[1]+1/n[2])*si[i,i])
    upper<-xdiff[i]+t*sqrt((1/n[1]+1/n[2])*si[i,i])
    ci[,i]<-rbind(lower,upper)
  }
  return(ci)
}

#' Two Independent Population
#' @author Muhammad Sirojuddin
#' @name Two_Independent_Population
#' @param xbar dadekno siji px2
#' @param s dadekno list(s1,s2)
#' @param n c(n1,n2)  
#' @export
twoindependent<-function(xbar,s,n,a=0.05,mu=0){
  result<-list(
    "spooled"=spooled(s,n),
    "t2"=t2(xbar,s,n,mu),
    "c2"=c2(xbar,n,a),
    "ci"=ci(xbar,s,n,a)
  )
  return(result)
  }

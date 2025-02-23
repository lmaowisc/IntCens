profile1=function(Hf,beta,yini,delta,gamma,index,Z,eps=1e-4,maxiter=100){

  y=yini
  yp=0

  i=1
  while(sum(abs(y-yp))>eps&i<=maxiter){
    yp=y
    csd=CSD.f2(Hf=Hf,beta=beta,y=y,index=index,delta=delta,gamma=gamma,Z=Z)
    G=csd$G
    Q=csd$Q
    y=CM(G=G,Q=Q)
    i=i+1
  }
  if (i==(maxiter+1)){
   # print("Non-convergence!y\n")
  }
  # print(i)
  return(y)
}

NR2=function(Hf,beta,y,delta,gamma,index,Z,eps=1e-4,maxiter=20){

  n=length(delta)
  p=length(beta)
  H=Hf[[1]]
  H1=Hf[[2]]
  H2=Hf[[3]]


  betap=rep(0,length(beta))-1
  k=1
  while(sum(abs(beta-betap))>eps&k<=maxiter){
    score=rep(0,p)
    info=matrix(0,p,p)
    epn=exp(Z%*%beta)

    for (i in 1:n){
      if (delta[i]==1){
        j=index[i,1]
        x=y[j]
        Fd=FdG.f(epn[i],x,H)
        if (Fd<=1e-10){Fd=1e-10}


        xi1=xi1G.f(epn[i],Z[i,],x,H1)
        xi2=xi2G.f(epn[i],Z[i,],x,H1,H2)


        score=score+xi1/Fd
        info=info+Fd^(-2)*xi1%*%t(xi1)-xi2/Fd
      }
      else{
        if (gamma[i]==1){
          j1=index[i,1]
          j2=index[i,2]
          x1=y[j1]
          x2=y[j2]
          Fd1=FdG.f(epn[i],x1,H)
          Fd2=FdG.f(epn[i],x2,H)
          if (Fd2<=Fd1+1e-10){Fd2=Fd1+1e-10}

          xi11=xi1G.f(epn[i],Z[i,],x1,H1)
          xi12=xi1G.f(epn[i],Z[i,],x2,H1)

          xi21=xi2G.f(epn[i],Z[i,],x1,H1,H2)
          xi22=xi2G.f(epn[i],Z[i,],x2,H1,H2)


          score=score+(Fd2-Fd1)^(-1)*(xi12-xi11)

          info=info+(Fd2-Fd1)^(-2)*(xi12-xi11)%*%t(xi12-xi11)-(Fd2-Fd1)^(-1)*(xi22-xi21)
        }
        else{
          j=index[i,2]
          x=y[j]
          Fd=FdG.f(epn[i],x,H)

          xi1=xi1G.f(epn[i],Z[i,],x,H1)
          xi2=xi2G.f(epn[i],Z[i,],x,H1,H2)

          #           if (sum(V[i]>=tt-1e-10)==0){
          #             ext=1
          #           }
          #           else{
          #             epnt=exp(Z[i]*betat)
          #             ext=Ht(epnt*yt[sum(V[i]>=tt-1e-10)])
          #             #cat(ext,Fd)
          #           }
          #           if (is.na(ext)){
          #             print("is na")
          #           }
          if (Fd>=1-1e-10){Fd=1-1e-10}

          score=score-xi1/(1-Fd)

          info=info+(1-Fd)^(-2)*xi1%*%t(xi1)+(1-Fd)^(-1)*xi2
          #print(zeta3)
        }
      }
    }
    betap=beta
    beta=beta+solve(info)%*%score
    k=k+1
  }
  return(beta)
}

lik.f2=function(H,beta,y,delta,gamma,index,Z){

  n=length(delta)
  lik=0
  epn=exp(Z%*%beta)

  for (i in 1:n){
    if (delta[i]==1){
      j=index[i,1]
      x=y[j]
      Fd=FdG.f(epn[i],x,H)
      if (Fd<=1e-6){Fd=1e-6}
      lik=lik+log(Fd)
    }
    else{
      if (gamma[i]==1){
        j1=index[i,1]
        j2=index[i,2]
        x1=y[j1]
        x2=y[j2]
        Fd1=FdG.f(epn[i],x1,H)
        Fd2=FdG.f(epn[i],x2,H)
        if (Fd2<=Fd1+1e-6){Fd2=Fd1+1e-6}
        lik=lik+log(Fd2-Fd1)


      }

      else{
        j=index[i,2]
        x=y[j]
        Fd=FdG.f(epn[i],x,H)
        if (Fd>=1){Fd=ext-1e-10}
        lik=lik+log(1-Fd)
        #print(zeta3)
      }
    }
  }

  return(lik)
}


profile2=function(Hf,y,beta,delta,gamma,index,Z,eps=1e-4,maxiter=100){

  yp=0
  j=1
  while(sum(abs(y-yp))>eps&j<=maxiter){
    yp=y

    y=profile1(Hf=Hf,beta=beta,yini=y,delta=delta,gamma=gamma,index=index,Z=Z,eps=1e-3,maxiter=100)
    j=j+1
  }

  lik=lik.f2(Hf[[1]],beta,y,delta,gamma,index,Z)
  return(lik)

}


quad.profile2=function(n,Hf,beta,y,delta,gamma,index,Z,eps=1e-4,maxiter=100){
  p=length(beta)

  ptb=1/sqrt(n)

  #likelihood at beta
  y=profile1(Hf=Hf,y=y,beta=beta,delta=delta,gamma=gamma,index=index,Z=Z,eps=eps,maxiter=maxiter)
  lik0=lik.f2(Hf[[1]],beta,y,delta,gamma,index,Z)

  info=matrix(NA,p,p)

  for (i in 1:p){
    beta.tmp=beta
    beta.tmp[i]=beta.tmp[i]+ptb
    y.tmp=profile1(Hf=Hf,y=y,beta=beta.tmp,delta=delta,gamma=gamma,index=index,Z=Z,eps=eps,maxiter=maxiter)
    lik.tmp=lik.f2(Hf[[1]],beta.tmp,y.tmp,delta,gamma,index,Z)
    info[i,i]=-2*(lik.tmp-lik0)/ptb^2
    #cat(i,"info[i,i]=",info[i,i],"\n")
  }
  if (i>1){
    for (i in 2:p){
      for (j in 1:(i-1)){
        beta.tmp=beta
        beta.tmp[i]=beta.tmp[i]+ptb
        beta.tmp[j]=beta.tmp[j]+ptb
        y.tmp=profile1(Hf=Hf,y=y,beta=beta.tmp,delta=delta,gamma=gamma,index=index,Z=Z,eps=eps,maxiter=maxiter)
        lik.tmp=lik.f2(Hf[[1]],beta.tmp,y.tmp,delta,gamma,index,Z)
        val=(-2*(lik.tmp-lik0)/ptb^2-info[i,i]-info[j,j])/2
        info[i,j]=val
        info[j,i]=val
        #cat("info[i,j]=",info[i,j],"\n")
      }

    }
  }
  return(info)

}


#' @export
ICSurvICM=function(delta,gamma,U,V,Z,model,maxiter=500){


  if (model!="PO"&model!="PH"){

    G.l=trans_log(0)
    Hf=H.func(G.l[[1]],G.l[[2]],G.l[[3]])

    #dummy data
    n=length(delta)
    Z=as.matrix(rep(0,n))

    obj=order_rank(delta,gamma,U,V,Z)

    data=obj$data
    t=obj$t
    index=obj$index
    delta=data[,1]
    gamma=data[,2]
    U=data[,3]
    V=data[,4]
    Z=as.matrix(data[,5:ncol(data)])


    mp=length(t)
    y=(1:mp)/mp



    yp=0


    beta=rep(0,ncol(Z))
    betap=beta-1

    eps=1e-3

    j=1

    while(sum(abs(y-yp))>eps&j<=maxiter){
      yp=y

      y=profile1(Hf=Hf,beta=beta,yini=y,delta=delta,gamma=gamma,index=index,Z=Z,eps=1e-3,maxiter=200)

      j=j+1
    }
    if (j==maxiter+1){
      print(conv=F)
    }else{
      conv=T
      var=NULL
    }
  }else{
  if (model=="PO"){
    r=1
  }else{
    r=0
    }

  G.l=trans_log(r)
  Hf=H.func(G.l[[1]],G.l[[2]],G.l[[3]])

#   delta=data[,1]
#   gamma=data[,2]
#   U=data[,3]
#   V=data[,4]
#   Z=as.matrix(data[,5:ncol(data)])

  Z=as.matrix(Z)
  obj=order_rank(delta,gamma,U,V,Z)

  data=obj$data
  t=obj$t
  index=obj$index
  delta=data[,1]
  gamma=data[,2]
  U=data[,3]
  V=data[,4]
  Z=as.matrix(data[,5:ncol(data)])
  n=length(delta)

  mp=length(t)
  y=(1:mp)/mp



  yp=0


  beta=rep(0,ncol(Z))
  betap=beta-1

  eps=1e-4

  j=1

  while(sum(abs(beta-betap))>eps&j<=maxiter){
    yp=y

    y=profile1(Hf=Hf,beta=beta,yini=y,delta=delta,gamma=gamma,index=index,Z=Z,eps=1e-3,maxiter=200)
    betap=beta
    beta=NR2(Hf=Hf,beta=beta,y=y,delta=delta,gamma=gamma,index=index,Z=Z,eps=1e-4,maxiter=20)
    j=j+1
  }
  if (j==maxiter+1){
    print("conv=F")
  }else{
    conv=T
    var=solve(quad.profile2(n,Hf,beta,y,delta,gamma,index,Z,eps=1e-3,maxiter=maxiter))
  }
  }
    obj=list(beta=beta,y=y,t=t,var=var,j=j,model=model,conv=conv,varnames=colnames(Z),call=match.call())
    class(obj)="ICSurv"
    return(obj)
}

#' @export
icsurfit <- function(L, R, Z, model = "NP", maxiter = 500){

  # delta, gamma, U, V
  # EXAMPLES
  # L <- c(NA, 1, - Inf, 0, 2, 1)
  # R <- c(2, Inf, 5, 2, 4, NA)

  # Convert L & R
  # L: NA, <=0, - Inf to 0
  # R: NA to Inf
  L[is.na(L)] <- 0
  L[L <= 0] <- 0
  L[L == -Inf] <- 0
  R[is.na(R)] <- Inf

  # check L<= R
  if (any(L > R)){
    stop("L must be less than or equal to R")
  }

  # Censoring type: left, interval, right
  # 1: left, 2: interval, 3: right
  flag <- ifelse(L == 0, 1, ifelse(L > 0 & R < Inf, 2, 3))
  # delta: left censoring indicator
  delta <- (flag == 1) + 0
  # gamma: int censoring indicator
  gamma <- (flag == 2) + 0
  # U: R if left censoring, L if interval censoring, NA if right censoring
  U <- ifelse(flag == 3, NA, ifelse(flag == 2, L, R))
  # V: NA if left censoring, R if interval censoring, L if right censoring
  V <- ifelse(flag == 1, NA, ifelse(flag == 2, R, L))


  if (model!="PO"& model!="PH"){

    G.l=trans_log(0)
    Hf=H.func(G.l[[1]],G.l[[2]],G.l[[3]])

    #dummy data
    n=length(delta)
    Z=as.matrix(rep(0,n))

    obj=order_rank(delta,gamma,U,V,Z)

    data=obj$data
    t=obj$t
    index=obj$index
    delta=data[,1]
    gamma=data[,2]
    U=data[,3]
    V=data[,4]
    Z=as.matrix(data[,5:ncol(data)])


    mp=length(t)
    y=(1:mp)/mp



    yp=0


    beta=rep(0,ncol(Z))
    betap=beta-1

    eps=1e-3

    j=1

    while(sum(abs(y-yp))>eps&j<=maxiter){
      yp=y

      y=profile1(Hf=Hf,beta=beta,yini=y,delta=delta,gamma=gamma,index=index,Z=Z,eps=1e-3,maxiter=200)

      j=j+1
    }
    if (j==maxiter+1){
      print(conv=F)
    }else{
      conv=T
      var=NULL
    }
  }else{
    if (model=="PO"){
      r=1
    }else{
      r=0
    }

    G.l=trans_log(r)
    Hf=H.func(G.l[[1]],G.l[[2]],G.l[[3]])

    #   delta=data[,1]
    #   gamma=data[,2]
    #   U=data[,3]
    #   V=data[,4]
    #   Z=as.matrix(data[,5:ncol(data)])

    Z=as.matrix(Z)
    obj=order_rank(delta,gamma,U,V,Z)

    data=obj$data
    t=obj$t
    index=obj$index
    delta=data[,1]
    gamma=data[,2]
    U=data[,3]
    V=data[,4]
    Z=as.matrix(data[,5:ncol(data)])
    n=length(delta)

    mp=length(t)
    y=(1:mp)/mp



    yp=0


    beta=rep(0,ncol(Z))
    betap=beta-1

    eps=1e-4

    j=1

    while(sum(abs(beta-betap))>eps&j<=maxiter){
      yp=y

      y=profile1(Hf=Hf,beta=beta,yini=y,delta=delta,gamma=gamma,index=index,Z=Z,eps=1e-3,maxiter=200)
      betap=beta
      beta=NR2(Hf=Hf,beta=beta,y=y,delta=delta,gamma=gamma,index=index,Z=Z,eps=1e-4,maxiter=20)
      j=j+1
    }
    if (j==maxiter+1){
      print("conv=F")
    }else{
      conv=T
      var=solve(quad.profile2(n,Hf,beta,y,delta,gamma,index,Z,eps=1e-3,maxiter=maxiter))
    }
  }
  obj=list(beta=beta,y=y,t=t,var=var,j=j,model=model,conv=conv,varnames=colnames(Z),call=match.call())
  class(obj)="icsurvfit"
  return(obj)
}

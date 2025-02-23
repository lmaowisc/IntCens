
#### Functions to compute the log-likelihood and the CSD#######
#cdf
FdG.f=function(epn,y,H){
  return(1-H(epn*y))
}

eta1G.f=function(epn,y,H1){
  return(-epn*H1(epn*y))
}

eta2G.f=function(epn,y,H2){
  return(-epn*epn*H2(epn*y))
}

zeta1.f=function(Fd,eta1,eta2){
  return(Fd^(-2)*eta1^2-Fd^(-1)*eta2)
}

zeta21.f=function(Fd1,Fd2,eta1,eta2){
  return((Fd2-Fd1)^(-2)*eta1^2+(Fd2-Fd1)^(-1)*eta2)
}
zeta22.f=function(Fd1,Fd2,eta1,eta2){
  return((Fd2-Fd1)^(-2)*eta1^2-(Fd2-Fd1)^(-1)*eta2)
}
zeta3CR.f=function(Fd,eta1,eta2,ext){
  return((ext-Fd)^(-2)*eta1^2+(ext-Fd)^(-1)*eta2)
}


#derivative of beta
xi1G.f=function(epn,z,y,H1){
  return(-z*epn*y*H1(epn*y))
}

xi2G.f=function(epn,z,y,H1,H2){
  return(-z%*%t(z)*epn*y*(H1(epn*y)+y*epn*H2(epn*y)))
}

#############################################################

# Transformation functions

trans_log=function(gamma){

  transform=list()
  if (gamma==0){
    transform[[1]]=function(x){
      return(x)
    }
    transform[[2]]=function(x){
      return(1)
    }
    transform[[3]]=function(x){
      return(0)
    }
    transform[[4]]=transform[[3]]
    transform[[5]]=transform[[3]]
  }else{
    transform[[1]]=function(x){
      return(1/gamma*log(1+gamma*x))
    }
    transform[[2]]=function(x){
      return(1/(1+gamma*x))
    }
    transform[[3]]=function(x){
      return(-gamma/(1+gamma*x)^2)
    }

    transform[[4]]=function(x){
      return(2*gamma^2/(1+gamma*x)^3)
    }

    transform[[5]]=function(x){
      return(-6*gamma^3/(1+gamma*x)^4)
    }
  }
  return (transform);

}



trans_inv=function(gamma){

  if (gamma==0){

    Ginv=function(y){
      return (y)
    }
  }else{
    Ginv=function(y){
      return (1/gamma*(exp(gamma*y)-1))
    }

  }
  return (Ginv)

}


H.func=function(G,G1,G2){
  H=function(x){
    return(exp(-G(x)))
  }

  H1=function(x){
    return(-H(x)*G1(x))
  }

  H2=function(x){
    return(-H1(x)*G1(x)-H(x)*G2(x))
  }

  return(list(H=H,H1=H1,H2=H2))
}

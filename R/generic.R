
##generic functions
#' @export
print.icsurvfit=function(object,...){
cat("\n")
cat("Call:\n")
print(object$call)
cat("\n")
if (object$conv){

  beta=object$beta
  var=object$var
  y=object$y
  t=object$t
  model=object$model
  varnames=object$varnames

  if (model!="PO"&model!="PH"){
     result=cbind(t,exp(-y))
    rownames(result)=1:length(y)
    colnames(result)=c("Time","Survival rate")

    cat("NPMLE for interval-censored data:\n")
    cat("ICM algorithm converges in", object$j,"iterations.\n\n")

    cat("Estimated survival function\n")
    cat("---------------------------\n")
    print(result)
    cat("\n")
  }else{

    p=length(beta)
    se=sqrt(diag(var))
    pv=2*(1-pnorm(abs(beta/se)))

    if (model=="PO"){
      mod="proportional odds"
    }else{
      mod="proportional hazards"
    }


cat("NPMLE of",mod,"for interval-censored data:\n\n")
cat("ICM algorithm converges in", object$j,"iterations.\n\n")


table=cbind(Estimate = as.vector(beta),
             StdErr = se,
             z.value = as.vector(beta/se),
             p.value = as.vector(pv))

rownames(table)=varnames

cat("Maximum Likelihood Estimates for Regression parameters:\n\n")
printCoefmat(table, P.value=TRUE, has.Pvalue=TRUE)
cat("\n")
}
}else{
  cat("Convergence failure.\n")
}
}



#' @export
plot.icsurvfit=function(object,z=NULL,xlab="Time",
ylab="Survival rate",lty=1,frame.plot=F,add=F,
 ylim=c(0,1),...){

  beta=object$beta
  var=object$var
  y=object$y
  t=object$t
  model=object$model
  varnames=object$varnames

  if (model!="PO"&model!="PH"){

    if (add==F){
    plot(stepfun(t,c(1,exp(-y))),do.point=F,ylim=ylim,
         lty=lty,xlab=xlab,ylab=ylab,frame.plot=frame.plot,...)
    }else{
      lines(stepfun(t,c(1,exp(-y))),do.point=F,lty=lty,...)
    }
  }else{


    if (model=="PO"){
      if (add==F){
      plot(stepfun(t,c(1,(1+exp(sum(beta*z))*y)^(-1))),do.point=F,ylim=ylim,
           lty=lty,xlab=xlab,ylab=ylab,frame.plot=frame.plot,...)
    }else{
      lines(stepfun(t,c(1,(1+exp(sum(beta*z))*y)^{-1})),do.point=F,lty=lty,...)
    }


    }else{
      if (add==F){
        plot(stepfun(t,c(1,exp(-exp(sum(beta*z))*y))),do.point=F,ylim=ylim,
             lty=lty,xlab=xlab,ylab=ylab,frame.plot=frame.plot,...)
      }else{
        lines(stepfun(t,c(1,exp(-exp(sum(beta*z))*y))),do.point=F,lty=lty,...)
      }
    }

  }
}




#' Penalized weighted least-squares estimate for variable selection on correlated multiply imputed data
#'
#'This is a functions to estimate coefficients of wighted leat-squares model and
#'select variables for multiple imputed data sets
#',considering the correlation of multiple imputed observations.
#'
#'@param missdata A Matrix,missing data with variables X in the first p columns and response Y at the last column.
#'@param mice_time An intedevger, number of imputation.
#'@param penalty The method for variable selection,choose from "lasso" or "alasso".
#'@param lamda.vec Optimal tuning parameter for penalty,default seq(1,4,length.out=12).
#'@param Gamma Parameter for adjustment of the Adaptive Weights vector in adaptive LASSO,default c(0.5,1,1.5).
#'@return A Vsmi_est object, contians estcoef and index_sig , estcoef for estimate coefficients and index_sig for selected variable index.
#' @export
#'
#'@examples \donttest{
#'library(MASS)
#'library(mice)
#'library(qif)
#'entire<-generate_pwls_missing_data()
#'est_lasso<-PWLS(entire,penalty="lasso")
#'est_alasso <- PWLS(entire,penalty = "alasso")
#'}
PWLS <- function(missdata,mice_time=5,penalty="alasso",lamda.vec=seq(6,24,length.out=40),Gamma=c(0.5,1,2)){


  ##inner functions

  GEEsym.estimate=function(mydata,yuangeshu,D,b,maxiter=200,eps=1e-3){

    ## for each mydata, the first p columns are covariates X, and the last two is the outcome Y, and the last one is the index D
    rr<-NULL
    ## number of observations
    n=dim(mydata)[1]

    ## number of covariates
    p = dim(mydata)[2]-1

    ## Standardize covariates X and center outcome Y

    x = mydata[,1:p]
    y = mydata[,(p+1)]

    meany = mean(y)
    y= y - meany
    y<-matrix(y,nrow=n)
    #

    ## Rho moment estimation
    r=y-x%*%b

    r_time=n/yuangeshu

    r_list=list()
    r_string=""
    for (r_i in 1:r_time){
      exp11_chr=paste("r",as.character(r_i),"=rep(0,yuangeshu)",sep="")
      eval(parse(text=exp11_chr))
      r_name=paste("r",as.character(r_i),sep="")
      r_list=append(r_list,r_name)
      r_string=paste(r_string,r_name,sep=",")
    }

    for (r_i in 1:(r_time-1)){
      exp11_chr=paste("for(i in 1:yuangeshu){r",as.character(r_i),"[i]=r[(i-1)*D+",as.character(r_i),"]}",sep="")
      eval(parse(text=exp11_chr))}
    exp11_chr=paste("for(i in 1:yuangeshu){r",as.character(r_time),"[i]=r[i*D]}",sep="")
    eval(parse(text=exp11_chr))

    cbind_exp_chr=paste("rr<-cbind(",substring(r_string,2),")",sep="")
    eval(parse(text=cbind_exp_chr))

    rho=0
    for(j in 1:5){
      for(k in 1:5){
        if(j!=k){
          rho=rho+cor(rr[,j],rr[,k])
        }
      }
    }

    rho_sum=rho/(D*(D-1))
    ## working matrix
    #library(Matrix)
    V0=rho_sum*matrix(1,D,D)+(1-rho_sum)*diag(D)

    ##############GEE##############
    iter=0
    dif=1
    while(iter<=maxiter & dif>=eps){
      iter=iter+1

      ## beta
      b.old=b

      sum_xx=0
      sum_yy=0


      V0_ni=solve(V0)

      for(i in 1:yuangeshu){
        xx<-t(x[((i-1)*D+1):(i*D),])%*%V0_ni%*%x[((i-1)*D+1):(i*D),]
        sum_xx<-sum_xx+xx
      }
      for(i in 1:yuangeshu){
        yy<-t(x[((i-1)*D+1):(i*D),])%*%V0_ni%*%(y[((i-1)*D+1):(i*D),]-x[((i-1)*D+1):(i*D),]%*%b.old)
        sum_yy<-sum_yy+yy
      }



      xx_ni<-solve(sum_xx)

      b=b.old+xx_ni%*%sum_yy

      ## Rho moment estimation
      r=y-x%*%b

      r_time=n/yuangeshu

      r_list=list()
      r_string=""
      for (r_i in 1:r_time){
        exp11_chr=paste("r",as.character(r_i),"=rep(0,yuangeshu)",sep="")
        eval(parse(text=exp11_chr))
        r_name=paste("r",as.character(r_i),sep="")
        r_list=append(r_list,r_name)
        r_string=paste(r_string,r_name,sep=",")
      }

      for (r_i in 1:(r_time-1)){
        exp11_chr=paste("for(i in 1:yuangeshu){r",as.character(r_i),"[i]=r[(i-1)*D+",as.character(r_i),"]}",sep="")
        eval(parse(text=exp11_chr))}
      exp11_chr=paste("for(i in 1:yuangeshu){r",as.character(r_time),"[i]=r[i*D]}",sep="")
      eval(parse(text=exp11_chr))

      cbind_exp_chr=paste("rr<-cbind(",substring(r_string,2),")",sep="")
      eval(parse(text=cbind_exp_chr))

      rho=0
      for(j in 1:5){
        for(k in 1:5){
          if(j!=k){
            rho=rho+cor(rr[,j],rr[,k])
          }
        }
      }

      rho_sum=rho/(D*(D-1))

      ## working matrix
      V0=rho_sum*matrix(1,D,D)+(1-rho_sum)*diag(D)

      dif = max(abs(b-b.old))


    }

    return(list(beta=b,V=V0))

  }

  GEEsym_gl.lasso.estimate=function(mydata,yuangeshu,D,b,V0,lamda,maxiter=200,eps=1e-3){

    ## for each mydata, the first p columns are covariates X, and the last two is the outcome Y, and the last one is the index D
    rr<-NULL
    n=dim(mydata)[1]

    ## number of covariates
    p = dim(mydata)[2]-1

    ## Standardize covariates X and center outcome Y

    x = mydata[,1:p]
    y = mydata[,(p+1)]

    meany = mean(y)
    y= y - meany
    y<-matrix(y,nrow=n)

    ##############GEE##############


    iter=0
    dif=1

    V0_ni=solve(V0)


    sum_sigma=0
    for(i in 1:yuangeshu){
      sigma=t(y[((i-1)*D+1):(i*D),]-x[((i-1)*D+1):(i*D),]%*%b)%*%V0_ni%*%(y[((i-1)*D+1):(i*D),]-x[((i-1)*D+1):(i*D),]%*%b)/(n-p)
      sum_sigma=sum_sigma+sigma
    }
    sum_sigma<-c(sum_sigma)


    cc<-1/sqrt(b^2)
    cc<-diag(c(cc))

    while(iter<=maxiter & dif>=eps){
      iter=iter+1

      ## beta
      b.old=b

      sum_xx=0
      sum_yy=0


      V0_ni=solve(V0)

      for(i in 1:yuangeshu){
        xx<-t(x[((i-1)*D+1):(i*D),])%*%V0_ni%*%x[((i-1)*D+1):(i*D),]/sum_sigma###HAO XIANG YOU WEN TI
        sum_xx<-sum_xx+xx
      }
      sum_xx<-sum_xx++2*lamda*cc
      for(i in 1:yuangeshu){
        yy<-t(x[((i-1)*D+1):(i*D),])%*%V0_ni%*%(y[((i-1)*D+1):(i*D),]-x[((i-1)*D+1):(i*D),]%*%b.old)/sum_sigma###HAO XIANG YOU WEN TI
        sum_yy<-sum_yy+yy
      }
      sum_yy<-sum_yy-2*lamda*cc%*%b.old



      xx_ni<-solve(sum_xx)


      b=b.old+xx_ni%*%sum_yy

      ## Rho moment estimation
      r=y-x%*%b

      r_time=n/yuangeshu

      r_list=list()
      r_string=""
      for (r_i in 1:r_time){
        exp11_chr=paste("r",as.character(r_i),"=rep(0,yuangeshu)",sep="")
        eval(parse(text=exp11_chr))
        r_name=paste("r",as.character(r_i),sep="")
        r_list=append(r_list,r_name)
        r_string=paste(r_string,r_name,sep=",")
      }

      for (r_i in 1:(r_time-1)){
        exp11_chr=paste("for(i in 1:yuangeshu){r",as.character(r_i),"[i]=r[(i-1)*D+",as.character(r_i),"]}",sep="")
        eval(parse(text=exp11_chr))}
      exp11_chr=paste("for(i in 1:yuangeshu){r",as.character(r_time),"[i]=r[i*D]}",sep="")
      eval(parse(text=exp11_chr))

      cbind_exp_chr=paste("rr<-cbind(",substring(r_string,2),")",sep="")
      eval(parse(text=cbind_exp_chr))

      rho=0
      for(j in 1:5){
        for(k in 1:5){
          if(j!=k){
            rho=rho+cor(rr[,j],rr[,k])
          }
        }
      }

      rho_sum=rho/(D*(D-1))
      ## working matrix

      V0=rho_sum*matrix(1,D,D)+(1-rho_sum)*diag(D)


      ## sigma
      sum_sigma=0


      V0_ni=solve(V0)


      for(i in 1:yuangeshu){
        sigma=t(y[((i-1)*D+1):(i*D),]-x[((i-1)*D+1):(i*D),]%*%b)%*%V0_ni%*%(y[((i-1)*D+1):(i*D),]-x[((i-1)*D+1):(i*D),]%*%b)/(n-p)
        sum_sigma=sum_sigma+sigma
      }
      sum_sigma<-c(sum_sigma)

      ## CC
      cc<-1/sqrt(b^2)
      cc<-c(cc)
      cc[cc>1/(1e-10)] =1/(1e-10) ###kaolv shidangde fangkuai zuixiao xianzhi
      cc<-diag(cc)


      dif = max(abs(b-b.old))


    }
    b[abs(b)<=1e-4,] = 0

    return(list(beta=b,V=V0))
  }

  GEEsym_gl.alasso.estimate=function(mydata,yuangeshu,D,b,V0,lamda,gamma,maxiter=200,eps=1e-3){

    ## for each mydata, the first p columns are covariates X, and the last two is the outcome Y, and the last one is the index D
    rr<-NULL
    n=dim(mydata)[1]

    ## number of covariates
    p = dim(mydata)[2]-1

    ## Standardize covariates X and center outcome Y

    x = mydata[,1:p]
    y = mydata[,(p+1)]

    meany = mean(y)
    y= y - meany
    y<-matrix(y,nrow=n)

    ##############GEE##############
    w=abs(b)^(-gamma)

    iter=0
    dif=1

    V0_ni=solve(V0)


    sum_sigma=0
    for(i in 1:yuangeshu){
      sigma=t(y[((i-1)*D+1):(i*D),]-x[((i-1)*D+1):(i*D),]%*%b)%*%V0_ni%*%(y[((i-1)*D+1):(i*D),]-x[((i-1)*D+1):(i*D),]%*%b)/(n-p)
      sum_sigma=sum_sigma+sigma
    }
    sum_sigma<-c(sum_sigma)

    ####???###
    cc<-w/sqrt(b^2)
    cc<-diag(c(cc))

    while(iter<=maxiter & dif>=eps){
      iter=iter+1

      ## beta
      b.old=b

      sum_xx=0
      sum_yy=0


      V0_ni=solve(V0)

      for(i in 1:yuangeshu){
        xx<-t(x[((i-1)*D+1):(i*D),])%*%V0_ni%*%x[((i-1)*D+1):(i*D),]/sum_sigma###HAO XIANG YOU WEN TI
        sum_xx<-sum_xx+xx
      }
      sum_xx<-sum_xx+2*lamda*cc

      for(i in 1:yuangeshu){
        yy<-t(x[((i-1)*D+1):(i*D),])%*%V0_ni%*%(y[((i-1)*D+1):(i*D),]-x[((i-1)*D+1):(i*D),]%*%b.old)/sum_sigma###HAO XIANG YOU WEN TI
        sum_yy<-sum_yy+yy
      }
      sum_yy<-sum_yy-2*lamda*cc%*%b.old


      blog <- try(solve(sum_xx), silent = TRUE)
      if ('try-error' %in% class(blog)){
        #iter=maxiter+1
      }else{
        xx_ni<-solve(sum_xx)

        b=b.old+xx_ni%*%sum_yy


        ## Rho moment estimation ##if chale qitacishu ze gengai rx de shumu liru r6 r7 dengdeng
        r=y-x%*%b

        r_time=n/yuangeshu

        r_list=list()
        r_string=""
        for (r_i in 1:r_time){
          exp11_chr=paste("r",as.character(r_i),"=rep(0,yuangeshu)",sep="")
          eval(parse(text=exp11_chr))
          r_name=paste("r",as.character(r_i),sep="")
          r_list=append(r_list,r_name)
          r_string=paste(r_string,r_name,sep=",")
        }

        for (r_i in 1:(r_time-1)){
          exp11_chr=paste("for(i in 1:yuangeshu){r",as.character(r_i),"[i]=r[(i-1)*D+",as.character(r_i),"]}",sep="")
          eval(parse(text=exp11_chr))}
        exp11_chr=paste("for(i in 1:yuangeshu){r",as.character(r_time),"[i]=r[i*D]}",sep="")
        eval(parse(text=exp11_chr))

        cbind_exp_chr=paste("rr<-cbind(",substring(r_string,2),")",sep="")
        eval(parse(text=cbind_exp_chr))

        rho=0
        for(j in 1:5){
          for(k in 1:5){
            if(j!=k){
              rho=rho+cor(rr[,j],rr[,k])
            }
          }
        }

        rho_sum=rho/(D*(D-1))
        ## working matrix

        V0=rho_sum*matrix(1,D,D)+(1-rho_sum)*diag(D)


        ## sigma
        sum_sigma=0


        V0_ni=solve(V0)


        for(i in 1:yuangeshu){
          sigma=t(y[((i-1)*D+1):(i*D),]-x[((i-1)*D+1):(i*D),]%*%b)%*%V0_ni%*%(y[((i-1)*D+1):(i*D),]-x[((i-1)*D+1):(i*D),]%*%b)/(n-p)
          sum_sigma=sum_sigma+sigma
        }
        sum_sigma<-c(sum_sigma)

        ## CC
        cc<-1/sqrt(b^2)
        cc<-c(cc)
        cc[cc>1/(1e-10)] =1/(1e-10)
        cc=matrix(cc,nrow=p)
        cc<-diag(c(w*cc))


        dif = max(abs(b-b.old))

      }
    }
    b[abs(b)<=1e-4,] = 0

    return(list(beta=b,V=V0))
  }

  g1 <-  mi_data <-  mi_entire <- NULL
  missdata=as.matrix.data.frame(missdata)
  ##PWLS
  pp=ncol(missdata)-1
  n=nrow(missdata)

  mice_data_list=list()
  mice_data_string=""

  if (requireNamespace("mice",quietly = TRUE)){
    for(mice_i in 1:mice_time)
    {
      exp1_chr=paste("mi_entire_",as.character(mice_i),"<-mice::complete(mice::mice(missdata,seed=100,printFlag=FALSE,m=",as.character(mice_time),"),action=",as.character(mice_i),")",sep="")
      eval(parse(text=exp1_chr))
      exp2_chr=paste("mi_entire_",as.character(mice_i),"$id<-1:n",sep="")
      eval(parse(text=exp2_chr))

      mice_name=paste("mi_entire_",as.character(mice_i),sep="")
      mice_data_list=append(mice_data_list,mice_name)
      mice_data_string=paste(mice_data_string,mice_name,sep=",")
    }
  }else{
    return(NULL)
  }



  rbind_exp_chr=paste("mi_entire<-rbind(",substring(mice_data_string,2),")",sep="")
  eval(parse(text=rbind_exp_chr))

  mi_entire_ord<-mi_entire[order(mi_entire$id),]
  mi_entire_ord$id<-gl(n,mice_time,mice_time*n)

  mi_entire_ord<-as.matrix(mi_entire_ord[,1:(pp+1)])

  bx=mi_entire_ord[,1:pp]
  by=mi_entire_ord[,(pp+1)]

  meanbx = apply(bx, 2, mean)
  bx = scale(bx, meanbx, FALSE)

  meanby = mean(by)
  by= by - meanby

  by<-matrix(by,nrow=mice_time*n)


  b = qr.solve(t(bx)%*%bx)%*%t(bx)%*%by
  est<-GEEsym.estimate(mydata=mi_entire_ord,yuangeshu=n,D=mice_time,b=b,maxiter=200,eps=1e-6)
  bb<-est$beta

  #######select tuning paramete

  if(penalty=="alasso"){

    QGCVLAM = matrix(0,length(Gamma),length(lamda.vec))
    SSE=matrix(0,length(Gamma),length(lamda.vec))
    jieguo=NULL
    xiechazhen=NULL

    for(jj in 1:length(Gamma)){
      for (j in 1:length(lamda.vec)) {
        estt<-GEEsym_gl.alasso.estimate(mydata=mi_entire_ord,yuangeshu=n,D=mice_time,b=bb,gamma=Gamma[jj],V0=est$V,lamda=lamda.vec[j],maxiter=200,eps=1e-6)
        bbp<-estt$beta
        jieguo=cbind(jieguo,bbp)
        xiechazhen=c(xiechazhen,estt$V)

        mx<-mi_entire_ord[,-(pp+1)]
        re<-mi_entire_ord[,(pp+1)]-mx%*%bbp
        ni<-solve(estt$V)


        sse=0
        ccid= which(!is.na(missdata[,(pp+1)]))

        for(k in ccid){
          sss<-re[k*mice_time]^2
          sse<-sse+sss
        }
        SSE[jj,j]=sse

        ####QGCV
        sum_xrx=0

        for(i in which(is.na(missdata[,(pp+1)])) ){
          xrx<-t(re[((i-1)*mice_time+1):(i*mice_time),])%*%ni%*%re[((i-1)*mice_time+1):(i*mice_time),]
          sum_xrx<-sum_xrx+xrx
        }

        wdev<-sum_xrx
        p=length(which(bbp!=0))*sum(abs(bbp))/sum(abs(b))

        N=length(which(is.na(missdata[,(pp+1)])))*mice_time*mice_time/sum(estt$V)

        QGCV=wdev/(length(which(is.na(missdata[,(pp+1)])))*(1-(p/N))^2)
        #


        QGCVLAM[jj,j]=QGCV


      }
    }



    geshu=0
    for(ii in 1:ncol(jieguo)){
      geshu[ii]=length(which(jieguo[,ii]!=0))
    }
    #####a new tuning parameter seletion method
    w=length(ccid)/n

    wic=(1-w)*QGCVLAM+(-2*log(SSE)+geshu*log(length(ccid)))*w
    zhi=data.frame(WIC=c(t(wic)),gamma=rep(Gamma,each=length(lamda.vec)),lambda=rep(lamda.vec,length(Gamma)))#???3=length(gamma)

    zhi$gamma_label=factor(rep(Gamma,each=length(lamda.vec)))
  }else if(penalty=="lasso"){
    QGCVLAM = matrix(0,1,length(lamda.vec))
    SSE=matrix(0,1,length(lamda.vec))
    jieguo=NULL
    xiechazhen=NULL
    for (j in 1:length(lamda.vec)) {
      estt<-GEEsym_gl.lasso.estimate(mydata=mi_entire_ord,yuangeshu=n,D=mice_time,b=bb,V0=est$V,lamda=lamda.vec[j],maxiter=200,eps=1e-6)

      bbp<-estt$beta
      jieguo=cbind(jieguo,bbp)
      xiechazhen=c(xiechazhen,estt$V)

      mx<-mi_entire_ord[,-(pp+1)]
      re<-mi_entire_ord[,(pp+1)]-mx%*%bbp
      ni<-solve(estt$V)

      ######bic
      sse=0
      ccid= which(!is.na(missdata[,(pp+1)]))

      for(k in ccid){
        sss<-re[k*mice_time]^2
        sse<-sse+sss
      }
      SSE[1,j]=sse

      ####QGCV
      sum_xrx=0

      for(i in which(is.na(missdata[,(pp+1)])) ){
        xrx<-t(re[((i-1)*mice_time+1):(i*mice_time),])%*%ni%*%re[((i-1)*mice_time+1):(i*mice_time),]
        sum_xrx<-sum_xrx+xrx
      }

      wdev<-sum_xrx
      p=length(which(bbp!=0))*sum(abs(bbp))/sum(abs(b))

      N=length(which(is.na(missdata[,(pp+1)])))*mice_time*mice_time/sum(estt$V)

      QGCV=wdev/(length(which(is.na(missdata[,(pp+1)])))*(1-(p/N))^2)
      #


      QGCVLAM[1,j]=QGCV


    }



    geshu=0
    for(ii in 1:ncol(jieguo)){
      geshu[ii]=length(which(jieguo[,ii]!=0))
    }
    #####a new tuning parameter seletion method
    w=length(ccid)/n
    wic=(1-w)*QGCVLAM+(-2*log(SSE)+geshu*log(length(ccid)))*w
    zhi=data.frame(WIC=c(t(wic)),lambda=lamda.vec)#???3=length(gamma)
  }else{
    message("Penalty should be alasso or lasso")
    return(NULL)
  }


  ind=which.min(zhi$WIC)
  bbp=matrix(jieguo[,ind],nrow=pp)
  chazhenguji=xiechazhen[ind]


  index_juti<-0
  k=1
  for(i in 1:dim(bbp)[1]){
    if(abs(bbp[i,1])>0){
      index_juti[k]=i
      k=k+1
    }
  }
  estcoef=as.vector(bbp)

  Vsmi_est <- function(x, y) {
    obj <- list(estcoef = x, index_sig = y)
    class(obj) <- "vsmi_est"
    return(obj)
  }


  obj <- Vsmi_est(estcoef, index_juti)
  print.vsmi <- function(obj) {
    message("Estimate coefficients:\n", obj$estcoef,"\n")
    message("Selected variable index:\n", obj$index_sig, "\n")
  }
  print.vsmi(obj)

  return(obj)
}


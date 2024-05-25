
#' Penalized estimating equations for generalized linear models with multiple imputation
#'
#'This is a function to impute missing data, estimate coefficients of generalized linear models and select variables for multiple imputed data sets,
#'considering the correlation of multiple imputed observations.
#'
#'
#'@param missdata A Matrix,missing data with variables X in the first p columns and response Y at the last column.
#'@param mice_time an integer, number of imputation.
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
#'
#'data_with_missing <- generate_pee_missing_data(outcome="binary")
#'est.alasso <-PEE(data_with_missing,penalty="alasso")
#'est.lasso <-PEE(data_with_missing,penalty="lasso")
#'
#'count_data_with_missing <- generate_pee_missing_data(outcome="count")
#'count_est.alasso <-PEE(data_with_missing,penalty="alasso")
#'count_est.lasso <-PEE(data_with_missing,penalty="lasso")
#'}
PEE <- function(missdata,mice_time=5,penalty,lamda.vec =seq(1,4,length.out=12),Gamma=c(0.5,1,1.5)){

  ##inner functions
  logit_QL_gl.alasso.estimate <- function(mydata,yuangeshu,D,b,R0,penalty,lamda,gamma,maxiter=200,eps=1e-3){
    rr <- NULL
    ## for each mydata, the first p columns are covariates X, and the last two is the outcome Y, and the last one is the index D
    n=dim(mydata)[1]

    ## number of covariates
    p = dim(mydata)[2]-1

    ## Standardize covariates X and center outcome Y
    ## note that for continuous x, it should be processed but not for categorical x
    ## for continuous Y, same and for categorical Y not
    x = mydata[,1:p]
    y = mydata[,(p+1)]

    y<-matrix(y,nrow=n)

    ## beta
    b.old=b

    sum_xx=0
    sum_yy=0

    A_chu=c(exp(x%*%b.old)/(1+exp(x%*%b.old))^2)
    A0=list()
    for(i in 1:yuangeshu){
      A0[[i]]=diag(A_chu[((i-1)*D+1):(i*D)])
    }

    if (penalty=="alasso"){
      w=abs(b)^(-gamma)
    }else{
      w=rep(1,p)
    }

    iter=0
    dif=1
    V0=list()
    for(i in 1:yuangeshu){
      V0[[i]]=sqrt(A0[[i]])%*%R0%*%sqrt(A0[[i]])
    }

    V0_ni=list()
    for(i in 1:yuangeshu){
      if(det(V0[[i]])>1e-5){
        V0_ni[[i]]=qr.solve(V0[[i]])
      } else {
        T<-svd(V0[[i]])
        index=which(T$d>1e-5)
        T$d[index]<-1/T$d[index]
        T$d[-index]<-0
        V0_ni[[i]]<-T$v %*% diag(T$d) %*% t(T$u)
      }
    }

    cc<-w/sqrt(b^2)
    cc<-diag(c(cc))

    while(iter<=maxiter & dif>=eps){
      iter=iter+1

      ## beta
      b.old=b

      sum_xx=0
      sum_yy=0

      V0_ni=list()
      blog <- try(det(V0[[1]]), silent = TRUE)
      if ('try-error' %in% class(blog)){
        iter=maxiter+1
      }else{
        for(i in 1:yuangeshu){
          if(det(V0[[i]])>1e-5){
            V0_ni[[i]]=qr.solve(V0[[i]])
          } else {
            T<-svd(V0[[i]])
            index=which(T$d>1e-5)
            T$d[index]<-1/T$d[index]
            T$d[-index]<-0
            V0_ni[[i]]<-T$v %*% diag(T$d) %*% t(T$u)
          }
        }
        mu=exp(x%*%b.old)/(1+exp(x%*%b.old))
        mujian=exp(x%*%b.old)/(1+exp(x%*%b.old))^2
        mudao=matrix(0,n,p)
        for(i in 1:n){
          mudao[i,]=mujian[i,]*x[i,]
        }

        for(i in 1:yuangeshu){
          xx<-t(mudao[((i-1)*D+1):(i*D),])%*%V0_ni[[i]]%*%mudao[((i-1)*D+1):(i*D),]
          sum_xx<-sum_xx+xx
        }
        sum_xx<-sum_xx+2*lamda*cc

        for(i in 1:yuangeshu){
          yy<-t(mudao[((i-1)*D+1):(i*D),])%*%V0_ni[[i]]%*%(y[((i-1)*D+1):(i*D),]- mu[((i-1)*D+1):(i*D),])
          sum_yy<-sum_yy+yy
        }

        sum_yy<-sum_yy-2*lamda*cc%*%b.old


        blog <- try(solve(sum_xx), silent = TRUE)
        if ('try-error' %in% class(blog)){
          iter=maxiter+1
        }else{
          xx_ni<-solve(sum_xx)
          b=b.old+xx_ni%*%sum_yy

          ## Rho moment estimation ##if chale qitacishu ze gengai rx de shumu liru r6 r7 dengdeng
          gailv=exp(x%*%b)/(1+exp(x%*%b))

          id1=which(y==0)
          id2=which(y==1)
          r=matrix(rep(0,n),nrow=n)
          r[id1,]=-sqrt(2*(-log(1-gailv[id1,])))
          r[id2,]=sqrt(2*-log(gailv[id2,]))

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

          A_chu=c(exp(x%*%b)/(1+exp(x%*%b))^2)
          A0=list()
          for(i in 1:yuangeshu){
            A0[[i]]=diag(A_chu[((i-1)*D+1):(i*D)])
          }

          V0=list()
          for(i in 1:yuangeshu){
            V0[[i]]=sqrt(A0[[i]])%*%((rho/(D*(D-1)))*matrix(1,D,D)+(1-(rho/(D*(D-1))))*diag(D)+0.01*diag(D))%*%sqrt(A0[[i]])
          }

          ## CC
          cc<-1/sqrt(b^2)
          cc<-c(cc)
          cc[cc>1/(1e-10)] =1/(1e-10) ###kaolv shidangde fangkuai zuixiao xianzhi
          cc=matrix(cc,nrow=p)
          cc<-diag(c(w*cc))

          dif = max(abs(b-b.old))
        }
      }
    }
    b[abs(b)<=eps,] = 0

    return(list(beta=b,V=R0,A=A0))
  }

  Poisson_QL_gl.alasso.estimate <- function(mydata,yuangeshu,D,b,R0,penalty,lamda,gamma,maxiter=200,eps=1e-3){
    rr <- NULL
    ## for each mydata, the first p columns are covariates X, and the last two is the outcome Y, and the last one is the index D
    #library(Matrix)
    #library(MASS)
    ## number of observations
    ## number of observations
    n=dim(mydata)[1]

    ## number of covariates
    p = dim(mydata)[2]-1

    ## Standardize covariates X and center outcome Y
    ## note that for continuous x, it should be processed but not for categorical x
    ## for continuous Y, same and for categorical Y not
    x = mydata[,1:p]
    y = mydata[,(p+1)]


    y<-matrix(y,nrow=n)
    x=cbind(rep(1,n),x)

    ## beta
    b.old=b

    sum_xx=0
    sum_yy=0

    A_chu=c(exp(x%*%b.old))
    A0=list()
    for(i in 1:yuangeshu){
      A0[[i]]=diag(A_chu[((i-1)*D+1):(i*D)])
    }

    if(penalty=="alasso"){
      w=abs(b)^(-gamma)
    }else{w=1}

    iter=0
    dif=1
    V0=list()
    for(i in 1:yuangeshu){
      V0[[i]]=sqrt(A0[[i]])%*%R0%*%sqrt(A0[[i]])
    }

    V0_ni=list()
    for(i in 1:yuangeshu){

      V0_ni[[i]]=solve(V0[[i]])

    }

    cc<-w/sqrt(b^2)
    cc<-diag(c(cc))

    while(iter<=maxiter & dif>=eps){
      iter=iter+1

      ## beta
      b.old=b

      sum_xx=0
      sum_yy=0

      V0_ni=list()
      blog <- try(det(V0[[1]]), silent = TRUE)
      if ('try-error' %in% class(blog)){
        iter=maxiter+1
      }else{
        for(i in 1:yuangeshu){
          V0_ni[[i]]=solve(V0[[i]])
        }
        mu=exp(x%*%b.old)
        mujian=exp(x%*%b.old)
        mudao=matrix(0,n,p+1)
        for(i in 1:n){
          mudao[i,]=mujian[i,]*x[i,]
        }

        for(i in 1:yuangeshu){
          xx<-t(mudao[((i-1)*D+1):(i*D),])%*%V0_ni[[i]]%*%mudao[((i-1)*D+1):(i*D),]
          sum_xx<-sum_xx+xx
        }
        sum_xx<-sum_xx+2*lamda*cc

        for(i in 1:yuangeshu){
          yy<-t(mudao[((i-1)*D+1):(i*D),])%*%V0_ni[[i]]%*%(y[((i-1)*D+1):(i*D),]- mu[((i-1)*D+1):(i*D),])
          sum_yy<-sum_yy+yy
        }

        sum_yy<-sum_yy-2*lamda*cc%*%b.old

        blog <- try(solve(sum_xx), silent = TRUE)
        if ('try-error' %in% class(blog)){
          iter=maxiter+1
        }else{
          xx_ni<-solve(sum_xx)
          b=b.old+xx_ni%*%sum_yy

          ## Rho moment estimation ##if chale qitacishu ze gengai rx de shumu liru r6 r7 dengdeng
          gailv=exp(x%*%b)


          r=sign(y-gailv)*sqrt(abs(2*(y*log((y+0.001)/(gailv+0.001))-(y-gailv))))

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



          A_chu=c(exp(x%*%b))
          A0=list()
          for(i in 1:yuangeshu){
            A0[[i]]=diag(A_chu[((i-1)*D+1):(i*D)])
          }


          V0=list()
          for(i in 1:yuangeshu){
            V0[[i]]=sqrt(A0[[i]])%*%R0%*%sqrt(A0[[i]])
          }

          ## CC
          cc<-1/sqrt(b^2)
          cc<-c(cc)
          cc[cc>1/(1e-10)] =1/(1e-10) ###kaolv shidangde fangkuai zuixiao xianzhi
          cc=matrix(cc,nrow=p+1)
          cc<-diag(c(w*cc))


          dif = max(abs(b-b.old))
        }
      }
    }
    b[abs(b)<=eps,] = 0


    return(list(beta=b,V=R0))
  }
  g1 <-  mi_data <-  mi_entire <- NULL
  missdata=as.matrix.data.frame(missdata)
  ##PEE
  n=nrow(missdata)
  pp=ncol(missdata)-1
  if(length(levels(as.factor(missdata[,(pp+1)])))<=2){
    outcome="binary"
  }else{
    outcome="count"
  }

  miss_id <- (1:(pp+1))[apply(missdata, 2, function(x) any(is.na(x)))]
  miss_p <-length(miss_id)

  m1=rep("pmm",miss_p)
  blk=c()
  pmx_row=rep(1,pp+1)
  for (pmx_i in miss_id){
    pmx_row[pmx_i]=0
    miss_name_chr=paste("V",as.character(pmx_i),sep="")
    blk=append(blk,miss_name_chr)
  }

  pmx_rbind=pmx_row
  for (pmx_i in head(miss_id,-1)){
    pmx_row[pmx_i]=1
    pmx_rbind=rbind(pmx_rbind,pmx_row)
  }
  pmx=as.matrix(pmx_rbind)
  rownames(pmx)<-NULL


  mice_data_list=list()
  mice_data_string=""

  if (requireNamespace("mice",quietly = TRUE)){
    for(mice_i in 1:mice_time)
    {
      exp1_chr=paste("mi_entire_",as.character(mice_i),"<-mice::complete(mice::mice(missdata,method=m1,printFlag=FALSE,predictorMatrix = pmx,blocks = blk,visitSequence = 'roman',seed=",as.character(mice_time),"00),action=1)",sep="")
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



  rbind_exp_chr=paste("mi_data<-rbind(",substring(mice_data_string,2),")",sep="")
  eval(parse(text=rbind_exp_chr))
  mi_entire_ord<-mi_data[order(mi_data$id),]
  mi_entire_ord$id<-gl(n,mice_time,mice_time*n)

  #estimate
  mi_entire_ord<-as.matrix(mi_entire_ord[,1:(pp+1)])

  bx=mi_entire_ord[,1:pp]

  meanbx = apply(bx, 2, mean)
  bx = scale(bx, meanbx, FALSE)

  dairu=as.data.frame(mi_entire_ord)
  dairu[,1:pp]=bx
  if(outcome=="binary"){
    exp_g1=paste("g1<-glm(V",as.character(pp+1),"~.,data=dairu,family = binomial(link = 'logit'))",sep="")
  }else{
    exp_g1=paste("g1<-glm(V",as.character(pp+1),"~.,data=dairu,family = poisson())",sep="")

  }
  eval(parse(text=exp_g1))
  co=matrix(g1$coefficients[-1],ncol=1)

  b = co


  #PQIF estimate

  if(outcome=="binary"){
    estimate_func <- logit_QL_gl.alasso.estimate
  }else{
    estimate_func <- Poisson_QL_gl.alasso.estimate
  }


  if(penalty=="lasso"){
    Gamma=c(0)
  }

  QGCVLAM = matrix(0,length(Gamma),length(lamda.vec))
  SSE=matrix(0,length(Gamma),length(lamda.vec))
  jieguo=NULL
  xiechazhen=NULL

  blog<-try(estimate_func(mydata=mi_entire_ord,yuangeshu=n,penalty=penalty,D=mice_time,b=co,gamma=Gamma[1],R0=diag(mice_time),lamda=lamda.vec[1],maxiter=200,eps=1e-4),silent = TRUE)
  if ('try-error' %in% class(blog)){
    message("error")
  }else{
    for(jj in 1:length(Gamma)){
      for (j in 1:length(lamda.vec)) {
        estt<-estimate_func(mydata=mi_entire_ord,yuangeshu=n,penalty=penalty,D=mice_time,b=co,gamma=Gamma[jj],R0=diag(mice_time),lamda=lamda.vec[j],maxiter=200,eps=1e-4)
        bbp<-estt$beta

        jieguo=cbind(jieguo,bbp)
        xiechazhen=c(xiechazhen,estt$V)
        mx<-mi_entire_ord[,-(pp+1)]
        gailv=exp(mx%*%bbp)/(1+exp(mx%*%bbp))

        my=mi_entire_ord[,pp+1]
        Q=matrix(0,n*mice_time,1)
        for(i in 1:n*mice_time){
          if(my[i]==1){
            Q[i]=log(gailv[i])
          }else {
            Q[i]=log(1-gailv[i])
          }
        }
        re=sign(my-gailv)*sqrt(-2*Q)

        ni<-solve(estt$V)

        ######bic
        sse=0
        ccid= which(!is.na(missdata[,pp+1]))

        for(k in ccid){
          sss<-mi_entire_ord[k*mice_time,pp+1]*(mx[k*mice_time,]%*%bbp)-log(1+exp(mx[k*mice_time,]%*%bbp))
          sse<-sse+sss
        }
        SSE[jj,j]=sse
      }
    }

    geshu=0
    for(ii in 1:ncol(jieguo)){
      geshu[ii]=length(which(jieguo[,ii]!=0))
    }

    geshu=matrix(geshu,nrow=length(lamda.vec))
    geshu=t(geshu)
    #####a new tuning parameter seletion method
    w=length(ccid)/200
    wic=(1-w)*QGCVLAM+(-2*log(SSE-min(SSE)+1)+geshu*log(length(ccid)))*w
    zhi=data.frame(WIC=c(t(wic)),gamma=rep(Gamma,each=length(lamda.vec)),lambda=rep(lamda.vec,length(Gamma)))

    zhi$gamma_label=factor(rep(Gamma,each=length(lamda.vec)))

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
  }



  #estimation using QIF

  estdat <- as.data.frame(mi_entire_ord)
  estdatY <- estdat[,pp+1]
  estdatX <- estdat[,which(bbp!=0)]

  if (requireNamespace("qif",quietly = TRUE)){
    if(outcome=="binary"){
      fit <- qif::qif(estdatY ~ as.matrix(estdatX)-1, id=rep(1:n,each=mice_time),family=binomial)
    }
    else{
      fit <- qif::qif(estdatY ~ as.matrix(estdatX)-1, id=rep(1:n,each=mice_time),family=poisson)
    }
  }else{
    return(NULL)
  }


  estxishu <- rep(0,pp)
  estxishu[which(bbp!=0)] <- fit$coefficients

  ###output of PEE
  estcoef=estxishu
  index_sig=index_juti


  Vsmi_est <- function(x, y) {
    obj <- list(estcoef = x, index_sig = y)
    class(obj) <- "vsmi_est"
    return(obj)
  }


  obj <- Vsmi_est(estxishu, index_juti)
  print.vsmi <- function(obj) {
    message("Estimate coefficients:\n", obj$estcoef,"\n")
    message("Selected variable index:\n", obj$index_sig, "\n")
  }
  print.vsmi(obj)

  return(obj)
}

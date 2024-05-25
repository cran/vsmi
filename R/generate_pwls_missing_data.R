#' Generate example data for PWLS
#'
#'This is a functoin to generate example missing data for PWLS
#'
#' @param p The dimension of the independent variable X,default 20.
#' @param n The Number of rows of generated data,default 200.
#' @param pt1 Missing rate of independent variable X,default 0.5.
#' @param pt2 Missing rate of response Y, default 0.5.
#' @param tbeta True value of the coefficient,default c(1,-1,1,-1,1,-1,-1,1).
#' @param miss_sig A 0-1 vector of length p, where 1 means that variable at the index is with missing,while 0 means that it without missing,defualt c(0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0)
#'
#' @return A Matrix,missing data with variables X in the first p columns and response Y at the last column.
#' @export

generate_pwls_missing_data <-
  function(p=20,n=200,pt1=0.5,pt2=0.5,
           tbeta=c(1,-1,1,-1,1,-1,-1,1),
           miss_sig=c(0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0)){
    rr <- f1 <-   mi_data <-  mi_entire <-g1 <- NULL
    if(length(miss_sig)!=p){
      message("length(mis_sig) should equal to p")
      return(NULL)
    }

    miss_id <- (1:p)[miss_sig==1]
    unmiss_id <- (1:p)[miss_sig!=1]
    miss_p=length(miss_id)

    if(miss_p>1/3*p){
      message("number of miss variable should less than (1/3)*p,Please change miss_sig")
      return(NULL)
    }
    f1<-NULL
    mean=matrix(c(rep(0,p)),nrow=p)
    sigma=diag(2.7,p)+matrix(0.3,p,p)

    if (requireNamespace("MASS",quietly = TRUE)){
      x<-MASS::mvrnorm(n=n,mean,sigma)
    }else{
      return(NULL)
    }


    emean=matrix(c(rep(0,n)),nrow=n)
    esigma=3*diag(n)

    if (requireNamespace("MASS",quietly = TRUE)){
      ee1<-MASS::mvrnorm(n=1,emean,esigma)
    }else{
      return(NULL)
    }



    p_beta <- length(tbeta)
    if(p_beta<p){
      tbeta <-c(tbeta,rep(0,p-p_beta))
    }
    y1 <- x%*%tbeta+ee1
    y1<-matrix(y1,nrow=n)



    a<-matrix(c(rep(0,miss_p)),nrow=1)

    for (f_i in 1:miss_p){
      exp_f <- paste("f",as.character(f_i),"<-function(aa){sum(exp(aa+x[,unmiss_id[f_i]])/(1+exp(aa+x[,unmiss_id[f_i]])))-pt1*n}",sep="")
      eval(parse(text=exp_f))
      exp_root <- paste("a[1,f_i]<-uniroot(f",as.character(f_i),",c(-10,10))$root",sep="")
      eval(parse(text=exp_root))

    }



    p_mat<-matrix(rep(0,n*miss_p),nrow=n)

    for(i in 1:n){
      for (f_i in 1:miss_p){
        p_mat[i,f_i]=exp(a[1,f_i]+x[i,unmiss_id[f_i]])/(1+exp(a[1,f_i]+x[i,unmiss_id[f_i]]))
      }

    }

    ################
    k<-matrix(rep(0,n*miss_p),nrow=n)
    for(i in 1:n){
      for (f_i in 1:miss_p){
        k[i,f_i]=rbinom(1,1,p_mat[i,f_i])
      }

    }

    #############
    newx<-x
    for(i in 1:n){
      for (k_i in 1:miss_p){
        if(k[i,k_i]==1){
          newx[i,miss_id[k_i]]=NA
        }
      }

    }



    a1<-matrix(c(rep(0,1)),nrow=1)


    add_list1=""
    add_list2=""
    for (ms_i in 1:miss_p){
      add_list1=paste(add_list1,'+x[,',as.character(unmiss_id[ms_i+miss_p]),"]",sep='')
      add_list2=paste(add_list2,'+x[i,',as.character(unmiss_id[ms_i+miss_p]),"]",sep='')
    }


    f1_exp=paste("f1<-function(aa){sum(exp(aa",add_list1,")/(1+exp(aa",add_list1,")))-pt2*n}",sep="")
    eval(parse(text=f1_exp))

    a1[1,1]<-uniroot(f1,c(-10,10))$root

    p1<-matrix(rep(0,n),nrow=n)

    for(i in 1:n){
      p1_exp=paste("p1[i,1]=exp(a1[1,1]",add_list2,")/(1+exp(a1[1,1]",add_list2,"))",sep="")
      eval(parse(text=p1_exp))
    }


    k1<-matrix(rep(0,n),nrow=n)
    for(i in 1:n){

      k1[i,1]=rbinom(1,1,p1[i,1])

    }
    newy<-y1
    for(i in 1:n){

      if(k1[i,1]==1){
        newy[i,1]=NA
      }

    }

    entire<-cbind(newx,newy)
    #class(entire) <- "MissingData"
    message("Generated data with missing for PWLS.\n")
    return(entire)
  }

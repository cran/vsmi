#' Generate example data for PEE
#'
#'This is a functoin to generate example missing data for PEE
#'
#'@param outcome The type of response variable Y, choose "binary" for binary response or "count" for poisson response,defualt "binary"
#'@param p The dimension of the independent variable X,default 20.
#'@param n The Number of rows of generated data,default 200.
#'@param pt1 Missing rate of independent variable X,default 0.5.
#'@param tbeta True value of the coefficient,default c(3/4,(-3)/4,3/4,(-3)/4,3/4,(-3)/4,(-3)/4,3/4).
#'@param miss_sig A 0-1 vector of length p, where 1 means that variable at the index is with missing,while 0 means that it without missing,defualt c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
#'
#'@return A Matrix,missing data with variables X in the first p columns and response Y at the last column.
#' @export

generate_pee_missing_data <-
  function(outcome="binary",p=20,n=200,pt1=0.5,
           tbeta=c(3/4,(-3)/4,3/4,(-3)/4,3/4,(-3)/4,(-3)/4,3/4),
           miss_sig = c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)){
    rr <- f1 <-   mi_data <-  mi_entire <- g1 <- NULL
    if(length(miss_sig)!=p){
      message("length(mis_sig) should equal to p")
      return(0)
    }

    miss_id <- (1:p)[miss_sig==1]
    unmiss_id <- (1:p)[miss_sig!=1]
    p_miss=length(miss_id)

    if(p_miss>1/2*p){
      message("number of miss variable should less than (1/2)*p,Please change miss_sig")
      return(0)
    }
    f1<-NULL
    mean <- matrix(c(rep(0,p)),nrow=p)
    sigma <- diag(0.9,p)+matrix(0.1,p,p)
    if (requireNamespace("MASS",quietly = TRUE)){
      x<-MASS::mvrnorm(n,mean,sigma)
    }else{
      return(NULL)
    }




    p_beta <- length(tbeta)
    if(p_beta<p){
      tbeta <-c(tbeta,rep(0,p-p_beta))
    }

    if (outcome=="binary"){
      p1 <- 1/(1+exp(-x%*%tbeta))
      y1 <- matrix(0,n,1) ##outcome
      for(ii in 1:nrow(p1)){
        y1[ii,1]=rbinom(1,1,p1[ii,1])
      }
    }else{
      r1<-x%*%tbeta
      r1<-matrix(r1,nrow=n)
      lam1<-exp(r1)
      y1=matrix(0,n,1)
      for(ii in 1:nrow(lam1)){
        y1[ii,1]=rpois(1,lam1[ii,1])
      }
    }



    a<-matrix(c(rep(0,p_miss)),nrow=1)
    ########
    for (f_i in 1:p_miss){
      exp_f=paste("f",as.character(f_i),
                  "<-function(aa){sum(exp(aa+x[,unmiss_id[f_i]])/(1+exp(aa+x[,unmiss_id[f_i]])))-pt1*n}",sep="")
      eval(parse(text=exp_f))
      exp_root=paste("a[1,f_i]<-uniroot(f",as.character(f_i),",c(-10,10))$root",sep="")
      eval(parse(text=exp_root))
    }


    ############
    ps<-matrix(rep(0,p_miss*n),nrow=n)

    for(i in 1:n){
      for(ps_i in 1:p_miss){
        ps[i,ps_i]=exp(a[1,ps_i]+x[i,unmiss_id[ps_i]])/(1+exp(a[1,ps_i]+x[i,unmiss_id[ps_i]]))
      }

    }

    ################
    k<-matrix(rep(0,p_miss*n),nrow=n)
    for(i in 1:n){
      for (k_i in 1:p_miss){
        k[i,k_i]=rbinom(1,1,ps[i,k_i])
      }

    }

    #############
    newx<-x
    for(i in 1:n){
      for (na_i in 1:p_miss){
        if(k[i,na_i]==1){
          newx[i,miss_id[na_i]]=NA
        }
      }

    }

    data_with_missing <-as.matrix(cbind(newx,y1)) ##
    #class(data_with_missing) <-"Missingdata"
    message("Generated data with missing for PEE.\n")
    return(data_with_missing)
  }

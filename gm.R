#' Gaussian Mirrors for Controlled Variable Selection
#'
#' 
#' @param y response variable
#' @param x predictor variables 
#' @param q target FDR level
#' @param family 
#'
#' @return result: the index of selected variable
#' @return 
#'
#' @author Xin Xing, \email{xin_xing@stat.harvard.edu}
#' @references \url{}
#' @keywords variable selection
#'
#' @examples
#' gm(y, x)
#'
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom    

gm  = function(y, x, family="gaussian", ncores = 1)

    n = length(y)
    p = dim(x)[2]
    x= scale(x)
    lambda_max = 1*max(abs(t(x) %*% y))/(n^(3/2))
    lambda_min = lambda_max / 1e2
    nlambda=100
    k = (0:(nlambda-1)) / nlambda
    lambda = lambda_max * (lambda_min/lambda_max)^k
    fit0.lasso = cv.glmnet(x, y, family=family, lambda=lambda, nfold=10,alpha=1)
    lambda_ast = fit0.lasso$lambda.min

    coef_ast = coef(fit0.lasso, s='lambda.min')[-1]
    if(sum(coef_ast!=0)>0.9*n){
        ord_idx = order(abs(coef_ast[which(coef_ast!=0)]))
        coef_ast[which(coef_ast!=0)[1:(sum(coef_ast!=0)-0.9*n)]]=0
    }
    sigmaz_vec = numeric(p)
    z = matrix(0, nrow=n, ncol=p)


    for(j in 1:p){
        jidx =  intersect(which(coef_ast!=0), c(1:p)[-j])
        xj =  x[,j]
        xnj = x[,jidx]
        xnj.svd = svd(xnj)
        xnj.idx=1:(min(0.95*length(jidx),n))
        pxnj = xnj.svd$u[, xnj.idx] 
        pp= (1 - t(xj) %*% pxnj  %*%  solve(t(pxnj) %*% pxnj) %*% t(pxnj) %*% xj/n)

        if(pp<0.05){
            xnj.idx=1:10
            pxnj = xnj.svd$u[, xnj.idx] 
            sigmaz_vec[j] = (1 - t(xj) %*% pxnj  %*%  solve(t(pxnj) %*% pxnj) %*% t(pxnj) %*% xj/n)/(1- sum(coef_ast!=0)/n)
        }else{
            sigmaz_vec[j] = pp/(1- sum(coef_ast!=0)/n)
        }
        z[,j] = rnorm(n,0,sqrt(sigmaz_vec[j]))
    }


    w_vec = numeric(p)
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    para_list = foreach(z_idx_i = 1:p,.packages='glmnet') %dopar% {

        z_idx = z_idx_i
        z_indicator = numeric(p)
        z_indicator[z_idx]=1
        z_indicator = as.logical(z_indicator)
        x = scale(x)
        xnew = cbind(x[,z_indicator] + z[,z_indicator], x[,z_indicator] - z[,z_indicator], x[,as.logical(1-z_indicator)])

        fit1.lasso = glmnet(xnew, y, family=family, lambda=lambda_ast)
        z_idx_len = length(z_idx)
        b11 = coef(fit1.lasso, s = 'lambda.min')[2:(z_idx_len+1)]    
        b12 = coef(fit1.lasso, s = 'lambda.min')[(z_idx_len+2):(2*z_idx_len+1)]
        w =  abs(b11+b12)- abs(b11-b12)
        w
    }
    w_vec =  unlist(para_list)
    stopCluster(cl)

w_vec_nonzero = w_vec[which(w_vec!=0)]
w_idx = order(w_vec_nonzero, decreasing =T)
est_fdp = numeric(length(w_vec_nonzero))

for( i in 1:length(which(w_vec>0))){
    cut = w_vec_nonzero[w_idx[i]]
    est_fdp[i] = sum(w_vec_nonzero < -cut)/(i+1)
}

if(length(which(est_fdp-0.1>0))==0){
    gm_idx= length(which(w_vec>0))
}else{
    gm_idx = which(est_fdp-0.1>0)[1]-1
}

gm_selected = which(w_vec!=0)[w_idx[1:gm_idx]]

return(results= gm_selected, gm_statistics= w)

}
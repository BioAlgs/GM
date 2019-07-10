gm_var <- function( w_all, w_all_boot, alpha )
{
  ## This function calculate the decision based on GM and return the variance of the FDP
  ####################  
  ####################
  numrep = 1000 ## Number of simulation to get the estimation of the variance
  p <- dim(w_all_boot)[2]
  num.boot <- dim(w_all_boot)[1]
  
  w_all_sort = sort( w_all, decreasing=T )
  ord = order(w_all, decreasing=T)
  
  fdp.est = array(1, p)
  
  for(j in 1:p)
  {
    if( w_all_sort[j] > 0 )
    {
      thresh = w_all_sort[ j ]
      fdp.est[j] = sum( w_all < -thresh )/( ( sum( w_all >= thresh) + (sum(w_all>thresh)==0 ) ) )
    }
  }
  
  rej = max( ( fdp.est <= alpha ) *c(1:p) )
  
  if( rej > 0)
  {
    thresh = w_all_sort[ rej ]
    rej_ind = which( w_all >= thresh ) 
    var_prob = var_fdp( w_all_boot, thresh, numrep=numrep, alpha )
    variance =var_prob$var_fdp
  }else{
    thresh = w_all_sort[1]
    rej = 0
    rej_ind = c()
    variance = 0

  }
  list(rej=rej, rej_ind=rej_ind, variance=variance, thresh=thresh )  
}


GM_reject <- function( w_all, alpha )
{
  ## This function takes w and return and alpha and return a threshold for rejection
  
  w_all_sort = sort( w_all, decreasing=T )
  ord = order(w_all, decreasing=T)
  
  fdp.est = array(1, p)
  
  for(j in 1:p)
  {
    if( w_all_sort[j] > 0 )
    {
      thresh = w_all_sort[ j ]
      fdp.est[j] = sum( w_all < -thresh )/( ( sum( w_all >= thresh) + (sum(w_all>thresh)==0 ) ) )
    }
  }
  
  rej = max( ( fdp.est <= alpha ) *c(1:p) )
  
  if( rej > 0)
  {
    thresh = w_all_sort[ rej ]
    rej_ind = which( w_all >= thresh ) 
  }else{
    thresh =  max( w_all_sort[1], 5)
    rej = 0
    rej_ind = c()
  }
  list(rej=rej, rej_ind=rej_ind,  thresh=thresh )  
}


var_fdp <- function( w_all_boot, thresh, numrep=numrep, alpha=alpha )
{
  p <- dim(w_all_boot)[2]
  num.boot <- dim(w_all_boot)[1]
  fdp.sim = array(0, numrep)
  var_part = array(0, numrep)
  mean_part = array(0, numrep)
  for(ii in 1:numrep )
  {
    ind_sim = sample( c(1:num.boot), p, replace=T )
    w_j_sim = w_all_boot[ ind_sim + num.boot*c(0:(p-1)) ]
    fdp.sim[ii] = sum( w_j_sim < -thresh )/( sum( w_j_sim > thresh) + (sum(w_j_sim>thresh)==0 ) )
    evve.res = EV_VE( w_j_sim, w_all_boot, numrep=numrep, alpha=alpha )
    var_part[ii] = evve.res$var_part
    mean_part[ii] = evve.res$mean_part
  }
  list( var_fdp = mean( var_part) + var(mean_part) )
}

var_fdp_fixed_thresh <- function( w_all_boot, thresh, numrep=numrep, alpha=alpha )
{
  p <- dim(w_all_boot)[2]
  num.boot <- dim(w_all_boot)[1]
  fdp.sim = array(0, numrep)
  var_part = array(0, numrep)
  mean_part = array(0, numrep)
  for(ii in 1:numrep )
  {
    ind_sim = sample( c(1:num.boot), p, replace=T )
    w_j_sim = w_all_boot[ ind_sim + num.boot*c(0:(p-1)) ]
    fdp.sim[ii] = sum( w_j_sim < -thresh )/( sum( w_j_sim > thresh) + (sum(w_j_sim>thresh)==0 ) )
  }
  list( var_fdp = var(fdp.sim)  )
}

EV_VE <- function( w_j_sim, w_all_boot, numrep=numrep, alpha=alpha )
{
  gm.dec <- GM_reject( w_j_sim, alpha )
  thresh.evve = gm.dec$thresh
  p <- dim(w_all_boot)[2]
  num.boot <- dim(w_all_boot)[1]
  fdp.sim.evve = array(0, numrep)
  
  if( thresh.evve==0 ){
    fdp.sim.evve=array(0, numrep)  
  }else{
    ## Next we run the simulation with the boostrap sampling to get a sequence of fdp.sim
    for(iii in 1:numrep )
    {
      ind_sim.evve = sample( c(1:num.boot), p, replace=T )
      w_j_sim.evve = w_all_boot[ ind_sim.evve + num.boot*c(0:(p-1)) ]
      fdp.sim.evve[iii] = sum( w_j_sim.evve < -thresh.evve )/( sum( w_j_sim.evve > thresh.evve) + (sum(w_j_sim.evve >thresh.evve)==0 ) )
    }
  }
  list(var_part = var(fdp.sim.evve), mean_part = mean(fdp.sim.evve) )
}



GM_fix_thresh <- function( w_all, w_all_boot, alpha, thresh )
{
  numrep = 1000 ## Number of simulation to get the estimation of the variance
  p <- dim(w_all_boot)[2]
  num.boot <- dim(w_all_boot)[1]
  


  rej = sum( w_all > thresh )
  
  if( rej > 0)
  {
    rej_ind = which( w_all > thresh ) 
    var_prob = var_fdp_fixed_thresh( w_all_boot, thresh, numrep=numrep, alpha )
    variance =var_prob$var_fdp
    prob_fdp = 0
  }else{
    rej = 0
    rej_ind = c()
    variance = 0
    prob_fdp = 0
    
  }
  list(rej=rej, rej_ind=rej_ind, variance=variance, prob = prob_fdp, thresh=thresh )  
}

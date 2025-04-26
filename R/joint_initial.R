#' Initialise Joint Decomposition Parameters
#'
#' Generates starting values (*A*, shared‐component parameters,
#' modality-specific parameters and *S*) by ICA-based procedures and
#' eigen-decomposition heuristics.
#'
#' @inheritParams multi_view_decomposition
#' @param R  Optional preset vector of sub-rank values (default `NULL`,
#'           chosen automatically).
#' @param maxIter ICA iteration cap for `icaimax()` (default `100`).
#'
#' @return A `list` with elements `A`, `theta_common`, `theta_spe`, and `S`
#'         that feed directly into [joint_update_approx()].
#'
#' @export
joint_initial <- function(Y,q,q_common,V,rho=0.95,R = NULL,maxIter = 100)
{ 
  Y_data = Y[[1]]
  for (j in 2:length(q)){
      Y_data =rbind(Y_data, Y[[j]])
  }
  ICcorr = icaimax(t(Y_data),nc = q_common,center=FALSE,maxit=maxIter)
  
  ####################################Common components
  S_ini_common = list()
  for (j in 1:length(q)){
  S_ini_common[[j]]= matrix(0,ncol=dim(ICcorr$S)[1],nrow=q_common)
  }
  
  S_conn = ICcorr$S
  A = list()
  S = list()
  for (j in 1:length(q)){
  A[[j]] = Y[[j]]%*%S_conn%*%solve(t(S_conn)%*%S_conn)
  S[[j]] = solve(t(A[[j]])%*%A[[j]])%*%t(A[[j]])%*%Y[[j]]
  }
  theta_ini_common = list()
  for (j in 1:length(q)){
  for(i in 1:q_common)
  {
  theta_ini_common[[j]] = list()
  theta_ini_common[[j]][[i]] = list()
  }}
  
  for(i in 1:q_common)
  {
    S_mat = list()
    for (j in 1:length(q)){
    S_mat[[j]] = Ltrinv(S[[j]][i,],V,F)
    S_mat[[j]] = S_mat[[j]] + diag( rep(mean(S[[j]][i,]),V ))
    Sl = Ltrinv( S_conn[,i],V,FALSE)
    Sl = Sl + diag( rep(mean(S_conn[,i]),V ))
    eigenSl = eigen(Sl)
    orderEigen = order(abs(eigenSl$values),decreasing = TRUE)
      Rl = 2
      while( TRUE )
      {
        if (Rl >=10){
          rho1 = rho-0.05
        }
        else{rho1 = rho}
        eigenset = orderEigen[1:Rl]
        imgeRL = eigenSl$vectors[,eigenset]%*% diag(eigenSl$values[eigenset])%*% t(eigenSl$vectors[,eigenset])
        # image( imgeRL ) 
        if(cor(Ltrans(imgeRL,FALSE),ICcorr$S[,i]) > rho1) break
        Rl = Rl + 1
      }
   
  
    theta_ini_common[[j]][[i]]$J_l = matrix(0,ncol = V, nrow = Rl)
    for(k in 1:Rl)
    {
 
      theta_ini_common[[j]][[i]]$J_l[k,] = eigenSl$vectors[,eigenset[k]]
    }
    sign = list()
    sign = sign(diag(theta_ini_common[[j]][[i]]$J_l%*%S_mat[[j]] %*% t(theta_ini_common[[j]][[i]]$J_l)))
    theta_ini_common[[j]][[i]]$lam_l = sign*sqrt(apply((S_mat[[j]] %*% t(theta_ini_common[[j]][[i]]$J_l))^2,2,sum))
    if( theta_ini_common[[j]][[i]]$lam_l[1]<0 ){theta_ini_common[[j]][[i]]$lam_l = -1*theta_ini_common[[j]][[i]]$lam_l}
    S_ini_common[[j]][i,] = Ltrans( t(theta_ini_common[[j]][[i]]$J_l)%*%diag(theta_ini_common[[j]][[i]]$lam_l)%*%theta_ini_common[[j]][[i]]$J_l,FALSE)
    }
  }
  ####################################Residual components
  theta_ini = list()
  for (j in 1:length(q)){
    theta_ini[[j]] = list()
    for(i in 1:q[j])
    {
  theta_ini[[j]][[i]] = list()
  }}
  
  projection = A[[1]]%*%S_ini_common[[1]]
  for (j in 2:length(q)){
      projection = rbind(projection, A[[j]]%*%S_ini_common[[j]])
   }
  if(q_common==1){
    residual = Y_data - projection
  }else{
  residual = Y_data - projection}
  
  theta_ini = list()
  for (j in 1:length(q)){
    theta_ini[[j]] = list()
    for(i in 1:q[j]){
      theta_ini[[j]][[i]]= list()
    }}
    S_ini = list()
  for (j in 1: length(q)){
  S_ini[[j]] = matrix(0,ncol=dim(ICcorr$S)[1],nrow=q[j])

  if(q[j]!=0){
    ICcorr = icaimax(t(residual[(1:nrow(Y[[j]])),]), nc = q[j], center=FALSE,maxit=maxIter)
    
  for(i in 1:q[j])
  {
    Sl = Ltrinv( ICcorr$S[,i],V,FALSE)
    Sl = Sl + diag( rep(mean(ICcorr$S[,i]),V ))
    eigenSl = eigen(Sl)
    orderEigen = order(abs(eigenSl$values),decreasing = TRUE)
    {
      Rl = 2
      while( TRUE )
      {
        eigenset = orderEigen[1:Rl]
        imgeRL = eigenSl$vectors[,eigenset]%*% diag(eigenSl$values[eigenset])%*% t(eigenSl$vectors[,eigenset])
        # image( imgeRL ) 
        if(cor(Ltrans(imgeRL,FALSE),ICcorr$S[,i]) > rho) break
        Rl = Rl + 1
      }
    }
    
    theta_ini[[j]][[i]]$lam_l = eigenSl$values[eigenset]
    if( theta_ini[[j]][[i]]$lam_l[1]<0 ){theta_ini[[j]][[i]]$lam_l = -1*theta_ini[[j]][[i]]$lam_l}
    theta_ini[[j]][[i]]$J_l = matrix(0,ncol = V, nrow = Rl)
    for(k in 1:Rl)
    {
      theta_ini[[j]][[i]]$J_l[k,] = eigenSl$vectors[,eigenset[k]]
    }
    S_ini[[j]][i,] = Ltrans( t(theta_ini[[j]][[i]]$J_l)%*%diag(theta_ini[[j]][[i]]$lam_l)%*%theta_ini[[j]][[i]]$J_l,FALSE)
  }}
  }
  
  S_ini_comb = list()
  A_ini = list()
for (j in 1:length(q)){
    S_ini_comb[[j]] = rbind(S_ini_common[[j]],S_ini[[j]])
    A_ini[[j]] = Y[[j]]%*%t(S_ini_comb[[j]])%*%solve(S_ini_comb[[j]]%*%t(S_ini_comb[[j]]))
  }
  
  
  # scale up
# for(l in 1:q){
#     # unit norm each column of A
#     # scaleL = sqrt(sum(A_ini[,l]^2))
#     scaleL = sd(A_ini[,l])
#     A_ini[,l] = A_ini[,l] / scaleL
#     S_ini[l,] = S_ini[l,] * scaleL
#     # scale X_l correspondingly
#     theta_ini[[l]]$X_l = theta_ini[[l]]$X_l * sqrt(scaleL)
#   }
  
  # since after preprocessing, A_tilde is orthogonal
  # compute A_tilde transpose/inverse (g-inverse of A-tilde)
  # if X has full column rank, (X'X)^(-1)X' is its g-inverse
  # why use g-inverse here: since A_ini is N*q
  # afterwards, in the update process, we actually can just use A_tilde transpose(after scale), or g-inverse
  M_ini = list()
for (j in 1:length(q)){
    M_ini[[j]] =  solve(t(A_ini[[j]])%*%A_ini[[j]])%*%t(A_ini[[j]]) # g-inverse of A
  for(l in 1:q_common){
    theta_ini_common[[j]][[l]]$M_l = M_ini[[j]][l,]
  }
    if (q[j]!=0){
  for(l in 1:q[j]){
      theta_ini[[j]][[l]]$M_l = M_ini[[j]][l,]
  }}
    }
  return(list(A=A_ini,theta_common = theta_ini_common,theta_spe = theta_ini, S = S_ini_comb))
}


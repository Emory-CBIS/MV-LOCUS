#' Joint Low-Rank Decomposition for Multiple Connectivity Modalities
#'
#' Implements the full iterative algorithm that extracts *q_common* joint
#' subnetworks and view-specific subnetworks (*q*) from a list of subject-by-edge
#' matrices.  Supports optional SCAD, L1, or hard-threshold sparsity penalties
#' and convergence checks on both the mixing matrices (*A*) and component
#' matrices (*S*).
#'
#' @param Y          `list` of length *K*; each element is an `N × p` matrix
#'                   (subjects × vectorised upper-triangular connectomes).
#' @param q_common   `integer(1)`. Number of shared subnetworks.
#' @param q          `integer` vector of length *M*. Numbers of view-specific
#'                     subnetworks.
#' @param V          `integer(1)`. Number of brain nodes (`K = V (V − 1)/2`).
#' @param MaxIteration Maximum outer iterations (default `5000`).
#' @param penalty    `"SCAD"`, `"L1"`, `"Hardthreshold"`, or `NULL`.
#' @param phi,psi,gamma,rho Numeric tuning parameters (see paper /
#'                   function details).
#' @param espli1,espli2 Absolute tolerance for change in *A* and *S*.
#' @param silent     Logical; suppress progress output if `TRUE`.
#'
#' @return A named `list` with components
#' * `A` – list of mixing matrices (one per modality)
#' * `S` – list of component loadings
#' * `S_sparse` – thresholded version of `S`
#' * `theta_common`, `theta_spe` – lists of parameter objects
#' * `Conver` – logical flag indicating convergence
#'
#' @examples
#' ## Not run:
#' # res <- joint_decomposition_multi(
#' #   Y            = abcd,
#' #   q_common     = 30,
#' #   q            = c(1, 1),
#' #   V            = 360,
#' #   penalty      = "SCAD",
#' #   phi          = 4.5,
#' #   psi          = 10,
#' #
#' # )
#' ## End(Not run)
#'
#' @seealso [joint_initial()], [joint_update_approx()]
#' @export

multi_view_decomposition <- function(Y, q, q_common, V, MaxIteration=5000, penalty="SCAD", phi = 0.9, psi = 1, gamma =3,
                  espli1=5e-4, espli2=5e-4, rho=0.95, silent=FALSE)
{
  eigen_cor    = 0
  # demean the data
  for (i in 1:length(q)){
  Y[[i]] = sweep(Y[[i]],2,apply(Y[[i]],2,mean),"-") }
  # preprocess the data if True'
  prec = list(); H_inv = list()
  for (i in 1:length(q)){
   prec[[i]] = Locus_preprocess(Y[[i]],q[i]+q_common)
    }

  for  (i in 1:length(q)){
    Y[[i]] = prec[[i]]$Y_new; H_inv[[i]] = prec[[i]]$H_inv}

  K = dim(Y[[1]])[2]          # Number of edges
  if(V != (sqrt(1+8*K)+1)/2)
  {
    stop("V is not correctly specified! Please double check the dimension of your input data.")
  }
          # Number of ICs

  # Initial Estimation
  theta_ini = joint_initial(Y,q,q_common,V,rho=rho)
  A = theta_ini$A;  S = theta_ini$S
  theta_common = theta_ini$theta_common
  theta_spe = theta_ini$theta_spe

  # Update Parameters
  Iter = 1
  while(Iter <= MaxIteration)
  {
      theta_new = joint_update_approx(Y,A,theta_common,theta_spe, q, q_common, eigen_cor = eigen_cor, psi = psi, penalt= penalty, lambda_ch = phi, gamma = gamma ,imput_method = "Previous",silent = silent,H_inv= H_inv, Iter)
    # else
    # {
    #   theta_new = Locus_update(Y,A,theta,penalt= penalty,lambda_ch = phi, gamma = 2.1,silent = silent)
    # }
    A_new = list()
    S_new = list()
    theta_spe_new = list()
    theta_common_new = list()
    for (j in 1:length(q)){
    S_new[[j]] = theta_new$S[[j]]/sd((theta_new$S[[j]]))*sd(S[[j]]) ;theta_common_new[[j]]  = theta_new$theta_common[[j]]
    theta_spe_new[[j]]  = theta_new$theta_spe[[j]]
    A_new[[j]] = theta_new$A[[j]]
    }
    errS = 0
    errA = 0
    for (j in 1:length(q)){
      errS = errS +norm(as.matrix(S_new[[j]]-S[[j]]))/norm(as.matrix(S[[j]]))
      errA = errA + norm(as.matrix(A_new[[j]]-A[[j]]))/norm(as.matrix(A[[j]]))
    }

    if(sum(is.na(c(errS,errA)))>0){return(list(Conver=FALSE))}

    if(!silent)
    {
      message(paste("Iter ",Iter,"; Percentage change on indpendent components: " , round(errS,3),"; Percentage change on mixing matrices: ",round(errA,3),".",sep=""))
    }

    #+matrix(rnorm(dimB[1]*dimB[2],0,1e-5),nrow = dimB[1],ncol = dimB[2])
    #

    A = A_new  ; S= S_new; theta_common = theta_common_new
    theta_spe = theta_spe_new

    A_final = list()
    if(errA < espli1 & errS < espli2)
    {
      if(!silent){cat("Converaged!")}
    for (j in 1:length(q)){
      A_final[[j]] = H_inv[[j]] %*% A_new[[j]]
      }
      return(list(Conver = TRUE, A=A_final,  S = S, S_sparse=theta_new$S_sparse, theta_common=theta_common_new, theta_spe = theta_spe_new
                  ))
    }
    Iter = Iter + 1
  }

  if(!silent){cat("Failed to converge!")}
A_final = list()
for (j in 1:length(q)){
  A_final[[j]] = H_inv[[j]] %*% A_new[[j]]
}
  return(list(Conver = TRUE, A=A_final,  S = S, S_sparse=theta_new$S_sparse, theta_common=theta_common_new, theta_spe = theta_spe_new
  ))
}

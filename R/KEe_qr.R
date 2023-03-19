r_tau<-function(tau){

  return(log(1/2)/log(1 - tau))

}


J_tau<-function(tau,x){
  if (0 < tau & tau <= 0.5){

    return(r_tau(tau) * (1 - x)^(r_tau(tau) - 1))

  }
  else{

    return(r_tau(1 - tau) * (x)^(r_tau(1 - tau) - 1))

  }

}


kernel_qr <- function(x, y, kernel, C, tau){

  n <- nrow(x)
  K <- kernelMatrix(kernel, x)
  a <- rep(0, n)
  A <- cbind(diag(n), -diag(n), rep(1,n))
  b <- c(rep(C*(tau-1),n), rep(-C*tau, n),0)

  alpha<-solve.QP(K, y, A, b)$solution
  q_hat <- K%*%alpha

  return(q_hat)
}

CDF_QR <- function(x, y, kernel, C, s){

  n <- nrow(x)
  q_hat <- matrix(0, n, s-1)

  for(i in 1:(s-1)){
    tau_s <- i/30
    q_hat_s <- kernel_qr(x, y, kernel, C, tau_s)
    q_hat[, i] <- q_hat_s
  }

  CDF_hat <- rep(0, n)
  for(i in 1:n){
    CDF_hat[i] <- which.min(abs(q_hat[i,]-y[i])) / s
  }

  return(CDF_hat)
}

OKE_qr <- function(x, y, lambda, kernel, C, s, tau, CDF_hat){

  n <- nrow(x)
  K <- kernelMatrix(kernel, x)

  W <- matrix(0, n, n)
  for(i in 1:n){
    W[i, i] = J_tau(tau, CDF_hat[i])
  }

  alpha = solve(K%*%W%*%K + lambda * n * K) %*% (K %*% W %*% y)

  return(alpha)
}



KEE_qr <- function(x, y, lambda, kernel, C, s, k_1, k_n, tau, CDF_hat){

  n <- nrow(x)
  K <- kernelMatrix(kernel, x)

  xi_tau_hat=matrix(0, n, (k_n - k_1 + 1))
  for(u in (n-k_n):(n-k_1)){
    tau_u <- u/(n+1)
    alpha <- OKE_qr(x, y, lambda, kernel, C, s, tau_u, CDF_hat)

    xi_tau_hat_u <- K%*%alpha
    xi_tau_hat[, u - n + k_n + 1] <- xi_tau_hat_u
  }

  gamma_hat <- rowSums(log(xi_tau_hat)) / (k_n - k_1) - log(xi_tau_hat[,1]) * (k_n - k_1 + 1) / (k_n - k_1)

  xi_tau_extreme <- rep(0, n)
  for(j in 1:n){
    xi_tau_extreme[j] = (((1 + n- k_n) / (n + 1))/(1 - tau))^(gamma_hat[j]) * xi_tau_hat[j, 1]
  }

  return(xi_tau_extreme)
}






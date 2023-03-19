solve.QP.compact <- function(Dmat, dvec, Amat, Aind, bvec, meq=0,
                             factorized=FALSE){  
  n     <- nrow(Dmat)
  q     <- ncol(Amat)
  anrow <- nrow(Amat)
  if( missing(bvec) )
    bvec <- rep(0,q)

  if( n != ncol(Dmat) )
    stop("Dmat is not symmetric!")
  if( n != length(dvec) )
    stop("Dmat and dvec are incompatible!")
  if( (anrow+1 != nrow(Aind)) || (q != ncol(Aind)) || (q != length(bvec)) )
    stop("Amat, Aind and bvec are incompatible!")
  Aindok <- .Fortran(.QP_aind,
                     as.integer(Aind), as.integer(anrow+1),
                     as.integer(q), as.integer(n),
                     ok=TRUE)$ok
  if( !Aindok )
    stop("Aind contains illegal indexes!")
  if( (meq > q) || (meq < 0 ) )
    stop("Value of meq is invalid!")

  iact  <- rep(0,q)
  nact  <- 0
  r     <- min(n,q)
  sol   <- rep(0,n)
  lagr  <- rep(0,q)
  crval <- 0
  work  <- rep(0,2*n+r*(r+5)/2+2*q+1)
  iter  <- rep(0,2)

  res1 <- .Fortran(.QP_qpgen1,
                   as.double(Dmat), dvec=as.double(dvec),
                   as.integer(n), as.integer(n),
                   sol=as.double(sol), lagr=as.double(lagr),
                   crval=as.double(crval),
                   as.double(Amat), as.integer(Aind), as.double(bvec),
                   as.integer(anrow), as.integer(q), as.integer(meq),
                   iact=as.integer(iact), nact=as.integer(nact),
                   iter=as.integer(iter), 
                   work=as.double(work), ierr=as.integer(factorized))

  if( res1$ierr == 1)
    stop("constraints are inconsistent, no solution!")
  else if( res1$ierr == 2)
    stop("matrix D in quadratic function is not positive definite!")
    
  list(solution=res1$sol,
       value=res1$crval,
       unconstrained.solution=res1$dvec,
       iterations=res1$iter,
       Lagrangian = res1$lagr,
       iact=res1$iact[1:res1$nact])   
}


solve.QP <- function(Dmat, dvec, Amat, bvec, meq=0, factorized=FALSE){  

  n     <- nrow(Dmat)
  q     <- ncol(Amat)
  if( missing(bvec) )
    bvec <- rep(0,q)

  if( n != ncol(Dmat) )
    stop("Dmat is not symmetric!")
  if( n != length(dvec) )
    stop("Dmat and dvec are incompatible!")
  if( n != nrow(Amat) )
    stop("Amat and dvec are incompatible!")
  if( q != length(bvec) )
    stop("Amat and bvec are incompatible!")
  if( (meq > q) || (meq < 0 ) )
    stop("Value of meq is invalid!")
  
  iact  <- rep(0,q)
  nact  <- 0
  r     <- min(n,q)
  sol   <- rep(0,n)
  lagr  <- rep(0,q)
  crval <- 0
  work  <- rep(0,2*n+r*(r+5)/2+2*q+1)
  iter  <- rep(0,2)

  res1 <- .Fortran(.QP_qpgen2,
                   as.double(Dmat), dvec=as.double(dvec),
                   as.integer(n), as.integer(n),
                   sol=as.double(sol), lagr=as.double(lagr),
                   crval=as.double(crval),
                   as.double(Amat), as.double(bvec), as.integer(n),
                   as.integer(q), as.integer(meq),
                   iact=as.integer(iact), nact=as.integer(nact),
                   iter=as.integer(iter), work=as.double(work),
                   ierr=as.integer(factorized))

  if( res1$ierr == 1)
    stop("constraints are inconsistent, no solution!")
  else if( res1$ierr == 2)
    stop("matrix D in quadratic function is not positive definite!")

  list(solution=res1$sol,
       value=res1$crval,
       unconstrained.solution=res1$dvec,
       iterations=res1$iter,
       Lagrangian = res1$lagr,
       iact=res1$iact[1:res1$nact])   
}

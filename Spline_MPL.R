library("R.matlab")
library(Matrix)
library(survivalMPLdc)
library(splines2)
library(readr)
path <- ('D:/科研项目/Informative censoring/simulation revision/Spline/mat')

pathname1 <- file.path(path,'T.mat')
pathname2 <- file.path(path,'Z.mat')
pathname3 <- file.path(path,'Delta.mat')
pathname4 <- file.path(path,'Beta0.mat')
pathname5 <- file.path(path,'Phi0.mat')

T <- data.frame(readMat(pathname1))
cova <- data.frame(readMat(pathname2))
cova <- data.matrix(cova)


Delta <- data.frame(readMat(pathname3))
eta <- 1-Delta


beta0 <- data.frame(readMat(pathname4))
phi0 <- data.frame(readMat(pathname5))

beta0 <- as.numeric(beta0[,1])
phi0 <- as.numeric(phi0[,1])

# Discretizing a sample of n survival times into sub-intervals with no. of 'binCount' subjects in each sub-interval
discrBinNA <- function( survdat, binCount, tie )
{
  X <- survdat[ , 1 ] # survival times
  len <- length( X )  # sample size
  binid <- seq( 1, len, binCount )
  nbins <- length( binid )
  sortX <- sort( X, index.return = TRUE )       # sort the data into ascending order
  sX <- sortX$x
  binedg <- sX[ binid ]
  
  if( tie == 'No' )
  {
    
    if( binCount == 1 )
    {
      binwv <- binedg[ 1:nbins ] - c( 0, binedg[ 1:( nbins - 1 ) ] ) #binwv: bin widths
      binedg <- c( 0, binedg )
      binID <- rank( X ) #binID: bin IDs
    }
    else
    {
      
      i <- 0
      while( i < ( nbins - 1 ) )
      {
        
        i <- i + 1
        ntied <- sum( binedg[ i ] == binedg )
        if( ( i + ntied ) <= nbins )
        {
          binedg <- c( binedg[ 1:i ], binedg[ ( i + ntied ):nbins ] )
        }
        else
        {
          binedg <- c( binedg[ 1:i ] )
        }
        nbins <- length( binedg )
        
      }
      
      if( binedg[ nbins ] == max( X ) )
      {
        nbins <- nbins - 1
        binedg <- binedg[ 1:nbins ]
      }
      binedg[ 1 ] <- min( X )
      binedg[ nbins + 1 ] <- max( X ) + ( 1e-2 )*( binedg[ nbins ] - binedg[ nbins - 1 ] )
      
      binedg[ nbins ] <- binedg[ nbins ]*( 1 - ( 1e-5 ) )
      
      binwv <- binedg[ 2:( nbins + 1 ) ] - binedg[ 1:nbins ]
      binID <- rep( 0, len )
      for ( i in 1:len )
      {
        for ( j in 1:( nbins ) )
        {
          if( X[ i ] < binedg[ j + 1 ] & X[ i ] >= binedg[ j ] )
          {
            binID[ i ] <- j
          }
        }
      }
      
    }
  }
  else
  {
    i <- 0
    while( i < ( nbins - 1 ) )
    {
      i <- i + 1
      ntied <- sum( binedg[ i ] == binedg )
      if( ( i + ntied ) <= nbins )
      {
        binedg <- c( binedg[ 1:i ], binedg[ ( i + ntied ):nbins ] )
      }
      else
      {
        binedg <- c( binedg[ 1:i ] )
      }
      nbins <- length( binedg )
      
    }
    if( binedg[ nbins ] == max( X ) )
    {
      nbins <- nbins - 1
      binedg <- binedg[ 1:nbins ]
    }
    binedg[ 1 ] <- min( X )
    binedg[ nbins + 1 ] <- max( X ) + ( 1e-2 )
    binedg[ nbins ] <- binedg[ nbins ]*( 1 - ( 1e-5 ) )
    binwv <- binedg[ 2:( nbins + 1 ) ] - binedg[ 1:nbins ]
    binID <- rep( 0, len )
    for ( i in 1:len )
    {
      for ( j in 1:( nbins ) )
      {
        if( X[ i ] < binedg[ j + 1 ] & X[ i ] >= binedg[ j ] )
        {
          binID[ i ] <- j
        }
      }
    }
  }
  classify <- cbind( binwv )
  return( list( discretize = classify, ID = binID, binedg = binedg ) )
}

# Initial estimates of theta (piecewise constant estimate of h_{0t}) based on independent censoring assumption
theta_initial <- function( del, psix, Psix, beta0, cova )
{
  n <- dim( psix )[ 1 ]
  m <- dim( psix )[ 2 ]
  eregt <- exp( cova%*%beta0 )
  theta <- colSums( matrix( del, n, m )*psix )/( colSums( matrix( eregt, n, m )*Psix ) + 1e-5 )
  return( theta )
}

# Initial estimates of gamma (piecewise constant estimate of h_{0c}) based on independent censoring assumption
gamma_initial <- function( eta, psix, Psix, phi0, cova )
{
  n <- dim( psix )[ 1 ]
  m <- dim( psix )[ 2 ]
  eregc <- exp( cova%*%phi0 )
  gamma <- colSums( matrix( eta, n, m )*psix )/( colSums( matrix( eregc, n, m )*Psix ) + 1e-5 )
  return( gamma )
}

# Penality functions

#For piecewise constant

# Second order difference
mat2 <- function( psix, X, eps )
{
  m <- dim( psix )[2]
  spsix <- psix[ sort( X, index.return = T )$ix, ]
  R <- t( diff( spsix, lag = 1, difference = 2 ) )%*%diff( spsix, lag = 1, difference = 2 )
  diag( R ) <- diag( R ) + eps
  return( R )
}

# First order difference
#For piecewise constant
mat1 <- function(psix, X, eps)
{
  m<-dim(psix)[2]
  spsix<-psix[sort(X, index.return=T)$ix,]
  R<-t(diff(spsix, lag=1, difference=1))%*%diff(spsix, lag=1, difference=1)
  diag(R)<-diag(R)+eps
  return(R)
}

#For m-spline
penalty_msp<-function(numSp, ordSp, IntKnt, bryKnt)
{
  R<-matrix(0, nrow=numSp, ncol=numSp)
  xknots <- c(rep(min(bryKnt), ordSp), IntKnt, rep(max(bryKnt), ordSp))
  for (ii in 1:numSp)
  {
    for (jj in ii:numSp){
      if (jj - ii<ordSp){
        kntset <- xknots[xknots>=xknots[jj] & xknots<=xknots[ii+ordSp]];
        kntsum <- 0;
        for (kk in 1:(length(kntset)-1))
        {
          kntsum <- kntsum + 
            mSpline(kntset[kk], knots=IntKnt, intercept = TRUE,degree=ordSp-1, Boundary.knots=bryKnt, 
                    derivs=ordSp-2)[ii]*
            mSpline(kntset[kk], knots=IntKnt, intercept = TRUE, degree=ordSp-1, 
                    Boundary.knots=bryKnt, derivs=ordSp-2)[jj]*(kntset[kk+1]-kntset[kk]);
        }
        R[ii, jj] <- kntsum;
      }
    }
  }
  R[lower.tri(R, diag = FALSE)] <- t(R)[lower.tri(R, diag = FALSE)]
  return(R)
}

discrsurvdat <- discrBinNA(cbind(T,Delta,eta),10,'NO')

groupsurvdat <- discrsurvdat$discretize
binwv <- groupsurvdat[,1] #the bin widths
ID <- discrsurvdat$ID  #the bin ID of each subject
binedg <- discrsurvdat$binedg
n <- nrow(cova)       #sample size
p <- ncol(cova)#number of covariates for each subject

numIntKnt <- length(binedg)-2
IntKnt <- binedg[2:(numIntKnt+1)]
bryKnt<-c(0, max(binedg))
ordSp <- 4; # 4-spline;1-piecewise constant
numSp <- numIntKnt+ordSp
m <- numSp

psix<-mSpline(T[,1], knots = IntKnt, degree = ordSp-1, intercept = TRUE, derivs=0, Boundary.knots=bryKnt )
Psix<-iSpline(T[,1], knots = IntKnt, degree = ordSp-1, intercept = TRUE, derivs=0, Boundary.knots=bryKnt )

RT <- penalty_msp(numSp, ordSp, IntKnt, bryKnt)

theta_ini <- theta_initial( del, psix, Psix, beta0, cova )
gamma_ini <- gamma_initial( del, psix, Psix, phi0, cova )



setwd("D:/科研项目/Informative censoring/simulation revision/Spline")
#write.csv(discrsurvdat,'discrsurvdat.csv')
write.csv(theta_ini,'theta_ini.csv')
write.csv(gamma_ini,'gamma_ini.csv')
write.csv(Psix,'Psi.csv')
write.csv(psix,'psix.csv')
write.csv(RT,'RT.csv')

write.table(Psix, "Psix.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

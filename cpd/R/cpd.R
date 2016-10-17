#' @import matrixcalc
NULL

#' repmat
#'
#' A function that is meant to recapitulate the repmat function of matlab.
#' This function is used here to make replicated versions of two point sets
#' so that 2D matrix math can be used to calculate the difference between
#' all possible combinations of points in the two point sets.
#'
#' @param x matrix to repeat/replicate
#' @param m integer number of rows of replication
#' @param n integer number of cols of replication
#' @param loopRows boolean indicating whether the each row should represent a looping sequence of the rows of x (loopRows==T produces a matrix that starts with x[1,], then x[2,]... while loopRows==F produces a matrix that starts with x[1,], then x[1,] repeated m times before switchin to x[2,])
#' @param loopCols boolean indicating whether the each col should represent a looping sequence of the cols of x (loopCols==T produces a matrix that starts with x[,1], then x[,2]... while loopRows==F produces a matrix that starts with x[,1], then x[,1] repeated n times before switchin to x[,2])
#'
#' @return matrix made by 'copying' the x matrix m-by-n times
#' @export
repmat <- function(x, m, n, loopRows=T, loopCols=T)
{
     if(loopRows & loopCols)
     {
          return(x[rep(1:nrow(x), times=m), rep(1:ncol(x), times=n)])
     }
     else if(loopRows & !loopCols)
     {
          return(x[rep(1:nrow(x), times=m), rep(1:ncol(x), each=n)])
     }
     else if(!loopRows & loopCols)
     {
          return(x[rep(1:nrow(x), each=m), rep(1:ncol(x), times=n)])
     }
     else
     {
          return(x[rep(1:nrow(x), each=m), rep(1:ncol(x), each=n)])
     }
}

#' get1(x, flip=F)
#'
#' Get a 1 column matrix of 1's based upon the dimensions of x. If the parameter
#' flip==F, then the resultant matrix has nrow(x) elements. If flip==T, then
#' the resultant matrix has ncol(x) elements.
#'
#' @param x the matrix to use for determining the dimension of the resultant matrix
#' @param flip boolean indicating which dimension to use to determine the number of elements in the resultant column matrix of 1's
#'
#' @export
get1 <- function(x, flip=F)
{
     if(flip)
     {
          matrix(rep(1,ncol(x)), ncol=1)
     }
     else
     {
          matrix(rep(1,nrow(x)), ncol=1)
     }
}

#' \%=\%
#'
#' Internal interface for the %=% assignment. This will be used to enable
#' Matlab-like assignment of variables from functions that return lists.
#'
#' E.g.,
#'   Matlab: [a, b] = dim(A)
#'   R: l(a, b) \%=\% dim(A)
#'
#' @param l left hand side of the assignment
#' @param r right hand side of the assignment
#'
#' @rdname equals
#'
#' @export
'%=%' <- function(l, r)
{
     UseMethod('%=%')
}

#' \%=\%.lbunch
#'
#' Internal function will be used to enable
#' Matlab-like assignment of variables from functions that return lists.
#'
#' E.g.,
#'
#'   Matlab: [m, n] = dim(A)
#'   R: l(m, n) \%=\% dim(A)
#'
#' @param l left hand side of the assignment
#' @param r right hand side of the assignment
#'
#' @rdname lbunch
#'
#' @export
#'
#' @examples A <- matrix(1:4, ncol=2); l(m, n) \%=\% dim(A);
'%=%.lbunch' <- function(l, r)
{
     Names = lapply(l, as.character)
     Envir = as.environment(-1)

     for (II in 1:length(Names)) {
          Name = Names[[II]]
          assign(Name, r[[II]], pos=Envir)
     }
}

#' l
#'
#' Internal function used with %=% to perform Matlab-like assignment
#' of variables from functions that return a list.
#'
#' @param ... variable to be gathered ans assigned in the list
#'
#' @export
l <- function(...)
{
     List = as.list(substitute(list(...)))[-1L]
     class(List) = 'lbunch'
     List
}

#' transformY(Y, B, tr)
#'
#' Function for transforming a matrix of points by a transform matrix
#' that performs scaling and/or rotation. Typically 2X2 for xy points
#'
#' @param Y matrix representing M points and D dimensions. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param B Transformation matrix (scaling and/or rotation)
#' @param tr single column matrix with the translation that will be applied for each dimension of the points in Y
#'
#' @return YT the transformed version of the points in Y
#' @export
transformY <- function(Y, B=diag(2), tr=matrix(c(0,0)))
{
     return(Y %*% t(B) + repmat(t(tr), m=nrow(Y)))
}

#' getDiffSquared(X, Y, B, tr)
#'
#' This function calculates the MXN matrix of L2 norms between the points in X
#' and the points in Y AFTER transformation with B and tr (i.e., || X - (BY+tr) ||^2)
#'
#' @param X matrix representing N points and D dimensions. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param Y matrix representing M points and D dimensions. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param B Transformation matrix (scaling and/or rotation)
#' @param tr single column matrix with the translation that will be applied for each dimension of the points in Y
#'
#' @export
getDiffSquared <- function(X, Y, B, tr)
{
     D <- ncol(X)
     M <- nrow(Y)
     N <- nrow(X)

     Y_t <- transformY(Y=Y, B=B, tr=tr)
     Ybig_t <- repmat(Y_t, m=N, loopRows=T) # index m
     Xbig <- repmat(X, m=M, loopRows=F) # index n
     temp <- matrix(rowSums((Xbig - Ybig_t)^2), ncol=N)
}

#' getP(X, Y, sigma2, w, B, tr)
#'
#' This function uses the function 'getDiffSquared' to estimate Pmn which is defined
#' in the associated manuscript for this package in Figure 3 depicting the
#' Coherent Point Drift algorithm. This implementation is based upon the paper by Andriy
#' Myronenko and Xubo Song titled "Point set registration: coherent point drift" in IEEE
#' Trans Pattern Anal Mach Intell, Dec 2010.
#'
#' @param X matrix representing N points and D dimensions. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param Y matrix representing M points and D dimensions. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param sigma2 double indicating the variance around each point in the GMM representation of the registration problem
#' @param w double value (0 \<= w \< 1) for accomodating noise and outliers
#' @param B Transformation matrix (scaling and/or rotation)
#' @param tr single column matrix with the translation that will be applied for each dimension of the points in Y
#'
#' @export
getP <- function(X, Y, sigma2, w, B, tr)
{
     D <- ncol(X)
     M <- nrow(Y)
     N <- nrow(X)

     temp <- getDiffSquared(X=X, Y=Y, B=B, tr=tr)
     ret <- exp( (-1/sigma2) * temp )
     colTots <- matrix(colSums(ret), ncol=N)
     colTots <- repmat(colTots, m=M)
     alpha <- ((2*pi*sigma2)^(D/2))*(w/(1-w))*(M/N)

     return(ret / (colTots + alpha))
}

#' getQ(P, X, Y, B, tr, sigma2)
#'
#' Calculate the objective function of the Coherent Point Drift algorithm
#'
#' @param P The M x N matrix of column normalized probabilities
#' @param X matrix representing N points and D dimensions. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param Y matrix representing M points and D dimensions. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param B Transformation matrix (scaling and/or rotation)
#' @param tr single column matrix with the translation that will be applied for each dimension of the points in Y
#' @param sigma2 double indicating the variance around each point in the GMM representation of the registration problem
#'
#' @return double value to be minimized during EM steps
#' @export
getQ <- function(P, X, Y, B, tr, sigma2)
{
     diffSquared <- getDiffSquared(X=X, Y=Y, B=B, tr=tr)
     ret <- (1/(2*sigma2))*sum(P * diffSquared)
     return(ret)
}

#' getRS(B)
#'
#' This function is used to calculate the x and y scaling and rotation
#' encoded in the transformation matrix B.
#'
#' @param B Transformation matrix (scaling and/or rotation)
#'
#' @return list(Sx, Sy, Theta) where Sx and Sy are the scaling in the x and y direction while Theta is the rotation.
#' @export
getRS <- function(B)
{
     Sx <- sign(B[1,1])*sqrt(B[1,1]^2+B[1,2]^2)
     Sy <- sign(B[2,2])*sqrt(B[2,1]^2+B[2,2]^2)
     Theta <- atan2(B[2,1], B[2,2]) # or atan2(-B[1,2], B[1,1])
     return(list(Sx=Sx, Sy=Sy, Theta=Theta))
}

#' makeR(Theta)
#'
#' This function creates a transformation matrix that will produce
#' a rotation of Theta radians in point sets to which it is applied.
#'
#' @param Theta double value representing the rotation angle in radians
#'
#' @export
makeR <- function(Theta)
{
     return(matrix(c(cos(Theta), -sin(Theta), sin(Theta), cos(Theta)), ncol=2, byrow=T))
}

#' makeB(Sx, Sy, Theta)
#'
#' This function creates a D x D transformation matrix that will produce
#' a rotation of Theta radians and x-y scaling of Sx and Sy in point
#' sets (number of dimensions = D) to which it is applied.
#'
#' @param Sx double x scaling factor
#' @param Sy double y scaling factor
#' @param Theta double value representing the rotation angle in radians
#'
#' @return D x D transformation matrix
#' @export
makeB <- function(Sx, Sy, Theta)
{
     a11 <- Sx*cos(-Theta)
     a12 <- Sx*sin(-Theta)
     a21 <- -Sy*sin(-Theta)
     a22 <- Sy*cos(-Theta)
     return( matrix(c(a11, a12, a21, a22), ncol=2, byrow=T))
}

#' plotXYYT(X, Y, YT, includeYT=T)
#'
#' Function for rudimentary plot showing the target point set, X (black),
#' and the initial location of the point set to be registered, Y (green).
#' If includeYT==T, then the location of the transformed point set, YT (red).
#'
#' @param X matrix representing N points and D dimensions. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param Y matrix representing M points and D dimensions. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param YT matrix representing the transformed version of Y
#' @param includeYT boolean indicating whether to plot YT
#'
#' @export
plotXYYT <- function(X, Y, YT, includeYT=T)
{
     xlim <- range(X[,1], Y[,1], YT[,1])
     ylim <- range(X[,2], Y[,2], YT[,2])
     plot(X, xlim=xlim, ylim=ylim)
     points(Y, col='green')
     if(includeYT)
     {
          points(YT, col='red', cex=0.7)
     }
}

#' Coherent Point Drift Registration (Affine Transformations)
#'
#' This implementation is based upon the paper by Andriy Myronenko and Xubo Song
#' titled "Point set registration: coherent point drift" in IEEE Trans Pattern Anal
#' Mach Intell, Dec 2010. This implementation builds upon this work in the referenced
#' manuscript by adding the ability to attenuate adjustments of scale, rotation, translation, and sigma^2
#' as the algorithm iterates. This helps to promote convergence in cases where certain
#' degrees of freedom are needed but do not dominate. This is done through the
#' rateSRTSigma parameter and is in addition to the w, sigma2, and tol parameters
#'
#' @param X matrix representing N points and D dimensions. These points are the targets of registration. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param Y matrix representing M points and D dimensions. These points are the those to be registered. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param w double value (0 <= w < 1) for accomodating noise and outliers
#' @param maxIter integer maximum number of iterations for the EM algorithm
#' @param tol double indicating the relative change/error in the value of the objective function, Q, at which the algorithm should stop and return a result
#' @param plot boolean indicating whether to visualize/plot the two point sets and the current state of the transformed point set after each maximization step
#' @param sigma2 double indicating the initial variance around each point for building the Gaussian Mixture Model from the point set
#' @param rateSRTSigma numeric vector with 4 values between 0 and 1 in sequence indicating how rapidly to adjust Scale, Rotation, Translation, and Sigma. Putting 0 forces no change in the parameter (e.g., Scale remains 1, Rotation remains 0, Translation remains 0, Sigma2 will stay the same as provided). The closer to 0, the slower that parameter changes.
#'
#' @return list(P, X, Y, YT, Sx, Sy, R, tr, iter, sigma2, eps10, nErr)
#' @export
#'
#' @examples demoCPD(rigid=FALSE)
cpd.affine <- function(X,Y, w, maxIter, tol, plot, sigma2, rateSRTSigma=c(1,1,1,1))
{
     # Gather some basic values
     l(N, D) %=% dim(X);
     l(M, D) %=% dim(Y);

     # Initialize variables
     B <- diag(2)
     tr <- matrix(c(0,0))
     YT <- Y

     # Iterate
     iter <- 0
     eps <- .Machine$double.eps
     nErr <- tol + 10
     nErr_old <- nErr
     sigma2_old <- sigma2
     Q <- .Machine$double.xmax;
     if(plot)
     {
          plotXYYT(X, Y, YT, includeYT=F)
     }
     while((iter < maxIter) && nErr_old > tol && (sigma2 > 10*eps))
     {
          Q_old <- Q;

          # Calculate Pmn with current information
          P <- getP(X=X, Y=YT, sigma2=sigma2, w=w, B=B, tr=tr);

          # Update the information
          Np <- sum(P);
          mux <- (1/Np) * t(X) %*% t(P) %*% get1(P, flip=F)
          muy <- (1/Np) * t(YT) %*% P %*% get1(P, flip=T)

          Xhat <- X - repmat(t(mux), m=N)
          Yhat <- YT - repmat(t(muy), m=M)

          B <- ( t(Xhat) %*% t(P) %*% Yhat ) %*% matrix.inverse( t(Yhat) %*% (diag(as.vector(P %*% get1(P, flip=T))) %*% Yhat))
          l(Sx, Sy, Theta) %=% getRS(B)
          Sx <- 1 + rateSRTSigma[1]*(Sx-1) # Attentuate rate of change in Sx
          Sy <- 1 + rateSRTSigma[1]*(Sy-1) # Attentuate rate of change in Sy
          Theta <- rateSRTSigma[2]*Theta # Attentuate rate of change in Theta
          B <- makeB(Sx, Sy, Theta)
          tr <- mux - B %*% muy
          tr <- rateSRTSigma[3]*tr # Attentuate rate of change in tr
          sigma2 <- abs((1/(Np * D))*(matrix.trace(t(Xhat) %*% (diag(as.vector(t(P) %*% get1(P, flip=F))) %*% Xhat)) - matrix.trace( t(Xhat) %*% t(P) %*% Yhat %*% t(B)))) # abs used here to avoid rounding errors that lead to negative values
          sigma2 <- (1/rateSRTSigma[4]) * sigma2 # Attentuate rate of convergence in sigma2
          if(is.nan(sigma2))
          {
               # Try to work through any issues of calculating sigma2 by just using the previous value
               sigma2 <- sigma2_old
          }

          # Calculate transformed points using new information
          YT <- transformY(Y=YT, B=B, tr=tr)

          # Finish up
          iter <- iter+1

          Q <- getQ(P=P, X=X, Y=YT, B=B, tr=tr, sigma2=sigma2)

          nErr <- abs((Q_old-Q)/Q)
          if(!is.nan(nErr))
          {
               nErr_old <- nErr
          }

          if(!is.nan(sigma2))
          {
               sigma2_old <- sigma2
          }

          cat('Sigma^2=', sigma2, "\tNormError=", nErr, "\tSx=", Sx, "\tSy", Sy, "\tR=", getRS(B)$Theta, "\ttr", tr, "\n")

          if(plot)
          {
               plotXYYT(X, Y, YT)
          }
     }
     return(list(P=P, X=X, Y=Y, YT=YT, B=B, tr=tr, iter=iter, sigma2=sigma2, eps10=10*eps, nErr=nErr))
}

#' Coherent Point Drift Registration (Affine Transformations, annealing version)
#'
#' This implementation is based upon the paper by Andriy Myronenko and Xubo Song
#' titled "Point set registration: coherent point drift" in IEEE Trans Pattern Anal
#' Mach Intell, Dec 2010. This implementation builds upon this work in the referenced
#' manuscript by adding the ability to attenuate adjustments of scale, rotation, translation, and sigma^2
#' as the algorithm iterates. This helps to promote convergence in cases where certain
#' degrees of freedom are needed but do not dominate. This is done through the
#' rateSRTSigma parameter and is in addition to the w, sigma2, and tol parameters
#'
#' In the annealing version, instead of calcualting the new sigma2 at each step, we
#' take the user provided sigma2 and slowly reduce it. This is done by multiplying
#' the sigma2 at each step by rateSRTSigma[4] (default 0.9)
#'
#' @param X matrix representing N points and D dimensions. These points are the targets of registration. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param Y matrix representing M points and D dimensions. These points are the those to be registered. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param w double value (0 <= w < 1) for accomodating noise and outliers
#' @param maxIter integer maximum number of iterations for the EM algorithm
#' @param tol double indicating the relative change/error in the value of the objective function, Q, at which the algorithm should stop and return a result
#' @param plot boolean indicating whether to visualize/plot the two point sets and the current state of the transformed point set after each maximization step
#' @param sigma2 double indicating the initial variance around each point for building the Gaussian Mixture Model from the point set
#' @param rateSRTSigma numeric vector with 4 values between 0 and 1 in sequence indicating how rapidly to adjust Scale, Rotation, Translation, and Sigma. Putting 0 forces no change in the parameter (e.g., Scale remains 1, Rotation remains 0, Translation remains 0, Sigma2 will stay the same as provided). The closer to 0, the slower that parameter changes.
#'
#' @return list(P, X, Y, YT, Sx, Sy, R, tr, iter, sigma2, eps10, nErr)
#' @export
#'
#' @examples demoCPD(rigid=FALSE)
cpd.affine.annealing <- function(X,Y, w, maxIter, tol, plot, sigma2, rateSRTSigma=c(1,1,1,0.9))
{
     # Gather some basic values
     l(N, D) %=% dim(X);
     l(M, D) %=% dim(Y);

     # Initialize variables
     B <- diag(2)
     tr <- matrix(c(0,0))
     YT <- Y

     # Iterate
     iter <- 0
     eps <- .Machine$double.eps
     nErr <- tol + 10
     nErr_old <- nErr
     sigma2_old <- sigma2
     Q <- .Machine$double.xmax;
     if(plot)
     {
          plotXYYT(X, Y, YT, includeYT=F)
     }
     while((iter < maxIter) && nErr_old > tol && (sigma2 > 10*eps))
     {
          Q_old <- Q;

          # Calculate Pmn with current information
          P <- getP(X=X, Y=YT, sigma2=sigma2, w=w, B=B, tr=tr);

          # Update the information
          Np <- sum(P);
          mux <- (1/Np) * t(X) %*% t(P) %*% get1(P, flip=F)
          muy <- (1/Np) * t(YT) %*% P %*% get1(P, flip=T)

          Xhat <- X - repmat(t(mux), m=N)
          Yhat <- YT - repmat(t(muy), m=M)

          B <- ( t(Xhat) %*% t(P) %*% Yhat ) %*% matrix.inverse( t(Yhat) %*% (diag(as.vector(P %*% get1(P, flip=T))) %*% Yhat))
          l(Sx, Sy, Theta) %=% getRS(B)
          Sx <- 1 + rateSRTSigma[1]*(Sx-1) # Attentuate rate of change in Sx
          Sy <- 1 + rateSRTSigma[1]*(Sy-1) # Attentuate rate of change in Sy
          Theta <- rateSRTSigma[2]*Theta # Attentuate rate of change in Theta
          B <- makeB(Sx, Sy, Theta)
          tr <- mux - B %*% muy
          tr <- rateSRTSigma[3]*tr # Attentuate rate of change in tr
          # sigma2 <- abs((1/(Np * D))*(matrix.trace(t(Xhat) %*% (diag(as.vector(t(P) %*% get1(P, flip=F))) %*% Xhat)) - matrix.trace( t(Xhat) %*% t(P) %*% Yhat %*% t(B)))) # abs used here to avoid rounding errors that lead to negative values
          sigma2 <- rateSRTSigma[4] * sigma2 # Attentuate rate of convergence in sigma2
          if(is.nan(sigma2))
          {
               # Try to work through any issues of calculating sigma2 by just using the previous value
               sigma2 <- sigma2_old
          }

          # Calculate transformed points using new information
          YT <- transformY(Y=YT, B=B, tr=tr)

          # Finish up
          iter <- iter+1

          Q <- getQ(P=P, X=X, Y=YT, B=B, tr=tr, sigma2=sigma2)

          nErr <- abs((Q_old-Q)/Q)
          if(!is.nan(nErr))
          {
               nErr_old <- nErr
          }

          if(!is.nan(sigma2))
          {
               sigma2_old <- sigma2
          }

          cat('Sigma^2=', sigma2, "\tNormError=", nErr, "\tSx=", Sx, "\tSy", Sy, "\tR=", getRS(B)$Theta, "\ttr", tr, "\n")

          if(plot)
          {
               plotXYYT(X, Y, YT)
          }
     }
     return(list(P=P, X=X, Y=Y, YT=YT, B=B, tr=tr, iter=iter, sigma2=sigma2, eps10=10*eps, nErr=nErr))
}

#' Coherent Point Drift Registration (Rigid Transformations)
#'
#' This implementation is based upon the paper by Andriy Myronenko and Xubo Song
#' titled "Point set registration: coherent point drift" in IEEE Trans Pattern Anal
#' Mach Intell, Dec 2010. This implementation builds upon this work in the referenced
#' manuscript by adding the ability to attenuate adjustments of scale, rotation, translation, and sigma^2
#' as the algorithm iterates. This helps to promote convergence in cases where certain
#' degrees of freedom are needed but do not dominate. This is done through the
#' rateSRTSigma parameter and is in addition to the w, sigma2, and tol parameters
#'
#' @references <https://www.ncbi.nlm.nih.gov/pubmed/20975122>
#'
#' @param X matrix representing N points and D dimensions. These points are the targets of registration. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param Y matrix representing M points and D dimensions. These points are the those to be registered. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param w double value (0 <= w < 1) for accomodating noise and outliers
#' @param maxIter integer maximum number of iterations for the EM algorithm
#' @param tol double indicating the relative change/error in the value of the objective function, Q, at which the algorithm should stop and return a result
#' @param plot boolean indicating whether to visualize/plot the two point sets and the current state of the transformed point set after each maximization step
#' @param sigma2 double indicating the initial variance around each point for building the Gaussian Mixture Model from the point set
#' @param rateSRTSigma numeric vector with 4 values between 0 and 1 in sequence indicating how rapidly to adjust Scale, Rotation, Translation, and Sigma. Putting 0 forces no change in the parameter (e.g., Scale remains 1, Rotation remains 0, Translation remains 0, Sigma2 will stay the same as provided). The closer to 0, the slower that parameter changes.
#'
#' @return list(P, X, Y, YT, Sx, Sy, R, tr, iter, sigma2, eps10, nErr)
#' @export
#'
#' @examples demoCPD(rigid=TRUE)
cpd.rigid <- function(X,Y, w, maxIter, tol, plot, sigma2, rateSRTSigma=c(1,1,1,1))
{
     # Gather some basic values
     l(N, D) %=% dim(X);
     l(M, D) %=% dim(Y);

     # Initialize variables
     Sxy <- 1
     R <- diag(2)
     tr <- matrix(c(0,0))
     YT <- Y;

     # Iterate
     iter <- 0;
     eps <- .Machine$double.eps
     nErr <- tol + 10
     nErr_old <- nErr
     sigma2_old <- sigma2
     Q <- .Machine$double.xmax
     if(plot)
     {
          plotXYYT(X, Y, YT, includeYT=T)
     }
     while((iter < maxIter) && nErr_old > tol && (sigma2 > 10*eps))
     {
          Q_old <- Q;

          P <- getP(X=X, Y=YT, sigma2=sigma2, w=w, B=Sxy*R, tr=tr);

          Np <- sum(P);
          mux <- (1/Np) * t(X) %*% t(P) %*% get1(P, flip=F)
          muy <- (1/Np) * t(YT) %*% P %*% get1(P, flip=T)

          Xhat <- X - repmat(t(mux), m=N)
          Yhat <- YT - repmat(t(muy), m=M)

          A <- t(Xhat) %*% t(P) %*% Yhat
          l(d, u, v) %=% svd(A)
          littlec <- det(u %*% t(u))
          C <- diag(c(rep(1, times=D-1), littlec))
          R <- u %*% C %*% t(v)
          R <- makeR(rateSRTSigma[2]*getRS(R)$Theta) # attenuate the rate of change of R
          Sxy <- matrix.trace(t(A) %*% R) / matrix.trace(t(Yhat) %*% diag(as.vector(P %*% get1(P, flip=T))) %*% Yhat)
          Sxy <- 1 + rateSRTSigma[1]*(Sxy-1) # attenuate the rate of change of R
          tr <- (mux - (Sxy * R) %*% muy)
          tr <- rateSRTSigma[3]*tr  # attenuate the rate of change of tr
          temp1 <- matrix.trace(t(Xhat) %*% (diag(as.vector(t(P) %*% get1(P, flip=F))) %*% Xhat))
          sigma2 <- abs((1/(Np * D))*(temp1 - Sxy * matrix.trace(t(A) %*% R))) # abs used here to avoid rounding errors that lead to negative values
          sigma2 <- (1/rateSRTSigma[4]) * sigma2 # attenuate the rate of convergence of sigma2

          if(is.nan(sigma2))
          {
               # Try to work through any issues of calculating sigma2 by just using the previous value
               sigma2 <- sigma2_old
          }

          YT <- transformY(Y=YT, B=Sxy*R, tr=tr)

          # Finish up
          iter <- iter+1

          Q <- getQ(P=P, X=X, Y=YT, B=Sxy*R, tr=tr, sigma2=sigma2)

          nErr <- abs((Q_old-Q)/Q)
          if(!is.nan(nErr))
          {
               nErr_old <- nErr
          }

          if(!is.nan(sigma2))
          {
               sigma2_old <- sigma2
          }

          cat('Sigma^2=', sigma2, "\tNormError=", nErr, "\ts=", Sxy, "\tR=", getRS(R)$Theta, "\ttr", tr, "\n")

          if(plot)
          {
               plotXYYT(X, Y, YT)
          }
     }
     return(list(P=P, X=X, Y=Y, YT=YT, Sx=Sxy, Sy=Sxy, R=R, tr=tr, iter=iter, sigma2=sigma2, eps10=10*eps, nErr=nErr))
}

#' Coherent Point Drift Registration (Rigid Transformations, annealing version)
#'
#' This implementation is based upon the paper by Andriy Myronenko and Xubo Song
#' titled "Point set registration: coherent point drift" in IEEE Trans Pattern Anal
#' Mach Intell, Dec 2010. This implementation builds upon this work in the referenced
#' manuscript by adding the ability to attenuate adjustments of scale, rotation, translation, and sigma^2
#' as the algorithm iterates. This helps to promote convergence in cases where certain
#' degrees of freedom are needed but do not dominate. This is done through the
#' rateSRTSigma parameter and is in addition to the w, sigma2, and tol parameters
#'
#' #' In the annealing version, instead of calcualting the new sigma2 at each step, we
#' take the user provided sigma2 and slowly reduce it. This is done by multiplying
#' the sigma2 at each step by rateSRTSigma[4] (default 0.9)
#'
#' @references <https://www.ncbi.nlm.nih.gov/pubmed/20975122>
#'
#' @param X matrix representing N points and D dimensions. These points are the targets of registration. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param Y matrix representing M points and D dimensions. These points are the those to be registered. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param w double value (0 <= w < 1) for accomodating noise and outliers
#' @param maxIter integer maximum number of iterations for the EM algorithm
#' @param tol double indicating the relative change/error in the value of the objective function, Q, at which the algorithm should stop and return a result
#' @param plot boolean indicating whether to visualize/plot the two point sets and the current state of the transformed point set after each maximization step
#' @param sigma2 double indicating the initial variance around each point for building the Gaussian Mixture Model from the point set
#' @param rateSRTSigma numeric vector with 4 values between 0 and 1 in sequence indicating how rapidly to adjust Scale, Rotation, Translation, and Sigma. Putting 0 forces no change in the parameter (e.g., Scale remains 1, Rotation remains 0, Translation remains 0, Sigma2 will stay the same as provided). The closer to 0, the slower that parameter changes.
#'
#' @return list(P, X, Y, YT, Sx, Sy, R, tr, iter, sigma2, eps10, nErr)
#' @export
#'
#' @examples demoCPD(rigid=TRUE)
cpd.rigid.annealing <- function(X,Y, w, maxIter, tol, plot, sigma2, rateSRTSigma=c(1,1,1,0.9))
{
     # Gather some basic values
     l(N, D) %=% dim(X);
     l(M, D) %=% dim(Y);

     # Initialize variables
     Sxy <- 1
     R <- diag(2)
     tr <- matrix(c(0,0))
     YT <- Y;

     # Iterate
     iter <- 0;
     eps <- .Machine$double.eps
     nErr <- tol + 10
     nErr_old <- nErr
     sigma2_old <- sigma2
     Q <- .Machine$double.xmax
     if(plot)
     {
          plotXYYT(X, Y, YT, includeYT=T)
     }
     while((iter < maxIter) && nErr_old > tol && (sigma2 > 10*eps))
     {
          Q_old <- Q;

          P <- getP(X=X, Y=YT, sigma2=sigma2, w=w, B=Sxy*R, tr=tr);

          Np <- sum(P);
          mux <- (1/Np) * t(X) %*% t(P) %*% get1(P, flip=F)
          muy <- (1/Np) * t(YT) %*% P %*% get1(P, flip=T)

          Xhat <- X - repmat(t(mux), m=N)
          Yhat <- YT - repmat(t(muy), m=M)

          A <- t(Xhat) %*% t(P) %*% Yhat
          l(d, u, v) %=% svd(A)
          littlec <- det(u %*% t(u))
          C <- diag(c(rep(1, times=D-1), littlec))
          R <- u %*% C %*% t(v)
          R <- makeR(rateSRTSigma[2]*getRS(R)$Theta) # attenuate the rate of change of R
          Sxy <- matrix.trace(t(A) %*% R) / matrix.trace(t(Yhat) %*% diag(as.vector(P %*% get1(P, flip=T))) %*% Yhat)
          Sxy <- 1 + rateSRTSigma[1]*(Sxy-1) # attenuate the rate of change of R
          tr <- (mux - (Sxy * R) %*% muy)
          tr <- rateSRTSigma[3]*tr  # attenuate the rate of change of tr
          temp1 <- matrix.trace(t(Xhat) %*% (diag(as.vector(t(P) %*% get1(P, flip=F))) %*% Xhat))
          # sigma2 <- abs((1/(Np * D))*(temp1 - Sxy * matrix.trace(t(A) %*% R))) # abs used here to avoid rounding errors that lead to negative values
          sigma2 <- rateSRTSigma[4] * sigma2 # attenuate the rate of convergence of sigma2

          if(is.nan(sigma2))
          {
               # Try to work through any issues of calculating sigma2 by just using the previous value
               sigma2 <- sigma2_old
          }

          YT <- transformY(Y=YT, B=Sxy*R, tr=tr)

          # Finish up
          iter <- iter+1

          Q <- getQ(P=P, X=X, Y=YT, B=Sxy*R, tr=tr, sigma2=sigma2)

          nErr <- abs((Q_old-Q)/Q)
          if(!is.nan(nErr))
          {
               nErr_old <- nErr
          }

          if(!is.nan(sigma2))
          {
               sigma2_old <- sigma2
          }

          cat('Sigma^2=', sigma2, "\tNormError=", nErr, "\ts=", Sxy, "\tR=", getRS(R)$Theta, "\ttr", tr, "\n")

          if(plot)
          {
               plotXYYT(X, Y, YT)
          }
     }
     return(list(P=P, X=X, Y=Y, YT=YT, Sx=Sxy, Sy=Sxy, R=R, tr=tr, iter=iter, sigma2=sigma2, eps10=10*eps, nErr=nErr))
}


#' demoCPD
#'
#' Run a small demo that creates to pointsets that differ by
#' small random noise and an equal offset in both the x and
#' y direction.
#'
#' @param rigid boolean indicating whether to demo affine registration or rigid registration
#'
#' @export
demoCPD <- function(rigid=F, annealing=F)
{
     set.seed(12345)
     n <- 100
     alpha <- 10 # scaling parameter
     xpoints <- matrix(alpha*runif(n*2), ncol=2)
     xpoints[1:n,2] <- xpoints[1:n,2] + alpha*2
     x2points <- matrix(alpha*runif(4), ncol=2)
     x2points[1:2,2] <- x2points[1:2,2] + alpha*2
     x2points <- rbind(xpoints, x2points)
     jitter1 <- matrix(alpha*rnorm((n+0)*2, -0.01, 0.01), ncol=2)
     jitter2 <- matrix(alpha*rnorm((n+2)*2, -0.01, 0.01), ncol=2)
     ypoints <- x2points + jitter2
     plot(xpoints, xlim=c(0, 4), ylim = c(0,4))
     points(ypoints, col='red')
     if(!rigid & !annealing)
     {
          # list(P, X, Y, YT, B, tr, iter, sigma2, eps10, nErr)
          temp <- cpd.affine(xpoints, ypoints+alpha*0.3, w=0.3, maxIter = 50, tol=1e-15, plot=T, sigma2=0.01, rateSRTSigma = c(0,0.1,1,0.0051))
     }
     else if(rigid & !annealing)
     {
          # list(P, X, Y, YT, B, tr, iter, sigma2, eps10, nErr)
          temp <- cpd.rigid(xpoints, ypoints+alpha*0.3, w=0.3, maxIter = 50, tol=1e-15, plot=T, sigma2=0.01, rateSRTSigma=c(0,0.1,1,0.0051))
     }
     else if(!rigid & annealing)
     {
          # list(P, X, Y, YT, B, tr, iter, sigma2, eps10, nErr)
          temp <- cpd.affine.annealing(xpoints, ypoints+alpha*0.3, w=0.3, maxIter = 50, tol=5e-3, plot=T, sigma2=100, rateSRTSigma=c(0,0.1,1,0.8))
     }
     else
     {
          # list(P, X, Y, YT, B, tr, iter, sigma2, eps10, nErr)
          temp <- cpd.rigid.annealing(xpoints, ypoints+alpha*0.3, w=0.3, maxIter = 50, tol=5e-3, plot=T, sigma2=100, rateSRTSigma=c(0,0.1,1,0.8))
     }

     plotXYYT(temp$X, temp$Y, temp$YT)
     print(paste0('Iterations: ', temp$iter))
     print(paste0('Error: ', temp$nErr))
     return(temp)
}

#' getCorrespondence(P, X, YT, threshDist)
#'
#' Calculate which cells correspond with which other cells
#'
#' @param P The M x N matrix of column normalized probabilities
#' @param X matrix representing N points and D dimensions. These points were the targets of registration. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param YT matrix representing M points and D dimensions. These points were the those that were registered/transformed. Each column represents a dimension (e.g., x, y, z...) and each row a point
#' @param threshDist double for calculating whether the corresponding points are within (inclusive) a certain distance.

#' @export
getCorrespondence <- function(P, X, YT, threshDist=NULL)
{
     # Get the maximum probability match for each X (i.e., get the col max's indicies)
     # We need to do it this way because the probabilities are normalized per column (vs row)
     bestYforX <- as.numeric(lapply(as.list(as.data.frame(P)), which.max))
     bestYforX <- data.frame(ym=bestYforX, xn=1:ncol(P))
     bestYforX$dist <- sqrt(rowSums((X[bestYforX$xn,]-YT[bestYforX$ym,])^2))
     # bestYforX$logProb <- log(P[as.matrix(bestYforX[,c('ym', 'xn')])])

     if(is.null(threshDist))
     {
          bestYforX$passDistThresh <- NA
     }
     else
     {
          bestYforX$passDist <- bestYforX$dist <= threshDist
     }
     return(bestYforX)
}


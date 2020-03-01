/******************************************************************************
Author, date: Rick Wicklin, NOV2012
Macro name:   NearestCorr
Description:  implement the Higham projection method to obtain the nearest 
              correlation matrix to a given symmetric correlation/covariance matrix 
Reference:    http://blogs.sas.com/content/iml/2012/11/28/
              computing-the-nearest-correlation-matrix.html
******************************************************************************/

%macro NearestCorr();
/* Project symmetric X onto S={positive semidefinite matrices}.
   Replace any negative eigenvalues of X with zero */
start ProjS(X);
   call eigen(D, Q, X); /* note X = Q*D*Q` */
   V = choose(D>.0001, D, .0001); /*pmbrown: this line changed as per comment 
                                   at bottom of sas blog*/
   W = Q#sqrt(V`);      /* form Q*diag(V)*Q` */
   return( W*W` );      /* W*W` = Q*diag(V)*Q` */
finish;
 
/* project square X onto U={matrices with unit diagonal}.
   Return X with the diagonal elements replaced by ones. */
start ProjU(X);
   n = nrow(X);
   Y = X;
   diagIdx = do(1, n*n, n+1);
   Y[diagIdx] = 1;      /* set diagonal elements to 1 */
   return ( Y );
finish;
 
/* Helper function: the infinity norm is the max abs value of the row sums */
start MatInfNorm(A);
   return( max(abs(A[,+])) );
finish;
 
/* Given a symmetric correlation matrix, A, 
   project A onto the space of positive semidefinite matrices.
   The function uses the algorithm of Higham (2002) to return
   the matrix X that is closest to A in the Frobenius norm. */
start NearestCorr(A);
   maxIter = 100; tol = 1e-8;  /* parameters...you can change these */
   iter = 1;      maxd   = 1;  /* initial values */ 
   Yold = A;  Xold = A;  dS = 0;
 
   do while( (iter <= maxIter) & (maxd > tol) );
     R = Yold - dS; /* dS is Dykstra's correction */
     X = ProjS(R);  /* project onto S={positive semidefinite} */
     dS = X - R;
     Y = ProjU(X);  /* project onto U={matrices with unit diagonal} */
 
     /* How much has X changed? (Eqn 4.1) */
     dx = MatInfNorm(X-Xold) / MatInfNorm(X);
     dy = MatInfNorm(Y-Yold) / MatInfNorm(Y);
     dxy = MatInfNorm(Y - X) / MatInfNorm(Y);
     maxd = max(dx,dy,dxy);
 
     iter = iter + 1;
     Xold = X; Yold = Y; /* update matrices */
   end;
   return( X ); /* X is positive semidefinite */
finish;
%mend NearestCorr;

***end********************************************************;

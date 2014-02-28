sgmViaIP <- function (A, B,m){
  require(Matrix)
  
  totv = nrow(A);
  n= totv-m;
  A11=A[1:m,1:m];
  A12=A[1:m,m+(1:n)];
  A21=A[m+(1:n),1:m];
  A22=A[m+(1:n),m+(1:n)];
  B11=B[1:m,1:m];
  B12=B[1:m,m+(1:n)];
  B21=B[m+(1:n),1:m];
  B22=B[m+(1:n),m+(1:n)];
  
  vecB12=B12
  dim(vecB12)<- c(m*n,1);

  vecA21=A21 ;
  dim(vecA21)<- c(m*n,1);
  eye.n<-diag(n)
  eye.nsq.2mn = diag(n^2+2*m*n)
  vec.1 <- matrix(1,1,n)
  
  M11 <- (eye.n%x% A22)
  M11 <-  M11 - (t(B22) %x% eye.n)
  M11=rBind( M11, 
      eye.n %x% A12, t(B21) %x% eye.n);
  M21= rBind( eye.n %x% vec.1 , vec.1 %x% eye.n);
  M=rBind(cBind( M11, eye.nsq.2mn , -eye.nsq.2mn),
          cBind(M21,  matrix(0,2*n,2*n^2+4*m*n))
          ); 
  f=cBind( matrix(0,1,n^2) , matrix(1,1, 2*n^2+4*m*n) );
  Mineq= rep(0,3*n^2+4*m*n);
  M=Matrix(M,sparse=TRUE);
  b=rBind(matrix(0,n^2,1) , vecB12 , vecA21 , matrix(1,2*n,1));
  
  sense_eq = rep('=',n^2+2*m*n+2*n);
 # sense_ineq = rep('<',0, 1];
  #sense = [sense_eq; sense_ineq];
   sense = sense_eq
  #Binary and Real Variables
  vtype1 = rep('B',n^2);
  vtype2 = rep('B',2*n^2+4*m*n);
  vtype = c(vtype1, vtype2);
  
  model = list()
  model$A = M;
  model$obj = f;
  model$modelsense = 'min';
  model$rhs = b;
  model$sense = sense;
  model$vtype = vtype;
 require(gurobi)
  result = gurobi(model);
  print(result$status)
  x = result$x;
  x = round(x);
  
  P=x;
  dim(P)<- c(n,n)
  P.full <- diag(totv)
  P.full[(m+1):totv,(m+1):totv] <- P
  temp=P%*%matrix(1:n, n, 1);
   alignment=c( 1:m , t(temp+m) );
  
   return (list(matching=alignment,matching.seeded=t(temp+m), 
                perm.mat=P.full,  perm.mat.unseeded = P))
}
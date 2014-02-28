sgmViaIP <- function (A, B,m){
  require(Matrix)
  
   n = nrow(A);
  test.v.count= n-m;
  A11=A[1:m,1:m];
  A12=A[1:m,m+(1:test.v.count)];
  A21=A[m+(1:test.v.count),1:m];
  A22=A[m+(1:test.v.count),m+(1:test.v.count)];
  B11=B[1:m,1:m];
  B12=B[1:m,m+(1:test.v.count)];
  B21=B[m+(1:test.v.count),1:m];
  B22=B[m+(1:test.v.count),m+(1:test.v.count)];
  
  vecB12=B12
  dim(vecB12)<- c(m*test.v.count,1);

  vecA21=A21 ;
  dim(vecA21)<- c(m*test.v.count,1);
  eye.test.v.count<-diag(test.v.count)
  eye.test.v.countsq.2mtest.v.count = diag(test.v.count^2+2*m*test.v.count)
  vec.1 <- matrix(1,1,test.v.count) 
  M11=rbind( eye.test.v.count%x% A22- t(B22) %x% eye.test.v.count, 
      eye.test.v.count %x% A12, t(B21) %x% eye.test.v.count);
  M21= rbind( eye.test.v.count %x% vec.1 , vec.1 %x% eye.test.v.count);
  M=rbind(cbind( M11, eye.test.v.countsq.2mtest.v.count , -eye.test.v.countsq.2mtest.v.count),
          cbind(M21,  matrix(0,2*test.v.count,2*test.v.count^2+4*m*test.v.count))
          ); 
  f=cbind( matrix(0,1,test.v.count^2) , matrix(1,1, 2*test.v.count^2+4*m*test.v.count) );
  Mineq= rep(0,3*test.v.count^2+4*m*test.v.count);
  M=Matrix(M,sparse=TRUE);
  b=rbind(matrix(0,test.v.count^2,1) , vecB12 , vecA21 , matrix(1,2*test.v.count,1));
  
  sense_eq = rep('=',test.v.count^2+2*m*test.v.count+2*test.v.count);
 # setest.v.countse_itest.v.counteq = rep('<',0, 1];
  #setest.v.countse = [sense_eq; setest.v.countse_itest.v.counteq];
   sense = sense_eq
  #Bitest.v.countary atest.v.countd Real Variables
  vtype1 = rep('B',test.v.count^2);
  vtype2 = rep('B',2*test.v.count^2+4*m*test.v.count);
  vtype = c(vtype1, vtype2);
  
  model = list()
  model$A = M;
  model$obj = f;
  model$modelsense = 'min';
  model$rhs = b;
  model$sense = sense;
  model$vtype = vtype;
  result = gurobi[model];
  print(result$status)
  x = result$x;
  x = round(x);
  
  P=x;
  dim(P)<- c(test.v.count,test.v.count)
  temp=P%*%matrix(1:test.v.count, test.v.count, 1);
   alignment=c( 1:m , t(temp+m) );
   return (list(matching=alignment, perm.mat = P))
}
function [alignment]=GM_exact(A,B,m)

% [alignment]=GM(A,B,m)
% A and B are adjacency matrices with first m vertices as seeds
% alignment is the graph matching; 
% eg if alignment is [1 2 3 25 68 ... ]
% then vertex1ofA->vertex1ofB ... vertex4ofA->vertex25ofB... 
% DEF Feb 18, 2014, ready

addpath('..\..\..\..\gurobi562\win64\matlab');

[mplusn,~]=size(A);
n=mplusn-m;
A11=A(1:m,1:m);
A12=A(1:m,m+1:m+n);
A21=A(m+1:m+n,1:m);
A22=A(m+1:m+n,m+1:m+n);
B11=B(1:m,1:m);
B12=B(1:m,m+1:m+n);
B21=B(m+1:m+n,1:m);
B22=B(m+1:m+n,m+1:m+n);

vecB12=zeros(m*n,1);
for i=1:n
    vecB12( (i-1)*m+1:i*m ) = B12(:,i);
end
vecA21=zeros(m*n,1);
for i=1:m
    vecA21( (i-1)*n+1:i*n ) = A21(:,i);
end

M11=[ kron(eye(n),A22)-kron(B22',eye(n)); 
      kron(eye(n),A12); kron(B21',eye(n))];
M21=[ kron(eye(n),ones(1,n)) ; kron(ones(1,n),eye(n))];
M=[ M11  eye(n^2+2*m*n) -eye(n^2+2*m*n) ;
    M21  zeros(2*n,2*n^2+4*m*n)]; 
f=[ zeros(1,n^2) ones(1,2*n^2+4*m*n) ];
Mineq= zeros(0,3*n^2+4*m*n);
M=sparse(M);
b=[zeros(n^2,1) ; vecB12 ; vecA21 ; ones(2*n,1)];

sense_eq = repmat('=',n^2+2*m*n+2*n, 1);
sense_ineq = repmat('<',0, 1);
sense = [sense_eq; sense_ineq];
%Binary and Real Variables
vtype1 = repmat('B',n^2, 1);
vtype2 = repmat('C',2*n^2+4*m*n,1);
vtype = [vtype1; vtype2];

model.A = M;
model.obj = f;
model.modelsense = 'min';
model.rhs = b;
model.sense = sense;
model.vtype = vtype;
result = gurobi(model);
result.status
x = result.x;
x=round(x);

P=zeros(n,n);
for i=1:n
    P(:,i)=x( (i-1)*n+1:i*n );
end
temp=P*[1:n]';
alignment=[ [1:m] temp'+m ];

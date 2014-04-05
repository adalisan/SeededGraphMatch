function [alignment,fwDISAGREE,exactDISAGREE,fwruntime,exactruntime]=GM_exact_with_FW_init(A,B,m)

%[alignment,fwDISAGREE,exactDISAGREE,fwruntime,exactruntime]=GMexactafterFW(A,B,m)
% A and B are adjacency matrices with first m vertices as seeds
% alignment is the graph matching; 
% eg if alignment is [1 2 3 25 68 ... ]
% then vertex1ofA->vertex1ofB ... vertex4ofA->vertex25ofB... 
% 
% This version runs Frank-Wolfe first, then uses to seed the exact IP solver
% alpha version Feb 20, 2014

addpath('F:\gurobi562\win64\matlab');

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


% begin Frank-Wolfe

tic

patience=20;
tol=.99;
P=ones(n,n)/n;
toggle=1;
iter=0;
while (toggle==1)&(iter<patience)
    iter=iter+1;
    Grad=A22*P*B22'+A22'*P*B22+A21*B21'+A12'*B12;
    ind=YiCaoHungarian(-Grad);
    T=eye(n);
    T=T(ind,:);
    c=trace(A22'*P*B22*P');
    d=trace(A22'*T*B22*P')+trace(A22'*P*B22*T');
    e=trace(A22'*T*B22*T');
    u=trace(P'*A21*B21'+P'*A12'*B12);
    v=trace(T'*A21*B21'+T'*A12'*B12);
    alpha=-(d-2*e+u-v)/(2*(c-d+e));
    f0=0;
    f1=c-e+u-v;
    falpha=(c-d+e)*alpha^2+(d-2*e+u-v)*alpha;
    if (alpha<tol)&(alpha>0)&(falpha>f0)&(falpha>f1)
        P=alpha*P+(1-alpha)*T;
    elseif (f0>f1)
        P=T;
    else
        toggle=0;
    end
end
ind=YiCaoHungarian(-P);
fwP=eye(n);
fwP=fwP(ind,:);
bigP=[eye(m) zeros(m,n) ; zeros(n,m) fwP]; 
fwDISAGREE=(1/2)*sum(sum(abs(A-bigP*B*bigP')));
vecfwP=zeros(n^2,1);
for i=1:n
    vecfwP( (i-1)*n+1:i*n ) = fwP(:,i);
end

fwruntime=toc


% begin exact

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
M=sparse(M);
b=[zeros(n^2,1) ; vecB12 ; vecA21 ; ones(2*n,1)];
shortb=[zeros(n^2,1) ; vecB12 ; vecA21];

st1= subplus(-(M11*vecfwP-shortb));
st2= subplus(M11*vecfwP-shortb);
fwst=[vecfwP; st1 ; st2 ];

sense_eq = repmat('=',n^2+2*m*n+2*n, 1);
sense_ineq = repmat('<',0, 1);
sense = [sense_eq; sense_ineq];
%Binary and Real Variables
vtype1 = repmat('B',n^2, 1);
vtype2 = repmat('C',2*n^2+4*m*n,1);
vtype = [vtype1; vtype2];

tic

model.A = M;
model.obj = f;
model.modelsense = 'min';
model.rhs = b;
model.sense = sense;
model.vtype = vtype;
model.start = fwst;
result = gurobi(model);
result.status
x = result.x;
x=round(x);

exactruntime=toc

P=zeros(n,n);
for i=1:n
    P(:,i)=x( (i-1)*n+1:i*n );
end

bigP=[eye(m) zeros(m,n) ; zeros(n,m) P]; 
exactDISAGREE=(1/2)*sum(sum(abs(A-bigP*B*bigP')));

temp=P*[1:n]';
alignment=[ [1:m] temp'+m ];

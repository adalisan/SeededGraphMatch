function [corr,fval,fval_after_proj,P,P_proj] = seedgraphmatchell1( A,B,m )
% [corr,P] = seedgraphmatchell1( A,B,m ) is the syntax.
%  A,B are (m+n)x(m+n) adjacency matrices, 
% loops/multiedges/directededges allowed.
% It is assumed that the first m vertices of A's graph
% correspond respectively to the first m vertices of B's graph,
% P is the doubly stochastic matrix that minimizes AP-PB in ell-1 sense.
% corr gives the vertex correspondences after projecting P to 
% a permutation matrix. For example, corr=[ 1 2 7 16 30 ...
% means that the vtx1ofA-->vtx1ofB, 2-->2, 3-->7, 4-->16, 5-->30 
%  example: EXECUTE the following:
% >> v=[ [1:5] 5+randperm(35)]; B=round(rand(40,40));A=B(v,v);
% >> [corr,P] = seedgraphmatchell1( A,B,5 ) ; [v corr]
% ready June 10, 2012

[totalmn,~]=size(A);
n=totalmn-m;

A12=A(1:m,m+1:m+n);
B12=B(1:m,m+1:m+n);
A21=A(m+1:m+n,1:m);
B21=B(m+1:m+n,1:m);
A22=A(m+1:m+n,m+1:m+n);
B22=B(m+1:m+n,m+1:m+n);

c=[zeros(n^2,1) ; ones(2*n^2+4*m*n,1)];

lb=zeros(3*n^2+4*m*n,1);

vecB12=[];
vecA21=[];
for i=1:n
    vecB12=[vecB12 ; B12(:,i)];
end
for i=1:m
    vecA21=[vecA21 ; A21(:,i)];
end
beq=[zeros(n^2,1) ; vecB12 ; vecA21 ; ones(2*n,1)];

F=[ kron(eye(n),A22) - kron(B22',eye(n)) ;
    kron(eye(n),A12);
    kron(B21',eye(n));
    kron(eye(n),ones(1,n));
    kron(ones(1,n),eye(n))];
Sk=[ eye(n^2+2*m*n) -eye(n^2+2*m*n) ; zeros(2*n,2*n^2+4*m*n)];
Aeq=[F Sk];
opts=optimset('MaxIter',100,'TolFun',1E-9);
[x ,~, exitflag]=linprog(c,[],[],Aeq,beq,lb,[],[],opts);
exitflag
P=[];
for i=1:n
    P=[P x((i-1)*n+1:i*n)];
end
P=[ eye(m) zeros(m,n) ; zeros(n,m) P];
corr=lapjv(-P,0.01);


%assign_of_test_v=lapjv(-P,0.01);
%corr= [1:m m+assign_of_test_v];
fval= sum(sum(abs( A*P-P*B)));


 P_proj=eye(m+n);
 P_proj=P_proj(corr,:);

fval_after_proj= sum(sum(abs(A*P_proj-P_proj*B)));

end


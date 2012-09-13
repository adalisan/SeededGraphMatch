function [ corr,iter ] = seedgraphmatchell2( A,B,m )

% [corr,iter] = seedgraphmatchell2( A,B,m ) is the syntax.
%  A,B are (m+n)x(m+n) adjacency matrices, 
% loops/multiedges/directededges allowed.
% It is assumed that the first m vertices of A's graph
% correspond respectively to the first m vertices of B's graph,
% corr gives the vertex correspondences  
% For example, corr=[ 1 2 7 16 30 ...
% means that the vtx1ofA-->vtx1ofB, 2-->2, 3-->7, 4-->16, 5-->30 
%  example: EXECUTE the following:
% >> v=[ [1:5] 5+randperm(400)]; B=round(rand(405,405));A=B(v,v);
% >> [corr,P] = seedgraphmatchell2( A,B,5 ) ; [v; corr]
% ready June 1, 2012   (Donniell's code)
% (Extends Vogelstein, Conroy et al method for nonseed graphmatch to seed)
% Original Code by Donniell Fishkind,
% Sancar Adali: replaced "YiCaoHungarian"  with faster "lapjv" 
% Same or similar with ConVogHard*

[totv,~]=size(A);
n=totv-m;

A12=A(1:m,m+1:m+n);
A21=A(m+1:m+n,1:m);
A22=A(m+1:m+n,m+1:m+n);
B12=B(1:m,m+1:m+n);
B21=B(m+1:m+n,1:m);
B22=B(m+1:m+n,m+1:m+n);


patience=20;
tol=.99;
P=ones(n,n)/n;
toggle=1;
iter=0;
while (toggle==1)&(iter<patience)
    iter=iter+1;
    Grad=A22*P*B22'+A22'*P*B22+A21*B21'+A12'*B12;
    ind=lapjv(-Grad,0.01);
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
corr=lapjv(-P,0.01);
corr=[ 1:m,  m+corr];


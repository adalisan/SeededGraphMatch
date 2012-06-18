function [ corr,iter ] = ConVogHard_rQAP2( A,B,m )

% [corr,iter] = ConVogHard( A,B,m ) is the syntax.
%  A,B are (m+n)x(m+n) adjacency matrices, 
% loops/multiedges/directededges allowed.
% It is assumed that the first m vertices of A's graph
% correspond respectively to the first m vertices of B's graph,
% corr gives the vertex correspondences  
% For example, corr=[ 1 2 7 16 30 ...
% means that the vtx1ofA-->vtx1ofB, 2-->2, 3-->7, 4-->16, 5-->30 
%  example: EXECUTE the following:
% >> v=[ [1:5] 5+randperm(400)]; B=round(rand(405,405));A=B(v,v);
% >> [corr,P] = ConVogHard( A,B,5 ) ; [v; corr]
% ready June 1, 2012   (Donniell's code)


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
while (toggle==1)&&(iter<patience)
    iter=iter+1;
    Grad=-2*A21*B21'-2*A12'*B12+2*P*(B21*B21')+2*(A12'*A12)*P+2*((A22'*A22*P)+ ...
        P*(B22*B22') - A22'*P*B22-A22*P*B22');
    ind=YiCaoHungarian(Grad);
    T=eye(n);
    T=T(ind,:);
    tempmat1 = (B21*B21'+B22*B22');
    tempmat2 = (A12'*A12+A22'*A22);
    c=trace(P'*P*tempmat1)+trace(tempmat2*(P*P'))-  ...
    trace(P'*A22'*P*B22)-trace(P'*A22*P*B22');
    d=trace((T'*P+P'*T)*tempmat1)+trace(tempmat2*(T*P'+P*T'))- ...
        trace(P'*(A22'*T*B22+(A22*T*B22')))-trace(T'*(A22'*P*B22+(A22*P*B22')));
    e=trace(T'*T*tempmat1)+trace(tempmat2*(T*T'))- ...
       trace(T'*A22'*T*B22)-trace(T'*A22*T*B22');
    u=-2*trace(P'*A21*B21'+P'*A12'*B12);
    v=-2*trace(T'*A21*B21'+T'*A12'*B12);
    c = -c;
    d = -d;
    e = -e;
    u = -u;
    v = -v;
    alpha=-(d-2*e+u-v)/(2*(c-d+e));
    f0=0;
    f1=c-e+u-v;
    falpha=(c-d+e)*alpha^2+(d-2*e+u-v)*alpha;
    if (alpha<tol)&&(alpha>0)&&(falpha>f0)&&(falpha>f1)
        P=alpha*P+(1-alpha)*T;
    elseif (f0>f1)
        P=T;
    else
        toggle=0;
    end
end
corr=YiCaoHungarian(-P);
corr=[ 1:m,  m+corr];


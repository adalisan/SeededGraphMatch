function [corr,P] = graphmatchHARDSEED( A,B,r )
% [corr,P] = graphmatchHARDSEED( A,B,r ) is the syntax.
%  A,B are (r+s)x(r+s) adjacency matrices, 
% loops/multiedges/directededges allowed.
% It is assumed that the first r vertices of A's graph
% correspond respectively to the first r vertices of B's graph,
% P is the doubly stochastic matrix that minimizes AP-PB in ell-1 sense.
% corr gives the vertex correspondences after projecting P to 
% a permutation matrix. For example, corr=[ 1 2 7 16 30 ...
% means that the vtx1ofA-->vtx1ofB, 2-->2, 3-->7, 4-->16, 5-->30 
%  example: EXECUTE the following:
% >> v=[ [1:5] 5+randperm(35)]'; B=round(rand(40,40));A=B(v,v);
% >> [corr,P] = graphmatchHARDSEED( A,B,5 ) ; [v corr]
% ready to go, Jan 8, 2012.

[n,~]=size(A);
s=n-r;

A12=A(1:r,r+1:r+s);
B12=B(1:r,r+1:r+s);
A21=A(r+1:r+s,1:r);
B21=B(r+1:r+s,1:r);
A22=A(r+1:r+s,r+1:r+s);
B22=B(r+1:r+s,r+1:r+s);

c=[zeros(s^2,1) ; ones(2*s^2+4*r*s,1)];

lb=zeros(3*s^2+4*r*s,1);

vecB12=[];
vecB21T=[];
for i=1:s
    vecB12=[vecB12 ; B12(:,i)];
    vecB21T=[vecB21T ; B21(i,:)'];
end
beq=[zeros(s^2,1) ; vecB12 ; vecB21T ; ones(2*s,1)];

F=[ kron(eye(s),A22) - kron(B22',eye(s)) ;
    kron(eye(s),A12);
    kron(eye(s),A21');
    kron(eye(s),ones(1,s));
    kron(ones(1,s),eye(s))];
Sk=[ eye(s^2+2*r*s) -eye(s^2+2*r*s) ; zeros(2*s,2*s^2+4*r*s)];
Aeq=[F Sk];

x=linprog(c,[],[],Aeq,beq,lb,[]);
P=[];
for i=1:s
    P=[P x((i-1)*s+1:i*s)];
end
P=[ eye(r) zeros(r,s) ; zeros(s,r) P];
corr=YiCaoHungarian(-P)';

end


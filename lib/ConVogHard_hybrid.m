function [ corr,iter ,fvals] = ConVogHard_hybrid( A,B,m )

% [corr,iter] = ConVogHard_hybrid( A,B,m ) is the syntax.
%  A,B are (m+n)x(m+n) adjacency matrices, 
% loops/multiedges/directededges allowed.
% m : the number of hard seeds to try (if it is a vector, all the values
% are used in turn)

% Uses rQAP/rQAP2 (sucessive relaxation) formulation of D.  Fishkind for SGM
% which is a minimization problem as opposed max. for rQAP (aka FAQ)
% It is assumed that the first m vertices of A's graph
% correspond respectively to the first m vertices of B's graph,
% corr gives the vertex correspondences  
% For example, corr=[ 1 2 7 16 30 ...
% means that the vtx1ofA-->vtx1ofB, 2-->2, 3-->7, 4-->16, 5-->30 
%  example: EXECUTE the following:
% >> v=[ [1:5] 5+randperm(400)]; B=round(rand(405,405));A=B(v,v);
% >> [corr,P] = ConVogHard_hybrid( A,B,5 ) ; [v; corr]
% ready Sep 24, 2012   (Sancar's code)


[totv,~]=size(A);
n=totv-m;

A12=A(1:m,m+1:m+n);
A21=A(m+1:m+n,1:m);
A22=A(m+1:m+n,m+1:m+n);
B12=B(1:m,m+1:m+n);
B21=B(m+1:m+n,1:m);
B22=B(m+1:m+n,m+1:m+n);


patience_1=20;
patience_2=20;
patience=patience_2;
tol_init = 1E-3;
tol=1E-3;
epsilon=0.01;
P=ones(n,n)/n;
toggle=1;
iter=0;
fvals = zeros(patience+1,4);
alpha_vals =zeros(1,patience+1);
r=1;

% Hack to overcome the lack of converge of rqap2 part of gradient
if (m<20) 
cutoff = patience_1;
else 
cutoff= 15- floor(m/27);
end

while (toggle==1)&&(iter<patience_1)
    
    iter=iter+1;
    r= 0.5- atan((iter-(patience/2)))/pi;
	if (iter>cutoff)     
	r=0;
	end    


    Grad=-2*A21*B21'-2*A12'*B12+2*r*P*(B21*B21')+2*r*(A12'*A12)*P+2*r*((A22'*A22*P)+ ...
        P*(B22*B22')) -2*( A22'*P*B22+A22*P*B22');
   % Grad=-2*A21*B21' - 2*A12'*B12 ...
   %      -2* A22'*P*B22 - 2*A22*P*B22';
    %Grad is the gradient of the function to be minimized
    ind=lapjv(Grad,0.01);
    T=eye(n);
    T=T(ind,:);
    tempmat1 = r*(B21*B21'+B22*B22');
    tempmat2 = r*(A12'*A12+A22'*A22);
   %c=trace(P'*P*tempmat1)+trace(tempmat2*(P*P'))  ...
   %-trace(P'*A22'*P*B22)-trace(P'*A22*P*B22');
      % c= -trace(P'*A22'*P*B22) - trace(P'*A22*P*B22');
      c=trace(P'*P*tempmat1)+trace(tempmat2*(P*P'))  ...
   -2*trace(P'*A22'*P*B22);
      



   % d=trace((T'*P+P'*T)*tempmat1)+trace(tempmat2*(T*P'+P*T'))- ...
   %    trace(P'*(A22'*T*B22+(A22*T*B22')))-trace(T'*(A22'*P*B22+(A22*P*B22')));
    
    
   % d= -trace(P'*(A22'*T*B22+(A22*T*B22'))) - trace(T'*(A22'*P*B22+(A22*P*B22')));
   d=trace((T'*P+P'*T)*tempmat1)+trace(tempmat2*(T*P'+P*T'))- ...
       2*trace(P'*(A22*T*B22'))-2*trace(T'*(A22*P*B22'));
    
    
    
   % e=trace(T'*T*tempmat1)+trace(tempmat2*(T*T'))- ...
   %    trace(T'*A22'*T*B22)-trace(T'*A22*T*B22');
   % e= - trace(T'*A22'*T*B22) - trace(T'*A22*T*B22');
   e=trace(T'*T*tempmat1)+trace(tempmat2*(T*T'))- ...
       2*trace(T'*A22'*T*B22);
   
   
   
    u= -2 * trace(P'*A21*B21'+P'*A12'*B12);
    v= -2 * trace(T'*A21*B21'+T'*A12'*B12);
   
    alpha=-(d-2*e+u-v)/(2*(c-d+e));
    f0=e+v;
    f1=c+u;
    %Function value at alpha
    falpha=(c-d+e)*alpha^2+(d-2*e+u-v)*alpha+e+v;
    %current function value (value at alpha=0)
    %fvals(iter,1)=f1;
    %fvals(iter,2)=falpha;
    %f0:  function value at alpha=0
    %fvals(iter,3)=f0;
    alpha_vals(iter+1) = alpha;
    %fvals(iter,4)= norm(Grad,2); 
    newfval= falpha;
    if ((alpha>0)&&(falpha<f0)&&(falpha<f1))
        P=alpha*P+(1-alpha)*T;
         newfval=falpha;
    elseif (f0<(f1-epsilon))
        P=T;
        newfval=f0;
    else
        toggle=0;
    end
   if ((abs(newfval-f1)/abs(f1))<tol_init) 
        toggle=0;    
    end
end
alpha_vals;
corr=lapjv(-P,0.01);
corr=[ 1:m,  m+corr];


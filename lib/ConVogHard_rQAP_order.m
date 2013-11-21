function [ corr,iter_final ,max_fvals,fval,fval_after_proj,P_whole,P_proj] = ConVogHard_rQAP_order( A,B,m,ordering,useJV )

% [ corr,iter ,max_fvals,fval,fval_after_proj,P_whole,P_proj] = ConVogHard_rQAP_order( A,B,m,ordering,useJV ) 
% is the syntax.
%  A,B are (m+n)x(m+n) adjacency matrices, 
% loops/multiedges/directededges allowed.
% m : the number of hard seeds to try (if it is a vector, all the values
% are used in turn)
% ordering : the indices of vertices such that the first m  are used as 
% seeds for m number of seeds and the remaining n vertices are matched (*)
% useJV : 1 if the faster JV algorithm is used for the linear assignment
% subproblem
%
% (*) It is assumed that the  m vertices  of A's graph listed in the vector
% ordering(1:m) correspond respectively to the m vertices  of B's graph with the same indices,
% the seeds are  the vertices in A's graph and in B's graph with indices in  ordering(1:m).
% corr gives the vertex correspondences 
% For example, corr=[ 1 2 7 16 30 ...
% means that the vtx_ordering(1)_of_A-->vtx_1_of_B, ordering(2)-->2, ordering(3)-->7, ordering(4)-->16, ordering(5)-->30 
% that is, the vertex of A with the index ordering(i) is matched to vertex of B with
% the index corr(i)
%  example: EXECUTE the following:
% >> v=[ [1:5] 5+randperm(400)]; B=round(rand(405,405));A=B(v,v);
% >> [corr,P] = ConVogHard_rQAP_order( A,B,5,1:(m+n),1 ) ; [v; corr]
% ready June 1, 2012   (Donniell's code)


[totv,~]=size(A);
n=totv-m;

A12=A(ordering(1:m),ordering(m+1:m+n));
A21=A(ordering(m+1:m+n),ordering(1:m));
A22=A(ordering(m+1:m+n),ordering(m+1:m+n));
B12=B(ordering(1:m),ordering(m+1:m+n));
B21=B(ordering(m+1:m+n),ordering(1:m));
B22=B(ordering(m+1:m+n),ordering(m+1:m+n));

%Maximum number of iterations
patience=50;


% change ratio in function value for stopping criterion
%If change ratio is lower, the iterations terminate
tol=1E-3;
epsilon=0.01;
P=ones(n,n)/n;
toggle=1;
iter=0;
fvals = zeros(patience+1,4);
alpha_vals =zeros(1,patience+1);
while (toggle==1)&(iter<patience)
    iter=iter+1;
    Grad=A22*P*B22'+A22'*P*B22+A21*B21'+A12'*B12;
    ind=[];
    if useJV 
        ind=lapjv(-Grad,0.01);
    else
        ind=YiCaoHungarian(-Grad);
    end
    T=eye(n);
    T=T(ind,:);
    c=trace(A22'*P*B22*P');
    d=trace(A22'*T*B22*P')+trace(A22'*P*B22*T');
    e=trace(A22'*T*B22*T');
    u=trace(P'*A21*B21'+P'*A12'*B12);
    v=trace(T'*A21*B21'+T'*A12'*B12);
    alpha=-(d-2*e+u-v)/(2*(c-d+e));
    f0=e+v;
    f1=c+u;
    %Function value at alpha
    falpha=(c-d+e)*alpha^2+(d-2*e+u-v)*alpha+e+v;
    %current function value (value at alpha=0)
     fvals(iter,1)=f1;
    fvals(iter,2)=falpha;
    %f0:  function value at alpha=0
    fvals(iter,3)=f0;
    fvals(iter,4)= norm(Grad,2);
    alpha_vals(iter+1) = alpha;
%     if (isnan(falpha))
%     Grad
%     c-d+e
%     alpha
%     end
    newfval=0;
    if ((alpha>0)&&(falpha>f0)&&(falpha>f1))
        P=alpha*P+(1-alpha)*T;
         newfval=falpha;
    elseif (f0>(f1+epsilon))
        P=T;
        newfval=f0;
    else
        toggle=0;
    end
    if ((abs(newfval-f1)/abs(f1))<tol) 
        toggle=0;    
    end
end
alpha_vals;
corr=[];
    if useJV 
        corr=lapjv(-P,0.01);
    else
        corr=YiCaoHungarian(-P);
        
    end
%corr=lapjv(-P,0.01);
%corr=lapjv(-P,0.01);


corr=[ ordering(1:m),  ordering(m+corr)];
P_whole=[ eye(m) zeros(m,n) ; zeros(n,m) P];
max_fvals= fvals;
norm_mat = A(ordering,ordering)*P_whole-P_whole*B(ordering,ordering);
fval=sum(sum(abs(norm_mat)));


 P_proj=eye(m+n);
 P_proj=P_proj(corr,:);
norm_mat_proj = A*P_proj-P_proj*B;
fval_after_proj= sum(sum(abs(norm_mat_proj)));

iter_final= iter;

end


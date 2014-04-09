function [truematch, truematch_ratio] = PermMat2TrueMatch(P, truealignment, m)

[n,~] =  size(P);
mplusn = m+n; 
n_vec = (1:n)';
temp=P*n_vec;
alignsoln=[ 1:m temp'+m ];
truematch = sum(truealignment(((m+1):mplusn)) == ...
        alignsoln((m+1):mplusn));
truematch_ratio = truematch/(n);
    

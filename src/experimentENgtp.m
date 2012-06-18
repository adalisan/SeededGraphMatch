function [meanperc,stdevperc] = experimentENgtp(k,m,n,p,q,r,numiter)

% syntax is: [meanperc,stdevperc] = experimentENgtp(k,m,n,p,q,r,numiter)
%
% There are k + m + n employees at Enron
% k + m are EnronExecutives, and the other n are EnronWorkerBees.
% We know the identities of k EnronExecutives and don't know anyone else.
% We want to identify the m unknown EnronExecutives.
%
% probability of communication between any two given Execs is p
% probability of communication between any two given WorkerBees is q
% prob of comm between any given Exec and any given WorkerBee is r
% (all pairs independent)
%
% the percentage of the m unknown Executives identified with 
% graphmatch AGAINST BEROULLI PARAMETERS using id'ed Execs as hard seeds 
% (ie graphmatched to MATRIX OF BERNOULLI PARAMETERS)
% is summarized with mean and standard deviation  "[meanperc,stdevperc]"
% over numiter such experiments.
%
% For example, execute 
% [meanperc,stdevperc] = experimentENgtp(50,50,100,.7,.4,.4,30)
% (note: this will take computer 2 to 3 min)


perc=zeros(1,numiter);
for i=1:numiter
    A11=rand(k+m,k+m)<p;
    A11=A11-triu(A11); A11=A11+A11';
    B11=p*ones(k+m,k+m)-p*eye(k+m);
    A22=rand(n,n)<q;
    A22=A22-triu(A22); A22=A22+A22';
    B22=q*ones(n,n)-q*eye(n);
    A12=rand(k+m,n)<r;
    B12=r*ones(k+m,n);
    A=[A11 A12; A12' A22];
    B=[B11 B12; B12' B22];
    
    confound=[1:k,k+randperm(m+n)];
    II=eye(k+m+n);
    Pconf=II(confound,:);
    A=Pconf*A*Pconf';
    
    bij=seedgraphmatchell2(A,B,k);   
    
    bij=bij*Pconf;
    
    for j=k+1:k+m
       if (k+1<=bij(j))&(bij(j)<=k+m)
           perc(i)=perc(i)+1/m;
       end
    end
end
meanperc=mean(perc);
stdevperc=std(perc);



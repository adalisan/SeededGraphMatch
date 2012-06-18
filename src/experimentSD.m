function [pc] = experimentSD(n,maxm,numiter)

pc=zeros(1,maxm+1);
for i=1:numiter
    i
    Bernoulli=rand(maxm+n,maxm+n);
    A=rand(maxm+n,maxm+n)<Bernoulli;
    A=A-triu(A);A=A+A';
    B=rand(maxm+n,maxm+n)<Bernoulli;
    B=B-triu(B);B=B+B';
    for j=0:maxm
        At=A(maxm-j+1:maxm+n,maxm-j+1:maxm+n);
        Bt=B(maxm-j+1:maxm+n,maxm-j+1:maxm+n);
        confound=[1:j,j+randperm(n)];Pconf=eye(n+j);Pconf=Pconf(confound,:); Bt=Pconf*Bt*Pconf';
        bij=seedgraphmatchell2(At,Bt,j);
        bij=bij*Pconf';
        pc(j+1)=pc(j+1)+sum(bij(j+1:j+n)==[j+1:j+n])/n;
    end
end
pc=pc/numiter;
close
plot([0:maxm],pc,'r*')



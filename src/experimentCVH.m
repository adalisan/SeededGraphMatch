function [pc,pc_rQAP2,pc_hybrid,iterofFW_rQAP,iterofFW_rQAP2,iterofFW_hybrid] = experimentCVH(n,maxm,numiter)

pc=zeros(1,maxm+1);
pc_rQAP2=zeros(1,maxm+1);
pc_hybrid=zeros(1,maxm+1);
iterofFW_rQAP =zeros(numiter,maxm+1);
iterofFW_rQAP2 =zeros(numiter,maxm+1);
iterofFW_hybrid =zeros(numiter,maxm+1);
fvals_rqap2 = zeros(26,3);
fvals_rqap = zeros(26,3);
fvals_hybrid = zeros(26,3);
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
        [bij,iterofFW_rQAP(i,(j+1)),fvals_rqap]=ConVogHard_rQAP(At,Bt,j);
        fvals_rqap(1:(iterofFW_rQAP(i,(j+1))+1),:)
        pc(j+1)=pc(j+1)+sum(bij(j+1:j+n)==[j+1:j+n])/n;
        
        %rQAP2 forumulation
         [cij,iterofFW_rQAP2(i,(j+1)),fvals_rqap2 ]=ConVogHard_rQAP2(At,Bt,j);
        pc_rQAP2(j+1)=pc_rQAP2(j+1)+sum(cij(j+1:j+n)==[j+1:j+n])/n;
        fvals_rqap2 (1:(iterofFW_rQAP2(i,(j+1))+1),:)
        
         [dij,iterofFW_hybrid(i,(j+1)),fvals_hybrid ]= ...
             ConVogHard_hybrid(At,Bt,j);
        pc_hybrid(j+1)=pc_hybrid(j+1)+sum(dij(j+1:j+n)==[j+1:j+n])/n;
        fvals_hybrid (1:(iterofFW_hybrid(i,(j+1))+1),:)
        
    end
    %save
end
pc=pc/numiter;
pc_rQAP2=pc_rQAP2/numiter;

pc_hybrid=pc_hybrid/numiter;
close
plot([0:maxm],pc,'r*')
hold on
plot([0:maxm],pc_rQAP2,'b.')
hold on 
plot([0:maxm],pc_hybrid,'g+')
title('rQAP vs rQAP2 vs hybrid')
xlabel('Number of Hard seeds')
ylabel('Fraction of correct matches')

mean(iterofFW_rQAP,1)
mean(iterofFW_rQAP2,1)
mean(iterofFW_hybrid,1)





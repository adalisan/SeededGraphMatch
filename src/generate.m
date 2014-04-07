function [A,B,truth]=generate(Bern,crln,m)

% [A,B,truth]=generate(Bern,crln,m)
% generates two random correlated Bernoulli graphs with 
% Bernoulli parameters Bern, correlation crln, and the first m vtc are seeds.
% Then the adjacency matrices are vertex-scrambled, with "truth" being true alignment  
% DEF Feb 18, 2014, ready

[nplusm,~]=size(Bern);
n=nplusm-m;
% 
% A=zeros(nplusm,nplusm);
% 
% B=zeros(nplusm,nplusm);
% for i=1:nplusm
%     for j=i+1:nplusm
%         Bern_tmp = double(Bern(i,j));
%         A(i,j)=    (rand<Bern_tmp);
%         A(j,i)=A(i,j);
%         A_tmp= double(A(i,j));
%         
%         B(i,j)=    (rand< (  (1.0-crln)*Bern_tmp+crln*A_tmp )  );
%         B(j,i)=B(i,j);
%     end
% end

 Bernoulli=rand(N);
        A=rand(N)<Bernoulli;
        A=A-triu(A);A=A+A';
        B=A;
        B=bitflip(A,rho_val);
        B=B-triu(B);B=B+B';



mix=[ [1:m] m+randperm(n) ];
A=A(mix,mix);
B=B(mix,mix);
truth=[ [1:m] m+randperm(n) ];
A=A(truth,truth);



function [A,B,truth]=generate(Bern,crln,m)

% [A,B,truth]=generate(Bern,crln,m)
% generates two random correlated Bernoulli graphs with 
% Bernoulli parameters Bern, correlation crln, and the first m vtc are seeds.
% Then the adjacency matrices are vertex-scrambled, with "truth" being true alignment  
% DEF Feb 18, 2014, ready

[nplusm,~]=size(Bern);
n=nplusm-m;

A=zeros(nplusm,nplusm);
B=zeros(nplusm,nplusm);
for i=1:nplusm
    for j=i+1:nplusm
        A(i,j)=    (rand<Bern(i,j));
        A(j,i)=A(i,j);
        A_temp = double(A);
        Bern_tmp = double(Bern)
        B(i,j)=    (rand< (  (1.0-crln)*Bern(i,j)+crln*A(i,j) )  );
        B(j,i)=B(i,j);
    end
end


mix=[ [1:m] m+randperm(n) ];
A=A(mix,mix);
B=B(mix,mix);
truth=mix;
A=A(truth,truth);



function [ Bmat ] = bitflip( Amat,q_prob )
%UNTITLED This function flips(0 to 1, or 1 to 0) each element 
%of Amat with probability q
n=size(Amat,1)
Bmat= double( xor(rand(n,n)<q_prob,Amat));
end


function [alignment]= PermMat2TrueMatch(P,true_alignment)

[mplusn,~] =  size(P);
temp=P*[1:mplusn]';
alignment=[ [1:m] temp'+m ];
mplusn = len(alignment);
sum(true_alignment[(m+1):mplusn] == alignment[(m+1):mplusn])

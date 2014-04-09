function [truematch,truematch_ratio]= PermMat2TrueMatch(P,true_alignment,m)
    [n,~] =  size(P);
    mplusn = m+n; 
    temp=P*[1:mplusn]';
    alignment=[ [1:m] temp'+m ];
    mplusn = len(alignment);
    truematch = sum(true_alignment[(m+1):mplusn] == alignment[(m+1):mplusn]);
    truematch_ration = truematch/(n)
    

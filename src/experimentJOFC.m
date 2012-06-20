%function [pc,pcjofc,pcjofc_diff,iterofFW_rQAP] = experimentJOFC(n,maxm,numiter)

pc=zeros(maxm+1,3);
pcjofc = zeros(maxm+1,3);
pcjofc_diff = zeros(maxm+1,3);

truematch = zeros(3,1);
truematch_diff = zeros(3,1);

iterofFW_rQAP =zeros(numiter,maxm+1);


fvals_rqap = zeros(26,3);


w_wals_vec = zeros(3,1)
w_vals_vec(1) = 0.5 
w_vals_vec(2) = 0.8
w_vals_vec(3) = 0.95

for i=1:numiter
    i
    Bernoulli=rand(maxm+n,maxm+n);
    A=rand(maxm+n,maxm+n)<Bernoulli;
    A=A-triu(A);A=A+A';
   % B=rand(maxm+n,maxm+n)<Bernoulli;
    B=bitflip(A,q);
    B=B-triu(B);B=B+B';
   
    
    for j=0:3:maxm
        At=A(maxm-j+1:maxm+n,maxm-j+1:maxm+n);
        Bt=B(maxm-j+1:maxm+n,maxm-j+1:maxm+n);
        [bij,iterofFW_rQAP(i,(j+1)),fvals_rqap]=ConVogHard_rQAP(At,Bt,j);
        fvals_rqap(1:(iterofFW_rQAP(i,(j+1))+1),:);
        pc(j+1)=pc(j+1)+sum(bij(j+1:j+n)==[j+1:j+n])/n;
        
        
        embed.dim =2;
        if (j<=2)
            continue
        else
            embed.dim = min(j,5) ;
      
        end
        
        putRdata('At',At)
        putRdata('Bt',Bt)
        insample_logic_vector = logical([ones(1,j) zeros(1,n) ]');
        putRdata('insample_logic_vec', insample_logic_vector)
        putRdata('embed.dim', embed.dim)
        
        
        putRdata('w_vals_vec',w_vals_vec)
        evalR('insample_logic_vec<-as.vector(insample_logic_vec)')
        evalR('insample_logic_vec<-rep(insample_logic_vec,2)')
        evalR('sink("debug.mat.txt")')
        %evalR('traceback()')
        evalR('print(insample_logic_vec)')
        evalR('sink()')
        %evalR('error.handle <- function()')
        evalR(['error.handle <- function(ex) {print(ex)}'])
        

        evalR('jofc.result <- try(JOFC.graph.custom.dist(G=At,Gp=Bt, in.sample.ind=insample_logic_vec,  d.dim=embed.dim,w.vals.vec=w_vals_vec,graph.is.directed=FALSE))')
        % evalR('eval(parse(jofc.call))')
        %evalR('jofc.result<- tryCatch ({eval(parse(jofc.call))},error=error.handle)')
        %  traceback   })')
        evalR('sink("debug.matlab.txt")')
        evalR('traceback()')
        evalR('print(jofc.result)')
        evalR('sink()')
        evalR('jofc.res.1<-jofc.result[[1]]')
        evalR('jofc.res.2<-jofc.result[[2]]')
        evalR('jofc.res.3<-jofc.result[[3]]')
        jofc_dist_mat_1=getRdata('jofc.res.1') 
        jofc_dist_mat_2=getRdata('jofc.res.2') 
        jofc_dist_mat_3=getRdata('jofc.res.3') 
        matching_1 = YiCaoHungarian(jofc_dist_mat_1)
        matching_2 = YiCaoHungarian(jofc_dist_mat_2)
        matching_3 = YiCaoHungarian(jofc_dist_mat_3)
        evalR('M.result.1<-try(solveMarriage(jofc.res.1))')
        evalR('M.result.2<-try(solveMarriage(jofc.res.2))')
        evalR('M.result.3<-try(solveMarriage(jofc.res.3))')
        
        
%         evalR('skip.iter<-FALSE')
%         evalR('	skip.iter <-inherits(M.result.1,"try-error")')
%         if (getRdata('skip.iter'))
%             'Skipping iteration'
%             continue
%             
%         end
%         
%         evalR('NumofTruePairing.1<-present(M.result.1)')
%         evalR('NumofTruePairing.2<-present(M.result.2)')
%         evalR('NumofTruePairing.3<-present(M.result.3)')
%         
%         
%         truematch(1) = getRdata('NumofTruePairing.1');
%         
%         truematch(2) = getRdata('NumofTruePairing.2');
%         
%         truematch(3) = getRdata('NumofTruePairing.3');
        
        truematch(1) = sum(matching_1==1:n);
        
        truematch(2) = sum(matching_2==1:n)
        truematch(3) = sum(matching_3==1:n);
        
        evalR('jofc.diff.dist.result <- try(JOFC.graph.diff(G=At,Gp=Bt, in.sample.ind=insample_logic_vec,  d.dim=embed.dim,w.vals.vec=w_vals_vec,T.param=2))')
        
        evalR('sink("debug.matlab.txt")')
        evalR('traceback()')
        evalR('print(jofc.diff.dist.result )')
        evalR('sink()')
        
        evalR('jofc.res.diff.1<-jofc.diff.dist.result [[1]]')
        evalR('jofc.res.diff.2<-jofc.diff.dist.result [[2]]')
        evalR('jofc.res.diff.3<-jofc.diff.dist.result [[3]]')
         jofc_dist_mat_1_diff=getRdata('jofc.res.diff.1') ;
        jofc_dist_mat_2_diff=getRdata('jofc.res.diff.2') ;
        jofc_dist_mat_3_diff=getRdata('jofc.res.diff.3') ;
        matching_1_diff = YiCaoHungarian(jofc_dist_mat_1_diff)
        matching_2_diff = YiCaoHungarian(jofc_dist_mat_2_diff)
        matching_3_diff = YiCaoHungarian(jofc_dist_mat_3_diff)
      
%         evalR('M.result.diff.1<-try(solveMarriage(jofc.res.diff.1))')
%         evalR('M.result.diff.2<-try(solveMarriage(jofc.res.diff.2))')
%         evalR('M.result.diff.3<-try(solveMarriage(jofc.res.diff.3))')
%         evalR('skip.iter<-FALSE')
%         evalR('	skip.iter <-inherits(M.result.diff.1,"try-error")')
% %         if (getRdata('skip.iter'))
%             'Skipping iteration'
%             continue
%             
% %         end
%         evalR('NumofTruePairing.diff.1<-present(M.result.diff.1)')
%         evalR('NumofTruePairing.diff.2<-present(M.result.diff.2)')
%         evalR('NumofTruePairing.diff.3<-present(M.result.diff.3)')
%         
%         
%         truematch_diff(1) = getRdata('NumofTruePairing.diff.1');
%         
%         truematch_diff(2) = getRdata('NumofTruePairing.diff.2');
%         
%         truematch_diff(3) = getRdata('NumofTruePairing.diff.3');
%                 
         truematch_diff(1) = sum(matching_1_diff==1:n);
        
        truematch_diff(2) = sum(matching_2_diff==1:n);
        truematch_diff(3) = sum(matching_3_diff==1:n)
    
        
        pcjofc(j+1,1) =  pcjofc(j+1,1) + truematch(1)/n;
        pcjofc(j+1,2) =  pcjofc(j+1,2) + truematch(2)/n;
        pcjofc(j+1,3) =  pcjofc(j+1,3) + truematch(3)/n;
        pcjofc_diff(j+1,1) =  pcjofc_diff(j+1,1) + truematch_diff(1)/n;
        pcjofc_diff(j+1,2) =  pcjofc_diff(j+1,2) + truematch_diff(2)/n;
        pcjofc_diff(j+1,3) =  pcjofc_diff(j+1,3) + truematch_diff(3)/n;
     
    end
end
pc=pc/numiter;
pcjofc =pcjofc /numiter;
pcjofc_diff =pcjofc_diff /numiter;












pc=zeros(maxm+1,3);
pcjofc = zeros(maxm+1,3);
pcjofc_diff = zeros(maxm+1,3);
pcjofc_custom = zeros(maxm+1,3);

truematch = zeros(3,1);
truematch_diff = zeros(3,1);
truematch_custom_dissimilar  = zeros(3,1);
iterofFW_rQAP =zeros(numiter,maxm+1);


fvals_rqap = zeros(26,3);

num_w =1;
w_wals_vec = zeros(num_w,1);
w_vals_vec(1) = 0.8;
%w_vals_vec(1) = 0.5 ;
%w_vals_vec(2) = 0.8;
%w_vals_vec(3) = 0.95;
q=0.1

for i=1:numiter
    i
    Bernoulli=rand(maxm+n,maxm+n);
    A=rand(maxm+n,maxm+n)<Bernoulli;
    A=A-triu(A);A=A+A';
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
        elseif (j<10)
            embed.dim=2;
        else
            embed.dim = min(j,5) ;
      
        end
        
        putRdata('At',At)
        putRdata('Bt',Bt)
        ret_At= getRdata('Bt');
        insample_logic_vector = logical([ones(1,j) zeros(1,n) ]');
        putRdata('insample_logic_vec', insample_logic_vector)
        putRdata('embed.dim', embed.dim)
        
        
        putRdata('w_vals_vec',w_vals_vec)
        evalR('insample_logic_vec<-as.vector(insample_logic_vec)')
        evalR('insample_logic_vec<-rep(insample_logic_vec,2)')
        evalR('sink("debug.mat.txt")')
        %evalR('traceback()')
        evalR('print(insample_logic_vec)')
        %evalR('print(At)')
        %evalR('print(Bt)')
        evalR('sink()')
        %evalR('error.handle <- function()')
        evalR(['error.handle <- function(ex) {print(ex)}'])
         evalR('jofc.result <- try(JOFC.graph(G=At,Gp=Bt, in.sample.ind=insample_logic_vec,  d.dim=embed.dim,w.vals.vec=w_vals_vec,graph.is.directed=FALSE))')
       % evalR('eval(parse(jofc.call))')
        %evalR('jofc.result<- tryCatch ({eval(parse(jofc.call))},error=error.handle)')
        %  traceback   })')
        evalR('sink("debug.matlab.txt")')
        evalR('traceback()')
        evalR('print(jofc.result)')
        evalR('sink()')
        evalR('jofc.res.1<-jofc.result[[1]]')
%        evalR('jofc.res.2<-jofc.result[[2]]')
 %       evalR('jofc.res.3<-jofc.result[[3]]')
        getRdata('jofc.res.1') ;
        evalR('M.result.1<-try(solveMarriage(jofc.res.1))')
 %       evalR('M.result.2<-try(solveMarriage(jofc.res.2))')
 %       evalR('M.result.3<-try(solveMarriage(jofc.res.3))')
        
        
        evalR('skip.iter<-FALSE')
        evalR('	skip.iter <-inherits(M.result.1,"try-error")')
        if (getRdata('skip.iter'))
            'Skipping iteration'
            continue
            
        end
        
        evalR('NumofTruePairing.1<-present(M.result.1)')
 %       evalR('NumofTruePairing.2<-present(M.result.2)')
 %       evalR('NumofTruePairing.3<-present(M.result.3)')
        
        
        truematch(1) = getRdata('NumofTruePairing.1');
        
 %       truematch(2) = getRdata('NumofTruePairing.2');
        
%        truematch(3) = getRdata('NumofTruePairing.3');
        
           
       evalR('jofc.diff.dist.result <- try(JOFC.graph.diff(G=At,Gp=Bt, in.sample.ind=insample_logic_vec,  d.dim=embed.dim,w.vals.vec=w_vals_vec,T.param=2))')
      
        evalR('sink("debug.matlab.txt")')
        evalR('traceback()')
        evalR('print(jofc.diff.dist.result )')
        evalR('sink()')
        
        evalR('jofc.res.diff.1<-jofc.diff.dist.result [[1]]')
 %       evalR('jofc.res.diff.2<-jofc.diff.dist.result [[2]]')
 %       evalR('jofc.res.diff.3<-jofc.diff.dist.result [[3]]')
%        getRdata('jofc.res.diff.3') ;
        evalR('M.result.diff.1<-try(solveMarriage(jofc.res.diff.1))')
%        evalR('M.result.diff.2<-try(solveMarriage(jofc.res.diff.2))')
%        evalR('M.result.diff.3<-try(solveMarriage(jofc.res.diff.3))')
        evalR('skip.iter<-FALSE')
        evalR('	skip.iter <-inherits(M.result.diff.1,"try-error")')
        if (getRdata('skip.iter'))
            'Skipping iteration'
            continue
            
        end
        evalR('NumofTruePairing.diff.1<-present(M.result.diff.1)')
 %       evalR('NumofTruePairing.diff.2<-present(M.result.diff.2)')
%        evalR('NumofTruePairing.diff.3<-present(M.result.diff.3)')
        
        
        truematch_diff(1) = getRdata('NumofTruePairing.diff.1');
        
 %       truematch_diff(2) = getRdata('NumofTruePairing.diff.2');
        
%        truematch_diff(3) = getRdata('NumofTruePairing.diff.3');
                
        pcjofc(j+1,1) =  pcjofc(j+1,1) + truematch(1)/n;
 %       pcjofc(j+1,2) =  pcjofc(j+1,2) + truematch(2)/n;
 %       pcjofc(j+1,3) =  pcjofc(j+1,3) + truematch(3)/n;
        pcjofc_diff(j+1,1) =  pcjofc_diff(j+1,1) + truematch_diff(1)/n;
%        pcjofc_diff(j+1,2) =  pcjofc_diff(j+1,2) + truematch_diff(2)/n;
 %       pcjofc_diff(j+1,3) =  pcjofc_diff(j+1,3) + truematch_diff(3)/n;
 
 
        evalR('jofc.result.custom <- try(JOFC.graph.custom.dist(G=At,Gp=Bt, in.sample.ind=insample_logic_vec,  d.dim=embed.dim,w.vals.vec=w_vals_vec,graph.is.directed=FALSE))')
        % evalR('eval(parse(jofc.call))')
        %evalR('jofc.result<- tryCatch ({eval(parse(jofc.call))},error=error.handle)')
        %  traceback   })')
        evalR('sink("debug.matlab.txt")')
        evalR('traceback()')
        evalR('print(jofc.result.custom)')
        evalR('sink()')
        evalR('jofc.res.1.custom<-jofc.result.custom[[1]]')
        evalR('jofc.res.2.custom<-jofc.result.custom[[2]]')
        evalR('jofc.res.3.custom<-jofc.result.custom[[3]]')
        jofc_dist_mat_1_custom_dissimilar=getRdata('jofc.res.1.custom') 
        jofc_dist_mat_2_custom_dissimilar=getRdata('jofc.res.2.custom') 
        jofc_dist_mat_3_custom_dissimilar=getRdata('jofc.res.3.custom') 
        matching_1_custom_dissimilar = YiCaoHungarian(jofc_dist_mat_1_custom_dissimilar)
        matching_2_custom_dissimilar = YiCaoHungarian(jofc_dist_mat_2_custom_dissimilar)
        matching_3_custom_dissimilar = YiCaoHungarian(jofc_dist_mat_3_custom_dissimilar)
%         evalR('M.result.1<-try(solveMarriage(jofc.res.1.custom))')
%         evalR('M.result.2<-try(solveMarriage(jofc.res.2.custom))')
%         evalR('M.result.3<-try(solveMarriage(jofc.res.3.custom))')
%         
        
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
        
        truematch_custom_dissimilar(1) = sum(matching_1_custom_dissimilar==1:n);
        
        truematch_custom_dissimilar(2) = sum(matching_2_custom_dissimilar==1:n);
        truematch_custom_dissimilar(3) = sum(matching_3_custom_dissimilar==1:n);
 
 
        pcjofc_custom(j+1,1) =  pcjofc_custom(j+1,1) + truematch_custom_dissimilar(1)/n;
         pcjofc_custom(j+1,2) =  pcjofc_custom(j+1,1) + truematch_custom_dissimilar(2)/n;
 
          pcjofc_custom(j+1,3) =  pcjofc_custom(j+1,1) + truematch_custom_dissimilar(3)/n;
 
 
 
 
 
 
        
    end
end
pc=pc/numiter;
pcjofc =pcjofc /numiter;
pcjofc_diff =pcjofc_diff /numiter;










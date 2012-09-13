
q= 0: 0.1:0.3;
q_len = length(q);
pc=zeros(maxm+1,q_len);
pcjofc_sp = zeros(maxm+1,q_len);
pcjofc_diff = zeros(maxm+1,q_len);
pcjofc_custom = zeros(maxm+1,q_len);

pcjofc_dice =  zeros(maxm+1,q_len);
pcjofc_jaccard = zeros(maxm+1,q_len);

pcjofc_invlogweighted = zeros(maxm+1,q_len);
pcjofc_ell1 = zeros(maxm+1,q_len);



truematch_shortest_path = zeros(q_len,1);
truematch_diffusion = zeros(q_len,1);
truematch_custom_dissimilar  = zeros(q_len,1);
truematch_dice = zeros(q_len,1);
truematch_jaccard = zeros(q_len,1);
truematch_invlogweighted= zeros(q_len,1);




iterofFW_rQAP =zeros(numiter,maxm+1);


fvals_rqap = zeros(26,q_len);

num_w =1;
w_wals_vec = zeros(num_w,1);
w_vals_vec(1) = 0.9;
%w_vals_vec(1) = 0.5 ;
%w_vals_vec(2) = 0.8;
%w_vals_vec(3) = 0.95;


for q_i= 1:length(q)
    
    q_val=q(q_i)
    for i=1:numiter
        i
        Bernoulli=rand(maxm+n,maxm+n);
        A=rand(maxm+n,maxm+n)<Bernoulli;
        A=A-triu(A);A=A+A';
        B=A;
        B=bitflip(A,q_val);
        B=B-triu(B);B=B+B';
        for j=0:3:maxm
            At=A(maxm-j+1:maxm+n,maxm-j+1:maxm+n);
            Bt=B(maxm-j+1:maxm+n,maxm-j+1:maxm+n);
            [bij,iterofFW_rQAP(i,(j+1)),fvals_rqap]=ConVogHard_rQAP(At,Bt,j);
            fvals_rqap(1:(iterofFW_rQAP(i,(j+1))+1),:);
            pc(j+1,q_i)=pc(j+1,q_i)+sum(bij(j+1:j+n)==[j+1:j+n]);
            
            
            embeddim=2;
            
            
            if (j<=5)
                continue
            elseif (j<12)
                embeddim=5;
            else
                embeddim =  10 ;
                
            end
            
            putRdata('At',At)
            putRdata('Bt',Bt)
            ret_At= getRdata('Bt');
            insample_logic_vector = logical([ones(1,j) zeros(1,n) ]');
            putRdata('insample_logic_vec', insample_logic_vector)
            putRdata('embed.dim', embeddim)
            
            
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
            
            
            
            evalR('jofc.result.shortest.path <- try(JOFC.graph.custom.dist(G=At,Gp=Bt, in.sample.ind=insample_logic_vec,  d.dim=embed.dim,w.vals.vec=w_vals_vec,graph.is.directed=FALSE,vert_diss_measure="default"))')
            evalR('sink("debug.matlab.txt")')
            evalR('traceback()')
            evalR('print(jofc.result.shortest.path)')
            evalR('sink()')
            
            evalR('jofc.res.1<-as.matrix(jofc.result.shortest.path[[1]])')
            Diss_mat=getRdata('jofc.res.1');
              matching = lapjv(Diss_mat,0.01);
             truematch_shortest_path(q_i) = sum(matching == 1:n);
             
             evalR('jofc.result.diffusion <- try(JOFC.graph.custom.dist(G=At,Gp=Bt, in.sample.ind=insample_logic_vec,  d.dim=embed.dim,w.vals.vec=w_vals_vec,graph.is.directed=FALSE,vert_diss_measure="diffusion"))')
            evalR('jofc.res.1<-jofc.result.diffusion[[1]]')
            Diss_mat=getRdata('jofc.res.1');
              matching = YiCaoHungarian(Diss_mat);
             truematch_diffusion(q_i) = sum(matching == 1:n);
            
             
             
            %%%%Dice dissimilarity result
            
            
            evalR('jofc.result <- try(JOFC.graph.custom.dist(G=At,Gp=Bt, in.sample.ind=insample_logic_vec,  d.dim=embed.dim,w.vals.vec=w_vals_vec,graph.is.directed=FALSE,vert_diss_measure="dice"))')
            
            % evalR('eval(parse(jofc.call))')
            %evalR('jofc.result<- tryCatch ({eval(parse(jofc.call))},error=error.handle)')
            %  traceback   })')
            evalR('sink("debug.matlab.dice.txt")')
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
            
            
            truematch_dice(q_i) = getRdata('NumofTruePairing.1');
            
            %       truematch(2) = getRdata('NumofTruePairing.2');
            
            %        truematch(3) = getRdata('NumofTruePairing.3');
            
            
            %%%%Jaccard dissimilarity result
            
            
            evalR('jofc.diff.dist.result  <- try(JOFC.graph.custom.dist(G=At,Gp=Bt, in.sample.ind=insample_logic_vec,  d.dim=embed.dim,w.vals.vec=w_vals_vec,graph.is.directed=FALSE,vert_diss_measure="jaccard"))')
            
            evalR('sink("debug.matlab.jaccard.txt")')
            evalR('traceback()')
            evalR('print(jofc.diff.dist.result )')
            evalR('sink()')
            
            evalR('jofc.res.diff.1<-jofc.diff.dist.result [[1]]')
            evalR('M.result.diff.1<-try(solveMarriage(jofc.res.diff.1))')
            evalR('skip.iter<-FALSE')
            evalR('	skip.iter <-inherits(M.result.diff.1,"try-error")')
            if (getRdata('skip.iter'))
                'Skipping iteration'
                continue
                
            end
            evalR('NumofTruePairing.diff.1<-present(M.result.diff.1)')
            
            truematch_jaccard(q_i) = getRdata('NumofTruePairing.diff.1');
            
            
            
            
            evalR('jofc.invlogweighted.dist.result <- try(JOFC.graph.custom.dist(G=At,Gp=Bt, in.sample.ind=insample_logic_vec,  d.dim=embed.dim,w.vals.vec=w_vals_vec,graph.is.directed=FALSE,vert_diss_measure="invlogweighted"))')
            
            evalR('sink("debug.matlab.invlog.txt")')
            evalR('traceback()')
            evalR('print(jofc.invlogweighted.dist.result )')
            evalR('sink()')
            
            evalR('jofc.res.invlogweighted.1<-jofc.invlogweighted.dist.result [[1]]')
            evalR('M.result.invlogweighted.1<-try(solveMarriage(jofc.res.invlogweighted.1))')
            evalR('skip.iter<-FALSE')
            evalR('	skip.iter <-inherits(M.result.invlogweighted.1,"try-error")')
            if (getRdata('skip.iter'))
                'Skipping iteration'
                continue
                
            end
            evalR('NumofTruePairing.invlogweighted.1<-present(M.result.invlogweighted.1)')
            
            truematch_invlogweighted(q_i) = getRdata('NumofTruePairing.invlogweighted.1');
            
            
            
            
            
            evalR('jofc.result.custom.ell1 <- try(JOFC.graph.custom.dist(G=At,Gp=Bt, in.sample.ind=insample_logic_vec,  d.dim=embed.dim,w.vals.vec=w_vals_vec,graph.is.directed=FALSE,vert_diss_measure="ell1"))')
            % evalR('eval(parse(jofc.call))')
            %evalR('jofc.result<- tryCatch ({eval(parse(jofc.call))},error=error.handle)')
            %  traceback   })')
            evalR('sink("debug.matlab.ell1.txt")')
            evalR('traceback()')
            evalR('print(jofc.result.custom.ell1)')
            evalR('sink()')
            evalR('jofc.res.1.custom.ell1<-jofc.result.custom.ell1[[1]]')
            %         evalR('jofc.res.2.custom<-jofc.result.custom[[2]]')
            %         evalR('jofc.res.3.custom<-jofc.result.custom[[3]]')
            jofc_dist_mat_1_custom_dissimilar=getRdata('jofc.res.1.custom.ell1')
            %         jofc_dist_mat_2_custom_dissimilar=getRdata('jofc.res.2.custom')
            %         jofc_dist_mat_3_custom_dissimilar=getRdata('jofc.res.3.custom')
            matching_1_custom_dissimilar = YiCaoHungarian(jofc_dist_mat_1_custom_dissimilar);
            
            
            truematch_custom_dissimilar(q_i) = sum(matching_1_custom_dissimilar==1:n);
            
            pcjofc_sp(j+1,q_i) = pcjofc_sp(j+1,q_i)+truematch_shortest_path(q_i);
            pcjofc_diff(j+1,q_i) = pcjofc_diff(j+1,q_i)+truematch_diffusion(q_i);
            
            pcjofc_dice(j+1,q_i) =  pcjofc_dice(j+1,q_i) + truematch_dice(q_i);
            pcjofc_jaccard(j+1,q_i) =  pcjofc_jaccard(j+1,q_i) + truematch_jaccard(q_i);
            
            pcjofc_invlogweighted(j+1,q_i) =  pcjofc_invlogweighted(j+1,q_i) + truematch_invlogweighted(q_i);
            pcjofc_ell1(j+1,q_i) =  pcjofc_ell1(j+1,q_i) + truematch_custom_dissimilar(q_i);
            
            
            
        end
    end
end
pc=pc/numiter;
pcjofc_sp =pcjofc_sp/numiter;
pcjofc_diff =pcjofc_diff/numiter;
pcjofc_dice =   pcjofc_dice/numiter;
pcjofc_jaccard =  pcjofc_jaccard /numiter;

pcjofc_invlogweighted = pcjofc_invlogweighted/numiter;
pcjofc_ell1 = pcjofc_ell1/numiter;










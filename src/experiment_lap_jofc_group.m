function [pc,pcjofc] = experiment(n,maxm,numiter)

pc=zeros(1,maxm+1);
pcjofc = zeros(1,maxm+1);
truematch = 0;

num_w =1;
w_wals_vec = zeros(num_w,1);
w_vals_vec(1) = 0.95;

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
        bij=graphmatchHARDSEED(At,Bt,j);
        pc(j+1)=pc(j+1)+sum(bij(j+1:j+n)'==[j+1:j+n])/n;
        embed.dim =2;
        if (j<=2)
            continue
        else
            embed.dim = ceil( j/2) ;
            
        end
        
        putRdata('At',At)
        putRdata('Bt',Bt)
        insample_logic_vector = logical([ones(1,j) zeros(1,n) ]');
        putRdata('insample_logic_vec', insample_logic_vector)
        putRdata('embed.dim', embed.dim)
        
        evalR('insample_logic_vec<-as.vector(insample_logic_vec)')
        evalR('insample_logic_vec<-rep(insample_logic_vec,2)')
        putRdata('w_vals_vec',w_vals_vec)
        evalR('sink("debug.mat.txt")')
        %evalR('traceback()')
        evalR('print(insample_logic_vec)')
        evalR('sink()')
        %evalR('error.handle <- function()')
        evalR(['error.handle <- function(ex) {print(ex)}'])
        evalR('jofc.result <-try(JOFC.graph.custom.dist(G=At,Gp=Bt, in.sample.ind=insample_logic_vec,  d.dim=embed.dim,w.vals.vec=w_vals_vec,graph.is.directed=FALSE,vert_diss_measure="default"))')
        % evalR('jofc.result <- try(jofc(G=At,Gp=Bt, in.sample.ind=insample_logic_vec,  d.dim=embed.dim,w.vals.vec=0.99,graph.is.directed=FALSE, oos=TRUE,use.weighted.graph=FALSE))')
        % evalR('eval(parse(jofc.call))')
        %evalR('jofc.result<- tryCatch ({eval(parse(jofc.call))},error=error.handle)')
        %  traceback   })')
        
        evalR('skip.iter<-FALSE')
        evalR('	skip.iter <-inherits(jofc.result,"try-error")')
        if (getRdata('skip.iter'))
            j
            evalR('sink("debug.matlab.txt")')
            evalR('traceback()')
            evalR('print(jofc.result)')
            evalR('sink()')
            
            'Skipping iteration'
            continue
            
        end
        evalR('jofc.res.1<-jofc.result[[1]]')
        getRdata('jofc.res.1') ;
        evalR('M.result<-try(solveMarriage(jofc.res.1))')
        evalR('	skip.iter <-inherits(M.result,"try-error")')
        if (getRdata('skip.iter'))
            j
            'M.result error'
            evalR('sink("debug.matlab.2.txt")')
            evalR('traceback()')
            evalR('print(M.result)')
            evalR('sink()')
            
            
            'Skipping iteration'
            continue
            
        end
        evalR('NumofTruePairing<-present(M.result)')
        truematch = getRdata('NumofTruePairing');
        pcjofc(j+1) =  pcjofc(j+1) + truematch/n;
        
    end
end
pc=pc/numiter;
pcjofc =pcjofc /numiter;
%close
%plot([0:maxm],pc,'r*')
%plot([0:maxm],pcjofc,'r*')
end


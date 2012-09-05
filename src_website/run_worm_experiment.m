function [fc,sd_fc,random_chance,n_vals,num_iter]=run_worm_experiment(n_vals,num_iter,binarize,symmetrize)
%Function run_worm_experiment
%[fc,sd_fc,random_chance,n_vals,num_iter] =
% run_worm_experiment(n_vals,num_iter,binarize,symmetrize)
% input arguments
% n_vals: number of hardseeds
% num_iter: Number of MC replicates
% binarize: if nonzero, turn the edge weight (numerical) matrices to (binary) adjacency matrices
% symmetrize: if nonzero, make the matrices symmetric.
% return arguments
% fc : fraction of correct matches
% sd_fc : standard error of fc
% random_chance : expected number of correct matches under chance

load('elegansGraph.mat')

GE=Achem;
GF=Agap;

N_dims=size(GE)
N_init= N_dims(1)

% Binarize and Symmetrize the adjacency matrices

if (binarize)
    if (symmetrize)
        GE = GE+GE';
        GF=GF+GF';
        GE= (GE>0);
        
        GF= (GF>0);
    else
        GE= (GE>0);
        
        GF= (GF>0);
    end
else
    GE=(GE+(GE'))/2;
    GF=(GE+(GF'))/2;
    
end



GE= (GE+zeros(N_init));
GF= (GF+zeros(N_init));



N = N_init
GE=GE(1:N,1:N);
GF=GF(1:N,1:N);
rowsum_E=sum(GE,2);

unconnected_verts_G1=find(rowsum_E==0);
rowsum_F=sum(GF,2);
unconnected_verts_G2=find(rowsum_F==0);

n_vals = n_vals((n_vals<N));
%num_iter = 50;
corr_match=zeros(length(n_vals),num_iter);
for n_i = 1:length(n_vals)
    for i=1:num_iter
        i
        ordering=randperm(N);
        matching=ConVogHard_rQAP_order(GE,GF,n_vals(n_i),ordering,1);
        corr_match(n_i,i) =  sum(matching(n_vals(n_i)+1:N)==ordering(n_vals(n_i)+1:N));
    end
end


pc=mean(corr_match,2);
fc=pc./(N-n_vals');
sd_pc = std(corr_match,0,2);
sd_fc= sd_pc./(N-n_vals');

'Connectome Finished'

% Plot match ratio
random_chance= 1./(N-n_vals');
%plot(n_vals,fc,'r-','LineWidth',2)
errorbar(n_vals,fc,2*sd_fc/sqrt(num_iter),'r-')
hold on

plot(n_vals,random_chance,'k-.','LineWidth',2)

title('C. Elegans','FontSize',20)

xlabel('$m$','Interpreter','latex','FontSize',20)
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)
xlim([-5 max(n_vals)+5])






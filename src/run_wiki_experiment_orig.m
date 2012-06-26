
function [fc,sd_fc,random_chance,n_vals,num_iter]=run_wiki_experiment_orig(N,n_vals,num_iter)

%Function run_wiki_experiment
%[fc,sd_fc,random_chance,n_vals,num_iter]=run_wiki_experiment(N,n_vals,num_iter)
% input arguments 
% N: Since the matrix of wiki is 1382x1382 matrix which is too large for computational
% reason  we extract the first N vertices of the wiki graph (NxN submatrix)
% n_vals: number of hardseeds
% num_iter: Number of MC replicates
% return arguments 
% fc : fraction of correct matches
% sd_fc : standard error of fc
% random_chance : expected number of correct matches under chance
load('./data/wiki_adj.mat')

'Loaded wiki adjacency matrix'


N_dims=size(G_EN_Adj)
N_init= N_dims(1)

%N= 400
GE=G_EN_Adj(1:N,1:N);
GF=G_FR_Adj(1:N,1:N);
rowsum_E=sum(GE,2);
find(rowsum_E==0)
rowsum_F=sum(GF,2);
find(rowsum_F==0)

%n_vals=[0 1 5 10 20 50 100 200 300 350 400 450];
n_vals = n_vals(find(n_vals<N))
%num_iter = 50;
corr_match=zeros(length(n_vals),num_iter);
for n_i = 1:length(n_vals)
    for i=1:num_iter
    i
    ordering=randperm(N);
    matching=ConVogHard_rQAP_order(GE,GF,n_vals(n_i),ordering);
    corr_match(n_i,i) =  sum(matching(n_vals(n_i)+1:N)==ordering(n_vals(n_i)+1:N));
    end
end


pc=mean(corr_match,2)
fc=pc./(N-n_vals')
sd_pc = std(corr_match,0,2)
sd_fc= sd_pc./(N-n_vals')


'Wiki Finished'
figure
random_chance= 1./(N-n_vals');
%plot(n_vals,fc,'r-')
hold on

title('Wiki article matching-First 400 Verts')
errorbar(n_vals,fc,2*sd_fc/sqrt(num_iter),'r-')
xlabel('Number of Hard seeds')
ylabel('Fraction of Correct Matches')
xlim([-5 N+5])


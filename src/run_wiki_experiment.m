
function [fc,sd_fc,n_vals,num_iter]=run_wiki_experiment(N,n_vals,num_iter)

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

print('Loaded wiki adjacency matrix')

try
[status seed] = system('od /dev/urandom --read-bytes=4 -tu | awk ''{print $2}''');
seed=str2double(seed);
rng(seed);
catch
    rng shuffle
    'If running in parallel, parallel simulation might have the same random seed'
    'Check the seeds for uniqueness'    
end
currseed= rng();
save('random_rng.mat','currseed')

defaultStream = RandStream.getDefaultStream();
savedState = defaultStream.State;
save('random_rng_state.mat','savedState')

N_dims=size(G_EN_Adj)
N_all= N_dims(1)

%N= 400
GE=G_EN_Adj;
GF=G_FR_Adj;


GE(1:(N_all+1):(N_all*N_all))=0;
GF(1:(N_all+1):(N_all*N_all))=0;

rowsum_E=sum(GE,2);
find(rowsum_E==0)
rowsum_F=sum(GF,2);
find(rowsum_F==0)

%n_vals=[0 1 5 10 20 50 100 200 300 350 400 450];

%num_iter = 50;
corr_match=zeros(length(n_vals),num_iter);
for i=1:num_iter
    i
    ordering=randperm(N_all);
    for n_i = 1:length(n_vals)
        totv = (N+n_vals(n_i));
        %Do this ugly very ugly hack to keep the same  test vertices
        %And have adjacency matrix for the seed vertices on the upper left of matrices
        seed_plus_test_vertex_indices = ordering([(N+1):totv 1:N]);
        
        GE=G_EN_Adj(seed_plus_test_vertex_indices,seed_plus_test_vertex_indices);
        GF=G_FR_Adj(seed_plus_test_vertex_indices,seed_plus_test_vertex_indices);
        %After the preciding two lines, the original indices of the
        %vertices are irrelevant neither ConVogHard_rQAP nor matching
        %treats the vertices as numbered from 1 to N+n_vals(i) where the
        %first n_vals(i) vertices are hard seeds
        matching=ConVogHard_rQAP(GE,GF,n_vals(n_i));
        corr_match(n_i,i) =  sum(matching((n_vals(n_i)+1):totv)==((n_vals(n_i)+1):totv));
    end
end


pc=mean(corr_match,2)
fc=pc./N
sd_pc = std(corr_match,0,2)
sd_fc= sd_pc./N


'Wiki Finished'
figure
%random_chance= 1./(N-n_vals');
%plot(n_vals,fc,'r-')
hold on

title('Wiki article matching-500 randomly sampled vertices plus seeds')
errorbar(n_vals,fc,2*sd_fc/sqrt(num_iter),'r-')
xlabel('Number of Hard seeds')
ylabel('Fraction of Correct Matches')
xlim([-5 N+5])
end


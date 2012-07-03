
load('./data/wiki_adj.mat')

'Loaded wiki adjacency matrix'


N_dims=size(G_EN_Adj)
N_all= N_dims(1)

%N= 400



n_vals=[0 1  20 200   600 ];

num_iter = 20;
corr_match=zeros(length(n_vals),num_iter);

N_vals=[50 60 70 80 90 100 125 150 175 200 225 250 275 300 350 400 450 500 550];
ex_time=zeros(length(N_vals),num_iter);
ex_time_seed =zeros(length(n_vals),length(N_vals));
for N_it=1:length(N_vals)
    N=N_vals(N_it);
    
    for i=1:num_iter
    GE=[];
    GF=[];
        tic;
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
    ex_time(N_it,i)=toc
    mean(ex_time,2)
    
    end
end
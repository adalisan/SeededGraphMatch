
load('./data/wiki_adj.mat')

print('Loaded wiki adjacency matrix')


N_dims=size(G_EN_Adj)
N_all= N_dims(1)

%N= 400
GE=G_EN_Adj;
GF=G_FR_Adj;
rowsum_E=sum(GE,2);
find(rowsum_E==0)
rowsum_F=sum(GF,2);
find(rowsum_F==0)

%n_vals=[0 1 5 10 20 50 100 200 300 350 400 450];

%num_iter = 50;
corr_match=zeros(length(n_vals),num_iter);

N_vals=[50 100 300 400 500 550];
ex_time=zeros(length(N_vals),num_iter);
for N_it=1:length(N_vals)
    N=N_vals(N_it)
    
for i=1:num_iter
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
        matching=ConVogHard_rQAP(GE,GF,n_vals(n_i));
        corr_match(n_i,i) =  sum(matching((n_vals(n_i)+1):totv)==((n_vals(n_i)+1):totv));
    end
    ex_time(N_it,i)=toc;
end
end
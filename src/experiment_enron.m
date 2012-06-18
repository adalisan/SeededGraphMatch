

evalR('setwd("C:/Users/Sancar/Documents/projects/DataFusion/data")')
evalR('load("AAA-187As-184x184.Rbin")')
evalR('sink("test_R_int.txt")')
evalR('print(str(AAA[[1]]))')
evalR('sink()')

evalR('GE=AAA[[131]]')
evalR('GF=AAA[[132]]')


GE=getRdata('GE');
GF=getRdata('GF');

N_dims=size(GE)
N_init= N_dims(1)

N= N_init%floor(N_init/4)
GE=GE(1:N,1:N);
GF=GF(1:N,1:N);
rowsum_E=sum(GE,2);
find(rowsum_E==0)
rowsum_F=sum(GF,2);
find(rowsum_F==0)
%n_vals=[1 5 10 20 50 100];
n_vals = [0 1 5 10 20 50 60 90 100 140]
%n_vals=[0 90  ];

num_iter = 60;
corr_match=zeros(length(n_vals),num_iter);
for n_i = 1:length(n_vals)
    for i=1:num_iter
    i
    ordering=randperm(N);
    matching=ConVogHard_rQAP_order(GE,GF,n_vals(n_i),ordering);
    corr_match(n_i,i) =  sum(matching(n_vals(n_i)+1:N)==ordering(n_vals(n_i)+1:N))
    end
end

pc=mean(corr_match,2)
sd_pc = std(corr_match,0,2)

fc= pc./(N-n_vals')
sd_fc= sd_pc./(N-n_vals')

'Enron Finished'
plot(n_vals,fc)


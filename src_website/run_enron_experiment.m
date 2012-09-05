
function [fc,sd_fc,random_chance,corr_match,pc,sd_pc, ...
    mismatched_verts_G1,mismatched_verts_G2,iter_in_test]= ...
    run_enron_experiment(n_vals,num_iter,time_stamp_G1,time_stamp_G2)
%Function run_enron_experiment
%[fc,sd_fc,random_chance,n_vals,num_iter]=run_enron_experiment(n_vals,num_iter,time_stamp_G1,time_stamp_G2)
% input arguments
% n_vals: number of hardseeds
% num_iter: Number of MC replicates
% time_stamp_G1,time_stamp_G2 : which graphs in the time series of graphs
% should we match
% return arguments
% fc : average fraction of correct matches
% sd_fc : sample standard deviation  of fc
% random_chance : expected fraction of  number of correct matches under chance
% corr_match: Number of correct matches for each number of hard seeds in each iteration
% pc: average number of correct matches
% sd_pc: sample standard deviation  of pc
% mismatched_verts_G1: indices of vertices that are mismatched for num_seeds_look_at=5;
% the collection of such vertices for each iteration are separated by NaNs
% mismatched_verts_G2: should be the same as mismatched_verts_G1
if (~exist('AAA','var'))
    load('enron.mat')
end




GE=AAA(:,:,time_stamp_G1);
GF=AAA(:,:,time_stamp_G2);

N_dims=size(GE);
N_init= N_dims(1);

N= N_init;
GE=GE(1:N,1:N);
GF=GF(1:N,1:N);
rowsum_E=sum(GE,2);
find(rowsum_E==0);
rowsum_F=sum(GF,2);
find(rowsum_F==0);

mismatched_verts_G1=[];
mismatched_verts_G2=[];

corr_match=zeros(length(n_vals),num_iter);

iter_in_test=zeros(N,length(n_vals));

% If there is a  list of "anomalous" vertices that are known to be not very
% well across time_stamps, and want to consider them as test_vertices for every replicate,
% put them in N_2_neigh_anomaly
N_2_neigh_anomaly = [];
for n_i = 1:length(n_vals)
    for i=1:num_iter
        i
        ordering = randperm(N);
        %Modify the ordering to keep  N_2_neigh_anomaly vertices among test
        %vertices to be matched
        if (length(N_2_neigh_anomaly)>0)
            keep_in_test = randsample( N_2_neigh_anomaly,min(N-n_vals(n_i),length(N_2_neigh_anomaly)));
            for l=1:length(keep_in_test)
                change_index= ordering==keep_in_test(l);
                tmp = ordering(n_vals(n_i)+l);
                ordering(n_vals(n_i)+l)= keep_in_test(l);
                ordering(change_index)= tmp;
            end
        end
        
        ordering(1:n_vals(n_i)) = sort(ordering(1:n_vals(n_i)));
        ordering(n_vals(n_i)+1:N) = sort(ordering(n_vals(n_i)+1:N));
        iter_in_test(ordering(n_vals(n_i)+1:N),n_i) = iter_in_test(ordering(n_vals(n_i)+1:N),n_i) + 1;
        matching=ConVogHard_rQAP_order(GE,GF,n_vals(n_i),ordering,1);
        corr_match(n_i,i) =  sum(matching(n_vals(n_i)+1:N)==ordering(n_vals(n_i)+1:N));
        
    end
end

pc=mean(corr_match,2);
sd_pc = std(corr_match,0,2);

fc= pc./(N-n_vals');
sd_fc= sd_pc./(N-n_vals');

'Enron Finished'

random_chance= 1./(N-n_vals');


title('Enron graph matching, graphs at t=130, 131, 132')
%legend('chance','130&131','131&132','130&132')
xlim([-5 145]);




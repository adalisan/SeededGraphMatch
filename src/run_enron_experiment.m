
function [fc,sd_fc,random_chance,corr_match,pc,sd_pc,mismatched_verts_G1,mismatched_verts_G2]=run_enron_experiment(n_vals,num_iter,time_stamp_G1,time_stamp_G2)
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
load('./data/enron.mat')



GE=AAA(:,:,time_stamp_G1);
GF=AAA(:,:,time_stamp_G2);

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
%n_vals = [0 1 5 10 20 50 60 90 100 140]
%n_vals=[0 90  ];
mismatched_verts_G1=[];
mismatched_verts_G2=[];

num_seeds_look_at=5;
%num_iter = 60;
seeds=zeros(num_iter,num_seeds_look_at);

corr_match=zeros(length(n_vals),num_iter);
%keep_in_test=[90 59 83 127];
for n_i = 1:length(n_vals)
    for i=1:num_iter
    i
    ordering = randperm(N);
%     for l=1:length(keep_in_test)
%        change_index=find(ordering==keep_in_test(l));
%        tmp = ordering(n_vals(n_i)+l);
%        ordering(n_vals(n_i)+l)= keep_in_test(l);
%        ordering(change_index)= tmp;
%     end
%     
    ordering(1:n_vals(n_i)) = sort(ordering(1:n_vals(n_i)));
    ordering(n_vals(n_i)+1:N) = sort(ordering(n_vals(n_i)+1:N));
    
    matching=ConVogHard_rQAP_order(GE,GF,n_vals(n_i),ordering);
    corr_match(n_i,i) =  sum(matching(n_vals(n_i)+1:N)==ordering(n_vals(n_i)+1:N));
    if (n_vals(n_i)==num_seeds_look_at)
        mismatched_ind=find(matching(n_vals(n_i)+1:N)~=ordering(n_vals(n_i)+1:N));
        mismatched_verts_G1=[mismatched_verts_G1 NaN matching(n_vals(n_i)+mismatched_ind)]
      
        mismatched_verts_G2=[mismatched_verts_G2 NaN ordering(n_vals(n_i)+mismatched_ind)]
        seeds(i,:) = ordering(1:n_vals(n_i));
    end
    end
end

pc=mean(corr_match,2)
sd_pc = std(corr_match,0,2)

fc= pc./(N-n_vals')
sd_fc= sd_pc./(N-n_vals')

'Enron Finished'

random_chance= 1./(N-n_vals');

figure


if time_stamp_G1==130 && time_stamp_G2==131
corr_match_run0=corr_match;
pc_run0 = pc;
sd_pc_run0 = sd_pc;
fc_run0 = fc;
sd_fc_run0 = sd_fc;
mismatched_verts_G1_run0 = mismatched_verts_G1;
mismatched_verts_G2_run0 = mismatched_verts_G2;
hold on
plot(n_vals,random_chance,'k-.')

errorbar(n_vals,fc_run0,2*sd_fc_run0/sqrt(num_iter),'r-')

elseif time_stamp_G1==131 && time_stamp_G2==132

corr_match_run1=corr_match;
pc_run1 = pc;
sd_pc_run1 = sd_pc;
fc_run1 = fc;
sd_fc_run1 = sd_fc;
mismatched_verts_G1_run1 = mismatched_verts_G1;
mismatched_verts_G2_run1 = mismatched_verts_G2;


hold on 
plot(n_vals,random_chance,'k-.')
errorbar(n_vals,fc_run1,2*sd_fc_run1/sqrt(num_iter),'b-')

elseif  time_stamp_G1==130 && time_stamp_G2==132
corr_match_run2=corr_match;
pc_run2 = pc;
sd_pc_run2 = sd_pc;
fc_run2 = fc;
sd_fc_run2 = sd_fc;
mismatched_verts_G1_run2 = mismatched_verts_G1;
mismatched_verts_G2_run2 = mismatched_verts_G2;

hold on 
plot(n_vals,random_chance,'k-.')
errorbar(n_vals,fc_run2,2*sd_fc_run2/sqrt(num_iter),'g-')

else
    %Not one of our test cases so no need for the rest of the function
    %plotting and so on
return
end



title('Enron graph matching, graphs at t=130, 131, 132')
%legend('chance','130&131','131&132','130&132')
xlim([-5 145]);




function [fc,sd_fc,random_chance,n_vals,num_iter]=run_general_experiment(GE,GF,n_vals,num_iter)
%Function run_general_experiment
%[fc,sd_fc,random_chance,n_vals,num_iter]=run_worm_experiment(n_vals,num_iter)
% input arguments 
% GE : First graph to be matched 
% GF : Second graph to be matched
% n_vals: number of hardseeds
% num_iter: Number of MC replicates
% return arguments 
% fc : fraction of correct matches
% sd_fc : standard error of fc
% random_chance : expected number of correct matches under chance


N_dims=size(GE)
N_init= N_dims(1)

naive_symmetrization = 1;

if (naive_symmetrization)
%GE= GE-triu(GE);
GE= (GE>0);

%GE=GE+GE';
%GF= GF-triu(GF);
GF= (GF>0); 
%GF=GF+GF';
else

%     sym_GE= zeros(N_init);
%     sym_GF= zeros(N_init);
%     sym_GE= sym_GE+triu(GE); 
%     sym_GF= sym_GF+triu(GF);
%     
%     sym_GE= sym_GE+sym_GE'; 
%     sym_GF= sym_GF+sym_GF';
%     
%     sym_GE= sym_GE+tril(GE);
%     sym_GF= sym_GF+tril(GF); 
  GE = GE+GE';
  GF=GF+GF';
  GE= (GE>0);

  GF= (GF>0); 
  
  
end

GE= (GE+zeros(N_init));
GF= (GF+zeros(N_init));



N = N_init
GE=GE(1:N,1:N);
GF=GF(1:N,1:N);
rowsum_E=sum(GE,2);
%There's no hope of matching the following vertices correctly
unconnected_verts_G1=find(rowsum_E==0)
rowsum_F=sum(GF,2);
unconnected_verts_G2=find(rowsum_F==0)

%n_vals=[0 1 5 10 20 50 75 100 150 200 300 350 400 450];
n_vals = n_vals(find(n_vals<N))
%num_iter = 50;
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
fc=pc./(N-n_vals')
sd_pc = std(corr_match,0,2)
sd_fc= sd_pc./(N-n_vals')


'Simulation Finished'
random_chance= 1./(N-n_vals');
plot(n_vals,fc,'r-')
hold on

plot(n_vals,random_chance,'b-.')



errorbar(n_vals,fc,2*sd_fc/sqrt(num_iter),'r-')
xlabel('Number of Hard seeds')
ylabel('Fraction of Correct Matches')
xlim([-5 N+5])





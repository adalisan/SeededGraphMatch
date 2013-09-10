
function [fc,sd_fc,fc_noseed,sd_fc_noseed,random_chance,n_vals,num_iter]=run_cnet_experiment_orig(n_vals,num_iter,plot_fig,subset)

%Function run_cnet_experiment_orig
%Carey's version of experiments
%[fc,sd_fc,random_chance,n_vals,num_iter]=run_cnet_experiment_orig(n_vals,num_iter)
% input arguments 
% n_vals: number of hardseeds
% num_iter: Number of MC replicates
%
% return arguments 
% fc : fraction of correct matches
% sd_fc : standard error of fc
% random_chance : expected number of correct matches under chance
load('./data/cnet_Ajt.mat')


'Loaded cnet adjacency matrix'


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

defaultStream = RandStream.getGlobalStream ;
savedState = defaultStream.State;
save('random_rng_state.mat','savedState')



GE=uint16(Ajt1);
GF=uint16(Ajt2);
N_dims=size(GE);
N_all=N_dims(1);
GE(1:(N_all+1):(N_all*N_all))=0;
GF(1:(N_all+1):(N_all*N_all))=0;

row_E = (sum(GE,2)==0 ) ;
col_E = (sum(GE,1)==0 ) ;
row_F = (sum(GF,2)==0 ) ;
col_F = (sum(GF,1)==0 ) ;

unconnected_verts_G1=(row_E &col_E');

unconnected_verts_G2= (row_F &col_F');
unconnected_verts = unconnected_verts_G1 | unconnected_verts_G2;
GE=GE(~unconnected_verts,~unconnected_verts);
GF=GF(~unconnected_verts,~unconnected_verts);
GE=double(GE);
GF=double(GF);

N_dims=size(GE);

N_all= N_dims(1);

if (subset<N_all)
    subset_ind = randperm(N_all,subset);
    GE=GE(subset_ind,subset_ind);
    GF=GF(subset_ind,subset_ind);
    
end
N_all=subset;

%n_vals=[0 1 5 10 20 50 100 200 300 350 400 450];
n_vals = n_vals((n_vals<N_all));
%num_iter = 50;
corr_match=zeros(length(n_vals),num_iter);
corr_match_no_seed=zeros(length(n_vals),num_iter);
 for i=1:num_iter
    i
    for n_i = 1:length(n_vals)
    	n_i
        % A random ordering for such that the first n_vals(n_i) are seeds
    ordering=randperm(N_all);
        
        matching=ConVogHard_rQAP_order(GE,GF,n_vals(n_i),ordering,1);
    test_v=(n_vals(n_i)+1):N_all;
        %Number of correct matches among test vertices
    corr_match(n_i,i) =  sum(matching(test_v)==ordering(test_v));
    test_v_ord=sort(ordering(test_v));
        
        %Match with no seeds
    matching_no_seed=ConVogHard_rQAP(GE(test_v_ord,test_v_ord),GF(test_v_ord,test_v_ord), ...
        0);
    corr_match_no_seed(n_i,i) =  sum(matching_no_seed==(1:(N_all-n_vals(n_i))));
    save('cnet-all-orig.mat')
    end
end


pc=mean(corr_match,2);
fc=pc./(N_all-n_vals');
sd_pc = std(corr_match,0,2);
sd_fc= sd_pc./(N_all-n_vals');

pc_noseed=mean(corr_match_no_seed,2);
fc_noseed=pc_noseed./(N_all-n_vals');
sd_pc_noseed = std(corr_match_no_seed,0,2);
sd_fc_noseed= sd_pc_noseed./(N_all-n_vals');

random_chance= 1./(N_all-n_vals');
'charitynet Finished'
if plot_fig

figure


hold on

title('Charitynet','FontSize',20)
errorbar(n_vals,fc,2*sd_fc/sqrt(num_iter),'r-','LineWidth',2)
hold on
errorbar(n_vals,fc_noseed,2*sd_fc_noseed/sqrt(num_iter),'b-','LineWidth',2)
xlabel('$m$','Interpreter','latex','FontSize',20)
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)
plot(n_vals,1./(N_all-n_vals),'k-.','LineWidth',2)


xlim([-5 max(n_vals)+5])
legend('seeded','no seeds','chance')
end
end


function [fc,sd_fc,fc_noseed,sd_fc_noseed,random_chance,n_vals,num_iter,fc_unwt,sd_fc_unwt] = ...
    run_worm_experiment(n_vals,num_iter,binarize,symmetrize,compare_wt_unwt,compare_unseed)
%Function run_worm_experiment
% [fc,sd_fc,fc_noseed,sd_fc_noseed,random_chance,n_vals,num_iter,fc_unwt,sd_fc_unwt=
% run_worm_experiment(n_vals,num_iter,binarize,symmetrize,compare_wt_unwt,compare_unseed)
% input arguments
% n_vals: number of hardseeds
% num_iter: Number of MC replicates
% binarize: if nonzero, turn the edge weight (numerical) matrices to (binary) adjacency matrices
% symmetrize: if nonzero, make the matrices symmetric.
%compare_wt_unwt if nonzero runs sim for both weighted and unweighted graphs
% compare_unseed: if nonzero, runs sim for the same test vertices with no
% seeds
% return arguments
% fc : fraction of correct matches
% sd_fc : standard error of fc
% random_chance : expected number of correct matches under chance

load('./data/elegansGraph.mat')


plot_graph =  0;



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



GE=Achem;
GF=Agap;

load('./data/elegansGraph_male.mat')

N_dims=size(GE);
N_init= N_dims(1)

GE_male= Achem;
GF_male= Agap;

N_dims_m =size(GE_male);
N_sec = N_dims_m[1];
matched_v =(ismember(Neuron_ordered,Achem_names));
common_v = sum(matched_v);
GE_m = zeros(max(N_init,N_sec));
GE_m(1:common_v,1:common_v) = GE_male;


N_dims_m =size(GF_male);
N_sec = N_dims_m[1];
matched_v =(ismember(Neuron_ordered,Agap_names));
common_v = sum(matched_v);
GF_m = zeros(max(N_init,N_sec));
GF_m(1:common_v,1:common_v) = GF_male;






% Binarize and Symmetrize the adjacency matrices
if  (compare_wt_unwt)
       if (symmetrize)
      GE=(GE+(GE'))/2;
       GF=(GF+(GF'))/2;
            GE_m=(GE_m+(GE_m'))/2;
       GF_m=(v+(GF_m'))/2;
       end
    
    GE_unwt =     double(GE>0);    
    GF_unwt =     double(GF>0);
    
    GE_m_unwt =     double(GE_m>0);    
    GF_m_unwt =     double(GF_m>0);
    
    
    %GE = GE/ norm(GE,'fro');
    %GF = GF/ norm(GF,'fro');
    
elseif (binarize)
    if (symmetrize)
        GE = GE+GE';
        GF=GF+GF';
        
           GE= double(GE>0);
        GF= double(GF>0);
        
          GE_m = GE_m+GE_m';
        GF_m=GF_m+GF_m';
        
     
        
          GE_m= double(GE_m>0);
        GF_m= double(GF_m>0);
    else
        GE= double(GE>0);
        GF= double(GF>0);
        
             GE_m= double(GE_m>0);
        GF_m= double(GF_m>0);
    end
else
    if (symmetrize)
    GE=(GE+(GE'))/2;
    GF=(GE+(GF'))/2;
    
    GE_m=(GE_m+(GE_m'))/2;
    GF_m=(GF_m+(GF_m'))/2;
    
    end
    GE = GE/ norm(GE,'fro');
    GF = GF/ norm(GF,'fro');
    
end


GE= full(GE);
GF= full(GF);

GE_unwt= full(GE_unwt);
GF_unwt= full(GF_unwt);

is_sym_GE =   sum((GE-GE')~=0)
is_sym_GF =   sum((GF-GF')~=0)

N = N_init;
GE=GE(1:N,1:N);
GF=GF(1:N,1:N);

%There's no hope of matching the following vertices more correctly than chance
row_E = (sum(GE,2)==0 ) ;
col_E = (sum(GE,1)==0 ) ;
row_F = (sum(GF,2)==0 ) ;
col_F = (sum(GF,1)==0 ) ;

unconnected_verts_G1=(row_E &col_E');

unconnected_verts_G2= (row_F &col_F');
unconnected_verts = unconnected_verts_G1 | unconnected_verts_G2;
GE=GE(~unconnected_verts,~unconnected_verts);
GF=GF(~unconnected_verts,~unconnected_verts);
GE_unwt=GE_unwt(~unconnected_verts,~unconnected_verts);
GF_unwt=GF_unwt(~unconnected_verts,~unconnected_verts);

N_dims=size(GE);
N = N_dims(1)
size(GE_unwt)
size(GF_unwt)

%n_vals=[0 1 5 10 20 50 75 100 150 200 300 350 400 450];
n_vals = n_vals((n_vals<N));
%num_iter = 50;
corr_match = zeros(length(n_vals),num_iter);
corr_match_unwt = zeros(length(n_vals),num_iter);
corr_match_no_seed= zeros(length(n_vals),num_iter);
for n_i = 1:length(n_vals)
    for i=1:num_iter
        'n_i  n_vals(n_i)  i'
        [n_i  n_vals(n_i)  i]
        ordering=randperm(N);
        matching=ConVogHard_rQAP_order(GE,GF,n_vals(n_i),ordering,1);
        corr_match(n_i,i) =  sum(matching(n_vals(n_i)+1:N)==ordering(n_vals(n_i)+1:N));
        if  (compare_wt_unwt)
            matching_unwt=ConVogHard_rQAP_order(GE_unwt,GF_unwt,n_vals(n_i),ordering,1);
            corr_match_unwt(n_i,i) =  sum(matching_unwt(n_vals(n_i)+1:N)==ordering(n_vals(n_i)+1:N));
        end
        test_v=(n_vals(n_i)+1):N;
        % For n_vals(1)=0 , we match with no seeds, so the all of the vertices
        % are matched and for different  i, we would get the same number of
        % matchings.  so first entry of sd_pc_noseed  
        % For other n_i, different i correspond to different seed 
        % selections, so different  corr_match_no_seed are expected for
        % different i.
        if  (compare_unseed)
            test_v_ord=sort(ordering(test_v));
             num_match= length(test_v_ord);
            %Match with no seeds
            matching_unseed=ConVogHard_rQAP(GE(test_v_ord,test_v_ord),GF(test_v_ord,test_v_ord), ...
                0);
            corr_match_no_seed(n_i,i) =  sum(matching_unseed==(1:num_match));
            
            
        end
        
    end
end


pc=mean(corr_match,2);
fc=pc./(N-n_vals');
sd_pc = std(corr_match,0,2);
sd_fc= sd_pc./(N-n_vals');

pc_unwt=mean(corr_match_unwt,2);
fc_unwt=pc_unwt./(N-n_vals');
sd_pc_unwt = std(corr_match_unwt,0,2);
sd_fc_unwt= sd_pc_unwt./(N-n_vals');
pc_noseed=mean(corr_match_no_seed,2);
fc_noseed=pc_noseed./(N-n_vals');
sd_pc_noseed = std(corr_match_no_seed,0,2);
sd_fc_noseed= sd_pc_noseed./(N-n_vals');

random_chance= 1./(N-n_vals');

'Connectome Finished'
if (plot_graph)
    % Plot match ratio
    
    %plot(n_vals,fc,'r-','LineWidth',2)
    errorbar(n_vals,fc,2*sd_fc/sqrt(num_iter),'r-')
    hold on
    
    plot(n_vals,random_chance,'k-.','LineWidth',2)
    
    
    %title('C. Elegans Connectome- Achem vs Agap (Simple Graph) ')
    title('C. Elegans','FontSize',20)
    
    xlabel('$m$','Interpreter','latex','FontSize',20)
    ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)
    xlim([-5 max(n_vals)+5])
    %legend('Unweighted','Weighted','Chance')
end




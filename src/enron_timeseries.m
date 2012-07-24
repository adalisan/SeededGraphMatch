function [fc_agg,sd_fc_agg,corr_match_agg,pc_agg,sd_pc_agg, ...
    mismatched_verts_G1_runs,mismatched_verts_G2_runs, iter_in_test_agg] = ...
    enron_timeseries( n_vals,num_iter,timestamps )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

total_run= (length(timestamps)-1)*length(timestamps)/2;
fc_agg=cell(total_run,1);
sd_fc_agg=cell(total_run,1);
corr_match_agg=cell(total_run,1);
pc_agg=cell(total_run,1);
sd_pc_agg=cell(total_run,1);
iter_in_test_agg = cell(total_run,1);
mismatched_verts_G1_runs=cell(total_run,1);
mismatched_verts_G2_runs=cell(total_run,1);
linetypes =['r-','g-' ,'b-','y-','l-'];
run_count=0;

for time_i=1:(length(timestamps)-1)
    for time_j=(time_i+1):length(timestamps)
        run_count=run_count+1
        [fc,sd_fc,random_chance,corr_match,pc,sd_pc, ...
            mismatched_verts_G1,mismatched_verts_G2,iter_in_test]= ...
            run_enron_experiment(n_vals,num_iter,timestamps(time_i),timestamps(time_j));
        
        mismatched_verts_G1_runs{run_count} = mismatched_verts_G1;
        mismatched_verts_G2_runs{run_count} = mismatched_verts_G2;
        
        fc_agg{run_count}= fc;
        sd_fc_agg{run_count} = sd_fc;
        corr_match_agg{run_count} = corr_match;
        pc_agg{run_count} = pc;
        sd_pc_agg{run_count} = sd_pc;
        iter_in_test_agg{run_count}=iter_in_test;
        plot(n_vals,random_chance,'k-.')
        errorbar(n_vals,fc,2*sd_fc/sqrt(num_iter),linetypes(run_count))
        
    end
end

end


[fc,sd_fc,fc_seed,sd_fc_seed,random_chance,n_vals,num_iter]=run_wiki_experiment_orig([0:2:20 20:5:95  100:10:240 250:50:1300] ,1,0);
%[fc,sd_fc,fc_seed,sd_fc_seed,random_chance,n_vals,num_iter]=run_wiki_experiment_orig([1300] ,1,0);

fname = strcat('./cache/wiki_', datestr(clock) ,int2str(uint32(randi(1E6,1))),'.mat')
save(fname)

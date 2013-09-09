[fc,sd_fc,n_vals,num_iter] = run_cnet_experiment(350, [4:4:20 30:10:95  100:20:280] ,10);
%[fc,sd_fc,fc_seed,sd_fc_seed,random_chance,n_vals,num_iter]=run_wiki_experiment_orig([1300] ,1,0);

fname = strcat('./cache/cnet/cnet_', datestr(clock) ,int2str(uint32(randi(1E6,1))),'.mat')
save(fname)

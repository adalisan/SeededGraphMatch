

%Run sim for directed graph (unsymmetrized adj or weight matrices)
[fc_dir,sd_fc_dir,fc_unseed_dir,sd_fc_unseed_dir,~,~,~,fc_unwt_dir,sd_fc_unwt_dir]=run_worm_experiment([ 0 5 245] ,1,0,0,1,1);



%Run sim for undirected graph (symmetrized adj or weight matrices)
[fc,sd_fc,fc_seed,sd_fc_seed,random_chance,n_vals,num_iter,fc_unwt,sd_fc_unwt]=run_worm_experiment([ 0 5 245],1,0,1,1,1);
fname = strcat('./cache/worm/worm_', datestr(clock) ,int2str(uint32(randi(1E3,1))),'.mat')


fname = strcat('./cache/worm/worm_pre_', datestr(clock) ,int2str(uint32(randi(1E3,1))),'.mat')
save(fname)
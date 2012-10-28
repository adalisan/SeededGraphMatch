

[fc,sd_fc,fc_seed,sd_fc_seed,random_chance,n_vals,num_iter,fc_unwt,sd_fc_unwt]=run_worm_experiment([0:2:20 20:5:95  100: 10:200 ] ,1,0,0,1);

fname = strcat('wiki_', datestr(clock) ,randi(1E6,1), '.mat')
save(fname)
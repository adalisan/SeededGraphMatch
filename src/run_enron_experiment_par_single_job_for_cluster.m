
%Run sim for directed graph (unsymmetrized adj or weight matrices)

n_vals_enron =[ 0 1 2 3 4 ...
     5     8     9    10    11    12    13 ...
    14    15    16    17    18    19    20 ...
    21    24    27    30    33    36    39 ...
    42    45    48    50    60    70    80 ...
    90   100   110   120   130   140];
[fc_dir,sd_fc_dir,rand_c,corr_m,pc_dir,sd_pc_dir]              =run_enron_experiment(n_vals_enron ,1,130,131,0);
[fc_undir,sd_fc_undir,rand_c,corr_m_undir,pc_undir,sd_pc_undir]=run_enron_experiment(n_vals_enron,1,130,131,1);

fname = strcat('./cache/enron/enron_', datestr(clock) ,int2str(uint32(randi(1E3,1))),'.mat')
save(fname)
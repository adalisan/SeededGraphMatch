mismatched_verts_G2_run0=mismatched_verts_G2_runs{1}{1}{2};
mismatched_verts_G2_run1=mismatched_verts_G2_runs{1}{2};
mismatched_verts_G2_run2=mismatched_verts_G2_runs{2};
mismatched_verts_G1_run0=mismatched_verts_G1_runs{1}{1}{2};
mismatched_verts_G1_run1=mismatched_verts_G1_runs{1}{2};
mismatched_verts_G1_run2=mismatched_verts_G1_runs{2};

mismatched_verts_G2_run0=mismatched_verts_G2_run0(find(~isnan(mismatched_verts_G2_run0)));
mismatched_verts_G2_run1=mismatched_verts_G2_run1(find(~isnan(mismatched_verts_G2_run1)));
mismatched_verts_G2_run2=mismatched_verts_G2_run2(find(~isnan(mismatched_verts_G2_run2)));
mismatched_verts_G1_run0=mismatched_verts_G1_run0(find(~isnan(mismatched_verts_G1_run0)));
mismatched_verts_G1_run1=mismatched_verts_G1_run1(find(~isnan(mismatched_verts_G1_run1)));
mismatched_verts_G1_run2=mismatched_verts_G1_run2(find(~isnan(mismatched_verts_G1_run2)));

[counts_0,values_0] = hist(mismatched_verts_G2_run0(:),1:184);
[counts_1,values_1] = hist(mismatched_verts_G2_run1(:),1:184);
[counts_2,values_2] = hist(mismatched_verts_G2_run2(:),1:184);
th_1=25;
th_2=45;
mismatch_detect= counts_0<th_1 &counts_1>th_2
anomaly_gr_truth=zeros(1,184);
anomaly_gr_truth(N_2_neigh_anomaly)=1;
cont_table=zeros(2,2)
cont_table(1,1)=sum(mismatch_detect&anomaly_gr_truth);
cont_table(1,2)=sum(mismatch_detect&~anomaly_gr_truth);
cont_table(2,1)=sum(~mismatch_detect&anomaly_gr_truth);
cont_table(2,2)=sum(~mismatch_detect&~anomaly_gr_truth);
[cont_table_2,chi_teststat,p_val]=crosstab(mismatch_detect,anomaly_gr_truth);


v_2=mismatched_verts_G1_runs{2};
v_2=v_2(~isnan(v_2));
v_2_counts= hist(v_2,1:184);
v_2_freq = v_2_counts'./sum(iter_in_test_agg{1},2);
[v_2_freq_sorted,v_2_detect]=sort(v_2_freq,'descend');

v_1=mismatched_verts_G1_runs{1};
v_1=v_1(~isnan(v_1));
v_1_counts= hist(v_1,1:184);
v_1_freq = v_1_counts'./sum(iter_in_test_agg{1},2);
[v_1_freq_sorted,v_1_detect]=sort(v_1_freq,'descend');



v_0 = N_2_neigh_anomaly';



num_detected=zeros(1,184) ;
num_corr_detect=0;
j=1;
for i=1:184
    if (sum(v_1_detect(1:23)==v_2_detect(i)))
     continue    
    end
    
    if (sum(v_0==v_2_detect(i))) 
        num_corr_detect= num_corr_detect+1;
    end
    num_detected(j)=num_corr_detect;
    j=j+1;
end

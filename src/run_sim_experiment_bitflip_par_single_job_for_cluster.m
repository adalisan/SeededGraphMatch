%Randomization from system
function [] = run_sim_bitflip_single_job( job_i)

[status seed] = system('od /dev/urandom --read-bytes=4 -tu | awk ''{print $2}''');
seed=str2double(seed);
rng(seed);




q= [0:0.05:0.5 ];
%q=[0 0.1 0.3 0.45 0.5];
q_len = length(q);



q_int=find(q==0.3);



N=300;
numiter=2;
slp_iter =-1;
n_vals=[0:39 40:2:98 100:5:195 200:5:275  ];
%n_vals=[0:1:5 6:2:20 20:5:35 ];
%n_vals= [0:20 22 24 26];
%n_vals=0:16
%n_vals=[ 0 1 5 15 17:2:25 30 35 40 60 80 100 300 450];
n_vals=n_vals(n_vals<N);

n_len = length(n_vals)
P_jv_found=[]
P_l1_found=[]
P_l2_found=[]

P_jv_cell=cell(n_len,numiter);

P_jv_pr_cell=cell(n_len,numiter);

P_l2_cell=cell(n_len,numiter);
P_l2_pr_cell=cell(n_len,numiter);
P_l1_cell=cell(n_len,numiter);
P_l1_pr_cell=cell(n_len,numiter);

running_time_FAQ=zeros(numiter,n_len,q_len);
running_time_SLP=zeros(numiter,n_len,q_len);

corr_match=zeros(n_len,numiter,q_len);

corr_match_ell2=zeros(n_len,numiter,q_len);
corr_match_unseed=zeros(n_len,numiter,q_len);
corr_match_slp=zeros(n_len,numiter,q_len);
corr_match_slp_ord=zeros(n_len,numiter,q_len);
corr_match_slp_YiCao=zeros(n_len,numiter,q_len);

obj_func_final_vals_JV=zeros(n_len,numiter,q_len);
obj_func_final_vals_ell2=zeros(n_len,numiter,q_len);
obj_func_final_vals_ell1=zeros(n_len,numiter,q_len);

obj_func_final_vals_proj_JV=zeros(n_len,numiter,q_len);
obj_func_final_vals_proj_ell2=zeros(n_len,numiter,q_len);
obj_func_final_vals_proj_ell1=zeros(n_len,numiter,q_len);
rng_state=cell(numiter,1);

for i=1:numiter
    rng shuffle
    rng_st=rng;
    rng_state{i}=rng_st.State;
    i
    
    for q_i= 1:q_len
        found=0
        q_val=q(q_i)
        
        Bernoulli=rand(N);
        A=rand(N)<Bernoulli;
        A=A-triu(A);A=A+A';
        B=A;
        B=bitflip(A,q_val);
        B=B-triu(B);B=B+B';
        
        ordering=randperm(N);
        for n_i = 1:n_len
            n_i
            n_val_for_i= n_vals(n_i)
            test_ind = (n_val_for_i+1):N;
            test_v =ordering(test_ind);
            tic;
            [matching, iter, ~,fval_JV,fval_proj_JV,P_jv,P_proj_jv]=ConVogHard_rQAP_order(A,B,n_val_for_i,ordering,1);
            
            %P_jv_cell{n_i,i}=P_jv(:,ordering);
            
            %P_jv_pr_cell{n_i,i}=P_proj_jv(:,ordering);
            %running_time_FAQ(i,n_i,q_i)=toc;
            matched_v=matching(test_ind);
            corr_match(n_i,i,q_i) =  sum(matched_v==test_v);
            obj_func_final_vals_JV(n_i,i,q_i)=fval_JV;
            obj_func_final_vals_proj_JV(n_i,i,q_i)=fval_proj_JV;
            
            %             [matching_ell2, iter,~,fval_ell2, fval_proj_ell2,P_l2,P_proj_l2]=seedgraphmatchell2_order(A,B,n_val_for_i,ordering,1);
            %
            %             P_l2_cell{n_i,i}=P_l2(:,ordering);
            %             P_l2_pr_cell{n_i,i}=P_proj_l2(:,ordering);
            %
            %             corr_match_ell2(n_i,i,q_i) =  sum(matching_ell2((n_val_for_i+1):N)==ordering((n_val_for_i+1):N));
            %             obj_func_final_vals_ell2(n_i,i,q_i)=fval_ell2;
            %             obj_func_final_vals_proj_ell2(n_i,i,q_i)=fval_proj_ell2;
            if (q_i==q_int && i<=slp_iter )
                tic;
                [matching_slp,fval_ell1,fval_proj_ell1,P_l1,P_proj_l1] = seedgraphmatchell1(A(ordering,ordering),B(ordering,ordering),n_val_for_i);
                %P_l1_cell{n_i,i} =  P_l1;
                %P_l1_pr_cell{n_i,i}= P_proj_l1;
               % running_time_SLP(i,q_i,n_i)=toc;
               matched_slp_v=matching_slp(test_ind);
                corr_match_slp(n_i,i,q_i) =  sum(matched_slp_v==test_ind);
                obj_func_final_vals_ell1(n_i,i,q_i)=fval_ell1;
                obj_func_final_vals_proj_ell1(n_i,i,q_i)=fval_proj_ell1;
                
                if (fval_ell1-(1E-5)>fval_JV)
                    
                    'Interesting! This should not happen. ell_1 fval is larger than ell_2'
                    
                end
            end
            
            A_sub = A(test_v,test_v);
            B_sub=  B(test_v,test_v);
            matching_unseed=ConVogHard_rQAP(A_sub,B_sub,0);
            corr_match_unseed(n_i,i,q_i) =  sum(matching_unseed==1:(N-n_val_for_i));
        end
    end
end

 pc = corr_match;
 
fc= pc./repmat((N-n_vals'),[1 numiter length(q)]);

pc_slp = corr_match_slp;

fc_slp= pc_slp./repmat((N-n_vals'),[1 numiter length(q)]);


pc_slp_ord = corr_match_slp_ord;

fc_slp_ord= pc_slp_ord./repmat((N-n_vals'),[1 numiter length(q)]);

pc_slp_YiCao = corr_match_slp_YiCao;

fc_slp_YiCao= pc_slp_YiCao./repmat((N-n_vals'),[1 numiter length(q)]);

pc_ell2 = corr_match_ell2;

fc_ell2= pc_ell2./repmat((N-n_vals'),[1 numiter length(q)]);

pc_unseed=corr_match_unseed;
fc_unseed= pc_unseed./repmat((N-n_vals'),[1 numiter length(q)]);
fname = strcat('sim_bitflip_', num2str(job_i) , '.mat')
save(fname)

end

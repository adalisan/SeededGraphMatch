
q= [0:0.05:0.5 ];
q_len = length(q);



q_int=find(q==0.3);
truematch_rqap = zeros(q_len,1);



N=20;
numiter=100;
slp_iter =100;
%n_vals=[0:25 30 35 40 45];
%n_vals=[0:1:5 6:2:20 20:5:35 ];
%n_vals= [0:20 22 24 26];
n_vals=0:16
n_vals=n_vals(n_vals<N);
P_jv_found=[]
P_l1_found=[]
P_l2_found=[]

P_jv_cell=cell(length(n_vals),numiter);

P_jv_pr_cell=cell(length(n_vals),numiter);

P_l2_cell=cell(length(n_vals),numiter);
P_l2_pr_cell=cell(length(n_vals),numiter);
P_l1_cell=cell(length(n_vals),numiter);
P_l1_pr_cell=cell(length(n_vals),numiter);

running_time_FAQ=zeros(numiter,length(n_vals),q_len);
running_time_SLP=zeros(numiter,length(n_vals),q_len);

corr_match=zeros(length(n_vals),numiter,q_len);

corr_match_ell2=zeros(length(n_vals),numiter,q_len);
corr_match_unseed=zeros(length(n_vals),numiter,q_len);
corr_match_slp=zeros(length(n_vals),numiter,q_len);
corr_match_slp_ord=zeros(length(n_vals),numiter,q_len);
corr_match_slp_YiCao=zeros(length(n_vals),numiter,q_len);

obj_func_final_vals_JV=zeros(length(n_vals),numiter,q_len);
obj_func_final_vals_ell2=zeros(length(n_vals),numiter,q_len);
obj_func_final_vals_ell1=zeros(length(n_vals),numiter,q_len);

obj_func_final_vals_proj_JV=zeros(length(n_vals),numiter,q_len);
obj_func_final_vals_proj_ell2=zeros(length(n_vals),numiter,q_len);
obj_func_final_vals_proj_ell1=zeros(length(n_vals),numiter,q_len);


for q_i= 1:length(q)
    found=0    
    q_val=q(q_i)
    for i=1:numiter
        i
        Bernoulli=rand(N);
        A=rand(N)<Bernoulli;
        A=A-triu(A);A=A+A';
        B=A;
        B=bitflip(A,q_val);
        B=B-triu(B);B=B+B';
        
        ordering=randperm(N);
        for n_i = 1:length(n_vals)
            n_i
            tic;
            [matching, iter, ~,fval_JV,fval_proj_JV,P_jv,P_proj_jv]=ConVogHard_rQAP_order(A,B,n_vals(n_i),ordering,1);
            
            P_jv_cell{n_i,i}=P_jv(:,ordering);
            
            P_jv_pr_cell{n_i,i}=P_proj_jv(:,ordering);
            running_time_FAQ(i,n_i,q_i)=toc;
            corr_match(n_i,i,q_i) =  sum(matching((n_vals(n_i)+1):N)==ordering((n_vals(n_i)+1):N));
            obj_func_final_vals_JV(n_i,i,q_i)=fval_JV;
            obj_func_final_vals_proj_JV(n_i,i,q_i)=fval_proj_JV;
           
            [matching_ell2, iter,~,fval_ell2, fval_proj_ell2,P_l2,P_proj_l2]=seedgraphmatchell2_order(A,B,n_vals(n_i),ordering,1);
           
            P_l2_cell{n_i,i}=P_l2(:,ordering);
            P_l2_pr_cell{n_i,i}=P_proj_l2(:,ordering);

            corr_match_ell2(n_i,i,q_i) =  sum(matching_ell2((n_vals(n_i)+1):N)==ordering((n_vals(n_i)+1):N));
            obj_func_final_vals_ell2(n_i,i,q_i)=fval_ell2;
            obj_func_final_vals_proj_ell2(n_i,i,q_i)=fval_proj_ell2;
            if (q_i==q_int && i<=slp_iter )
            tic;
            [matching_slp,fval_ell1,fval_proj_ell1,P_l1,P_proj_l1] = seedgraphmatchell1(A(ordering,ordering),B(ordering,ordering),n_vals(n_i));
            P_l1_cell{n_i,i} =  P_l1;
            P_l1_pr_cell{n_i,i}= P_proj_l1;
            running_time_SLP(i,n_i,q_i)=toc;
            corr_match_slp(n_i,i,q_i) =  sum(matching_slp((n_vals(n_i)+1):N)==(n_vals(n_i)+1):N);
            obj_func_final_vals_ell1(n_i,i,q_i)=fval_ell1;
            obj_func_final_vals_proj_ell1(n_i,i,q_i)=fval_proj_ell1;
            
            if (fval_ell1-(1E-5)>fval_JV)
               A_f=A(ordering,ordering)
               B_f=B(ordering,ordering)
               n_i
               n_vals(n_i)
               P_jv_found=P_jv
               P_l2_found=P_l2
               P_l1_found=P_l1
               P_int= P_l1_found-P_l2_found;
                P_int(abs(P_int)<5E-10)=0;
               found=1
               fval_ell1
               fval_JV
               
                break
            end
            %matching_slp_ord = seedgraphmatchell1_order(A,B,n_vals(n_i),ordering);
            %corr_match_slp_ord(n_i,i,q_i) =  sum(matching_slp_ord((n_vals(n_i)+1):N) ...
            %    ==ordering((n_vals(n_i)+1):N));
            
            %   matching_slp_YiCao = seedgraphmatchell1_YiCao(A(ordering,ordering),B(ordering,ordering),n_vals(n_i));
           %corr_match_slp_YiCao(n_i,i,q_i) =  sum(matching_slp((n_vals(n_i)+1):N)==(n_vals(n_i)+1):N);
            if found==1 
                break
            end
            
            
            end
            test_v =ordering((n_vals(n_i)+1):N);
            A_sub = A(test_v,test_v);
            B_sub=  B(test_v,test_v);
            matching_unseed=ConVogHard_rQAP(A_sub,B_sub,0);
            corr_match_unseed(n_i,i,q_i) =  sum(matching_unseed==1:(N-n_vals(n_i)));
            
        save('./sim_result-20.mat')   
        if found==1 
            break
        end
        end
        if found==1 
            break
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

main_colors = { 'r-' 'g-' 'b-'  'm-'   'y-' 'c-' 'k-.'};

%selected_runs = randi([1 75],[50 1]);


figure

figcolors= colormap(jet);
[num_colors,~]=size(figcolors);
incr=floor(num_colors/q_len);
for i= 1:q_len
    q_i=q(i);
avg_line=mean(fc(:,:,i),2);
sd_line = std(fc(:,:,i),1,2);
    %plot (n_vals(1:5),avg_line(1:5,:),colors{i},'LineWidth',2)
    errorbar (n_vals,avg_line,2*sd_line/sqrt(numiter),'Color',figcolors(i*incr,:),'LineWidth',2)
    hold on
end

xlabel('$m$','Interpreter','latex','FontSize',20)
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)
plot(n_vals,1./(N-n_vals),main_colors{length(main_colors)},'LineWidth',2)

avg_line=mean(fc_unseed(:,:,q_int),2);
sd_line = std(fc_unseed(:,:,q_int),1,2);

errorbar (n_vals,avg_line,2*sd_line/sqrt(numiter),'Color',figcolors(q_int*incr,:), ...
    'LineStyle','-.','LineWidth',1.5)
qvals= num2str(q');


legend(qvals)   
title('Simulation','FontSize',20)    
xlim([-0.5 max(n_vals)+0.5])
ylim([-0.1 1.1])


slp_plot=mean(fc_slp(:,1:slp_iter,q_int),2);
slp_plot_sd=std(fc_slp(:,1:slp_iter,q_int),1,2);

ell2_plot=mean(fc_ell2(:,1:slp_iter,q_int),2);
ell2_plot_sd=std(fc_ell2(:,1:slp_iter,q_int),1,2);


%slp_plot_YiCao=mean(fc_slp_YiCao(:,1:slp_iter,1),2);
%slp_plot_YiCao_sd=std(fc_slp_YiCao(:,1:slp_iter,1),1,2);


rqap_plot=mean(fc(:,1:slp_iter,q_int),2);
rqap_plot_sd=std(fc(:,1:slp_iter,q_int),1,2);
figure
errorbar(n_vals,slp_plot,slp_plot_sd/sqrt(slp_iter),'r-')

hold on
errorbar(n_vals,rqap_plot,rqap_plot_sd/sqrt(slp_iter),'b-')
hold on


errorbar(n_vals,ell2_plot,ell2_plot_sd/sqrt(slp_iter),'m-')
%errorbar(n_vals,slp_plot_YiCao,slp_plot_YiCao_sd/sqrt(slp_iter),'m--')
%hold on

%errorbar(n_vals,slp_plot_ord,slp_plot_ord_sd/sqrt(slp_iter),'g:')

xlabel('$m$','Interpreter','latex','FontSize',20)
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)
legend('SLP','FAQ/rQAP','FishkindFAQ')

[ro_JV,co_JV]=find((obj_func_final_vals_JV(:,:,q_int)+5E-3)<obj_func_final_vals_ell1(:,:,q_int));

[ro_JV_p,co_JV_p]=find((obj_func_final_vals_proj_JV(:,:,q_int)+5E-3)<obj_func_final_vals_proj_ell1(:,:,q_int));

[ro_ell2,co_ell2]=find((obj_func_final_vals_ell2(:,:,q_int)+5E-3)<obj_func_final_vals_ell1(:,:,q_int));

[ro_ell2_p,co_ell2_p]=find((obj_func_final_vals_proj_ell2(:,:,q_int)+5E-3)<obj_func_final_vals_proj_ell1(:,:,q_int));


indices_JV=sub2ind(size(squeeze(obj_func_final_vals_JV)),ro_JV,co_JV,q_int);
[obj_func_final_vals_JV(indices_JV) obj_func_final_vals_ell1(indices_JV)]
indices_JV_p=sub2ind(size(squeeze(obj_func_final_vals_JV)),ro_JV_p,co_JV_p,q_int);
[obj_func_final_vals_proj_JV(indices_JV_p) obj_func_final_vals_proj_ell1(indices_JV_p)]

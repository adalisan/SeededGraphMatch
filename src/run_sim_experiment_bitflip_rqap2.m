
%rho= [0:0.05:0.5 ];
rho=[0 0.1 0.3 0.45 0.5];
%rho=[0 0.1 0.3 0.45 0.5];
rho_len = length(rho);



rho_int=find(rho==0.3);



N=300;
numiter=5;
slp_iter =-1;
n_vals=[0:39 40:2:98 100:5:195 200:5:275  ];
%n_vals=[0:1:5 6:2:20 20:5:35 ];
%n_vals= [0:20 22 24 26];
%n_vals=0:16
%n_vals=[ 0 1 5 15 17:2:25 30 35 40 60 80 100 300 450];
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

running_time_FAQ=zeros(numiter,length(n_vals),rho_len);
running_time_SLP=zeros(numiter,length(n_vals),rho_len);

corr_match=zeros(length(n_vals),numiter,rho_len);

corr_match_ell2=zeros(length(n_vals),numiter,rho_len);
corr_match_unseed=zeros(length(n_vals),numiter,rho_len);
corr_match_slp=zeros(length(n_vals),numiter,rho_len);
corr_match_slp_ord=zeros(length(n_vals),numiter,rho_len);
corr_match_slp_YiCao=zeros(length(n_vals),numiter,rho_len);

obj_func_final_vals_JV=zeros(length(n_vals),numiter,rho_len);
obj_func_final_vals_ell2=zeros(length(n_vals),numiter,rho_len);
obj_func_final_vals_ell1=zeros(length(n_vals),numiter,rho_len);

obj_func_final_vals_proj_JV=zeros(length(n_vals),numiter,rho_len);
obj_func_final_vals_proj_ell2=zeros(length(n_vals),numiter,rho_len);
obj_func_final_vals_proj_ell1=zeros(length(n_vals),numiter,rho_len);


for rho_i= 1:rho_len
    found=0    
    rho_val=rho(rho_i)
    for i=1:numiter
        i
        Bernoulli=rand(N);
        A=rand(N)<Bernoulli;
        A=A-triu(A);A=A+A';
        B=A;
        B=bitflip(A,rho_val);
        B=B-triu(B);B=B+B';
        
        ordering=randperm(N);
        for n_i = 1:length(n_vals)
            n_i
            tic;
            [matching, iter]=ConVogHard_rQAP2_order(A,B,n_vals(n_i),ordering,1);
            
            %P_jv_cell{n_i,i}=P_jv(:,ordering);
            
            %P_jv_pr_cell{n_i,i}=P_proj_jv(:,ordering);
            running_time_FAQ(i,n_i,rho_i)=toc;
            corr_match(n_i,i,rho_i) =  sum(matching((n_vals(n_i)+1):N)==ordering((n_vals(n_i)+1):N));
            %obj_func_final_vals_JV(n_i,i,rho_i)=fval_JV;
            %obj_func_final_vals_proj_JV(n_i,i,rho_i)=fval_proj_JV;
           
%             [matching_ell2, iter,~,fval_ell2, fval_proj_ell2,P_l2,P_proj_l2]=seedgraphmatchell2_order(A,B,n_vals(n_i),ordering,1);
%            
%             P_l2_cell{n_i,i}=P_l2(:,ordering);
%             P_l2_pr_cell{n_i,i}=P_proj_l2(:,ordering);
% 
%             corr_match_ell2(n_i,i,rho_i) =  sum(matching_ell2((n_vals(n_i)+1):N)==ordering((n_vals(n_i)+1):N));
%             obj_func_final_vals_ell2(n_i,i,rho_i)=fval_ell2;
%             obj_func_final_vals_proj_ell2(n_i,i,rho_i)=fval_proj_ell2;
             if (rho_i==rho_int && i<=slp_iter )
            tic;
            [matching_slp,fval_ell1,fval_proj_ell1,P_l1,P_proj_l1] = seedgraphmatchell1(A(ordering,ordering),B(ordering,ordering),n_vals(n_i));
            %P_l1_cell{n_i,i} =  P_l1;
            %P_l1_pr_cell{n_i,i}= P_proj_l1;
            running_time_SLP(i,n_i,rho_i)=toc;
            corr_match_slp(n_i,i,rho_i) =  sum(matching_slp((n_vals(n_i)+1):N)==(n_vals(n_i)+1):N);
            obj_func_final_vals_ell1(n_i,i,rho_i)=fval_ell1;
            obj_func_final_vals_proj_ell1(n_i,i,rho_i)=fval_proj_ell1;
            
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
            %corr_match_slp_ord(n_i,i,rho_i) =  sum(matching_slp_ord((n_vals(n_i)+1):N) ...
            %    ==ordering((n_vals(n_i)+1):N));
            
            %   matching_slp_YiCao = seedgraphmatchell1_YiCao(A(ordering,ordering),B(ordering,ordering),n_vals(n_i));
           %corr_match_slp_YiCao(n_i,i,rho_i) =  sum(matching_slp((n_vals(n_i)+1):N)==(n_vals(n_i)+1):N);
            if found==1 
                break
            end
            
            
            end
            test_v =ordering((n_vals(n_i)+1):N);
            A_sub = A(test_v,test_v);
            B_sub=  B(test_v,test_v);
            matching_unseed=ConVogHard_rQAP(A_sub,B_sub,0);
            corr_match_unseed(n_i,i,rho_i) =  sum(matching_unseed==1:(N-n_vals(n_i)));
            
            save(strcat('./rqap2_sim_result_',num2str(N),'.mat'))
      
        end
        if found==1 
            break
        end
        
    end
end

pc = corr_match;
 
fc= pc./repmat((N-n_vals'),[1 numiter rho_len]);

pc_slp = corr_match_slp;
 
fc_slp= pc_slp./repmat((N-n_vals'),[1 numiter rho_len]);


pc_slp_ord = corr_match_slp_ord;
 
fc_slp_ord= pc_slp_ord./repmat((N-n_vals'),[1 numiter rho_len]);

pc_slp_YiCao = corr_match_slp_YiCao;
 
fc_slp_YiCao= pc_slp_YiCao./repmat((N-n_vals'),[1 numiter rho_len]);

pc_ell2 = corr_match_ell2;
 
fc_ell2= pc_ell2./repmat((N-n_vals'),[1 numiter rho_len]);

 pc_unseed=corr_match_unseed;
fc_unseed= pc_unseed./repmat((N-n_vals'),[1 numiter rho_len]);

main_colors = { 'r-' 'g-' 'b-'  'm-'   'y-' 'c-' 'k-.'};

%selected_runs = randi([1 75],[50 1]);


figure

figcolors= colormap(jet);
[num_colors,~]=size(figcolors);
incr=floor(num_colors/rho_len);
for i= 1:rho_len
    rho_i=rho(i);
avg_line=mean(fc(:,:,i),2);
sd_line = std(fc(:,:,i),1,2);
   % plot (n_vals,avg_line,'Color',figcolors(i*incr,:),'LineWidth',2)
    errorbar (n_vals,avg_line,2*sd_line/sqrt(numiter),'Color',figcolors(i*incr,:),'LineWidth',2)
    hold on
end

xlabel('$m$','Interpreter','latex','FontSize',20)
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)
plot(n_vals,1./(N-n_vals),main_colors{length(main_colors)},'LineWidth',2)

%avg_line=mean(fc_unseed(:,:,rho_int),2);
%sd_line = std(fc_unseed(:,:,rho_int),1,2);

%errorbar (n_vals,avg_line,2*sd_line/sqrt(numiter),'Color',figcolors(rho_int*incr,:), ...
%    'LineStyle','-.','LineWidth',1.5)
qvals= num2str(rho');


legend(qvals)   
title('Simulation','FontSize',20)    
xlim([-5 max(n_vals)+2])
ylim([-0.1 1.1])

figure
rho_i=rho(i);
avg_line=mean(fc(:,:,rho_int),2);
sd_line = std(fc(:,:,rho_int),1,2);
    %plot (n_vals(1:5),avg_line(1:5,:),colors{i},'LineWidth',2)
    errorbar (n_vals,avg_line,2*sd_line/sqrt(numiter),'Color','r','LineWidth',2)
    hold on

    
avg_line=mean(fc_unseed(:,:,rho_int),2);
sd_line = std(fc_unseed(:,:,rho_int),1,2);

    
errorbar (n_vals,avg_line,2*sd_line/sqrt(numiter),'Color','r', ...
    'LineStyle','-.','LineWidth',1.5)

title('Simulation   ($\rho=0.3$)','Interpreter','latex','FontSize',20)    
xlim([-0.5 max(n_vals)+0.5])
ylim([-0.1 1.1])
xlabel('$m$','Interpreter','latex','FontSize',20)
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)

legend('Seeded', 'Unseeded')






slp_plot=mean(fc_slp(:,1:slp_iter,rho_int),2);
slp_plot_sd=std(fc_slp(:,1:slp_iter,rho_int),1,2);

ell2_plot=mean(fc_ell2(:,1:slp_iter,rho_int),2);
ell2_plot_sd=std(fc_ell2(:,1:slp_iter,rho_int),1,2);


%slp_plot_YiCao=mean(fc_slp_YiCao(:,1:slp_iter,1),2);
%slp_plot_YiCao_sd=std(fc_slp_YiCao(:,1:slp_iter,1),1,2);


rqap_plot=mean(fc(:,1:slp_iter,rho_int),2);
rqap_plot_sd=std(fc(:,1:slp_iter,rho_int),1,2);
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
legend('SLP','rQAP2','FishkindFAQ')

[ro_JV,co_JV]=find((obj_func_final_vals_JV(:,:,rho_int)+5E-3)<obj_func_final_vals_ell1(:,:,rho_int));

[ro_JV_p,co_JV_p]=find((obj_func_final_vals_proj_JV(:,:,rho_int)+5E-3)<obj_func_final_vals_proj_ell1(:,:,rho_int));

[ro_ell2,co_ell2]=find((obj_func_final_vals_ell2(:,:,rho_int)+5E-3)<obj_func_final_vals_ell1(:,:,rho_int));

[ro_ell2_p,co_ell2_p]=find((obj_func_final_vals_proj_ell2(:,:,rho_int)+5E-3)<obj_func_final_vals_proj_ell1(:,:,rho_int));


 mean([running_time_FAQ(:,1,rho_int)  running_time_SLP(:,1,rho_int)],1)


std([running_time_FAQ(:,1,rho_int)  running_time_SLP(:,1,rho_int)],1,1)/sqrt(numiter)






%indices_JV=sub2ind(size(squeeze(obj_func_final_vals_JV)),ro_JV,co_JV,rho_int);
%[obj_func_final_vals_JV(indices_JV) obj_func_final_vals_ell1(indices_JV)]
%indices_JV_p=sub2ind(size(squeeze(obj_func_final_vals_JV)),ro_JV_p,co_JV_p,rho_int);
%[obj_func_final_vals_proj_JV(indices_JV_p) obj_func_final_vals_proj_ell1(indices_JV_p)]

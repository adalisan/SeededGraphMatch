% rho -- Perturbation parameter values
rho= 0:0.05:0.5 ;
rho_len = length(rho);

% run SLP only for rho_int value
rho_int=find(rho==0.3);



N=300;
numiter=5;

%run SLP up for iterations from 1 up to slp_iter
slp_iter =-1;

%n_vals  -- number of hard seeds
n_vals=[0:39 40:2:98 100:5:195 200:5:275  ];
n_vals=n_vals(n_vals<N);


running_time_FAQ=zeros(numiter,length(n_vals),rho_len);
running_time_SLP=zeros(numiter,length(n_vals),rho_len);

corr_match=zeros(length(n_vals),numiter,rho_len);

corr_match_ell2=zeros(length(n_vals),numiter,rho_len);
corr_match_unseed=zeros(length(n_vals),numiter,rho_len);
corr_match_slp=zeros(length(n_vals),numiter,rho_len);

obj_func_final_vals_JV=zeros(length(n_vals),numiter,rho_len);
obj_func_final_vals_ell2=zeros(length(n_vals),numiter,rho_len);
obj_func_final_vals_ell1=zeros(length(n_vals),numiter,rho_len);

obj_func_final_vals_proj_JV=zeros(length(n_vals),numiter,rho_len);
obj_func_final_vals_proj_ell2=zeros(length(n_vals),numiter,rho_len);
obj_func_final_vals_proj_ell1=zeros(length(n_vals),numiter,rho_len);


for rho_i= 1:length(rho)
    
    rho_val=rho(rho_i)
    for i=1:numiter
        i
        Bernoulli=rand(N);
        A=rand(N)<Bernoulli;
        A=A-triu(A);A=A+A';
        B=A;
        B=bitflip(A,rho_val);
        B=B-triu(B);B=B+B';
        
        %Create a random ordering such that the first m vertices are seeds
        ordering=randperm(N);
        for n_i = 1:length(n_vals)
            n_i
            tic;
            [matching, iter, ~,fval_JV,fval_proj_JV]=ConVogHard_rQAP_order(A,B,n_vals(n_i),ordering,1);
            
            running_time_FAQ(i,n_i,rho_i)=toc;
            corr_match(n_i,i,rho_i) =  sum(matching((n_vals(n_i)+1):N)==ordering((n_vals(n_i)+1):N));
            obj_func_final_vals_JV(n_i,i,rho_i)=fval_JV;
            obj_func_final_vals_proj_JV(n_i,i,rho_i)=fval_proj_JV;
            
            if (rho_i==rho_int & i<=slp_iter )
                tic;
                [matching_slp,fval_ell1,fval_proj_ell1,P_l1,P_proj_l1] = seedgraphmatchell1(A(ordering,ordering),B(ordering,ordering),n_vals(n_i));
                running_time_SLP(i,n_i,rho_i)=toc;
                corr_match_slp(n_i,i,rho_i) =  sum(matching_slp((n_vals(n_i)+1):N)==(n_vals(n_i)+1):N);
                obj_func_final_vals_ell1(n_i,i,rho_i)=fval_ell1;
                obj_func_final_vals_proj_ell1(n_i,i,rho_i)=fval_proj_ell1;
                
                
            end
            
            % Match with no seeds to highlight how much seeding helps.
            test_v =ordering((n_vals(n_i)+1):N);
            A_sub = A(test_v,test_v);
            B_sub=  B(test_v,test_v);
            matching_unseed=ConVogHard_rQAP(A_sub,B_sub,0);
            corr_match_unseed(n_i,i,rho_i) =  sum(matching_unseed==1:(N-n_vals(n_i)));
            
            save(strcat('./sim_result_',num2str(N),'.mat'))
        end
        
    end
end

pc = corr_match;

fc= pc./repmat((N-n_vals'),[1 numiter length(rho)]);

pc_slp = corr_match_slp;

fc_slp= pc_slp./repmat((N-n_vals'),[1 numiter length(rho)]);


pc_unseed=corr_match_unseed;
fc_unseed= pc_unseed./repmat((N-n_vals'),[1 numiter length(rho)]);

main_colors = { 'r-' 'g-' 'b-'  'm-'   'y-' 'c-' 'k-.'};





% Plotting commands


figure

figcolors= colormap(jet);
[num_colors,~]=size(figcolors);
incr=floor(num_colors/rho_len);
for i= 1:rho_len
    
    avg_line=mean(fc(:,:,i),2);
    sd_line = std(fc(:,:,i),1,2);
    % plot (n_vals,avg_line,'Color',figcolors(i*incr,:),'LineWidth',2)
    errorbar (n_vals,avg_line,2*sd_line/sqrt(numiter),'Color',figcolors(i*incr,:),'LineWidth',2)
    hold on
end

xlabel('$m$','Interpreter','latex','FontSize',20)
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)
plot(n_vals,1./(N-n_vals),main_colors{length(main_colors)},'LineWidth',2)

rhovals= num2str(rho');
legend(rhovals)
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



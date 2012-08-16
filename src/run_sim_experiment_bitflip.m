
q= [ 0.1   0.3  0.4 0.5 ];
q_len = length(q);

truematch_rqap = zeros(q_len,1);




N=60;
numiter=50;
n_vals=[0:25 30 35 40 45];

running_time_FAQ=zeros(numiter,length(n_vals),q_len);
running_time_SLP=zeros(numiter,length(n_vals),q_len);

corr_match=zeros(length(n_vals),numiter,q_len);
corr_match_unseed=zeros(length(n_vals),numiter,q_len);
corr_match_slp=zeros(length(n_vals),numiter,q_len);
for q_i= 1:length(q)
    
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
            
            tic;
            matching=ConVogHard_rQAP_order(A,B,n_vals(n_i),ordering);
            running_time_FAQ(i,n_i,q_i)=toc;
            corr_match(n_i,i,q_i) =  sum(matching((n_vals(n_i)+1):N)==ordering((n_vals(n_i)+1):N));
            if (q_i==2 && i<4)
            tic;
            matching_slp = seedgraphmatchell1(A(ordering,ordering),B(ordering,ordering),n_vals(n_i));
            running_time_SLP(i,n_i,q_i)=toc;
            corr_match_slp(n_i,i,q_i) =  sum(matching_slp((n_vals(n_i)+1):N)==(n_vals(n_i)+1):N);
            end
            test_v =ordering((n_vals(n_i)+1):N);
            A_sub = A(test_v,test_v);
            B_sub=  B(test_v,test_v);
            matching_unseed=ConVogHard_rQAP(A_sub,B_sub,0);
            corr_match_unseed(n_i,i,q_i) =  sum(matching_unseed==1:(N-n_vals(n_i)));
            
        save('sim_result-60.mat')   
        end
        
    end
end
pc = corr_match;
 
fc= pc./repmat((N-n_vals'),[1 numiter length(q)]);

colors = { 'r-' 'g-' 'b-'  'm-'   'y'  'k-.'};

%selected_runs = randi([1 75],[50 1]);


figure

for i= 1:q_len
    q_i=q(i);
avg_line=mean(fc(:,:,i),2);
sd_line = std(fc(:,:,i),1,2);
    %plot (n_vals(1:5),avg_line(1:5,:),colors{i},'LineWidth',2)
    errorbar (n_vals,avg_line,2*sd_line/sqrt(numiter),colors{i},'LineWidth',2)
    hold on
end

xlabel('$m$','Interpreter','latex','FontSize',20)
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)
plot(n_vals,1./(N-n_vals),colors{length(colors)},'LineWidth',2)
qvals= num2str(q');


legend(qvals)   
title('Simulation','FontSize',20)    
xlim([-2 max(n_vals)])







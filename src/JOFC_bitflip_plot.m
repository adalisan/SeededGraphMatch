
figure
rho_len= 6
numiter = 20
n_vals = [10 15 20 20 30 40 50 60 70 80 90 ];
figcolors= colormap(jet);
[num_colors,~]=size(figcolors);
incr=floor(num_colors/rho_len);
for i= 1:rho_len
 
    avg_line=mean(corr_results_agg(:,i,:),3);
    sd_line = std(corr_results_agg(:,i,:),1,3);
   % plot (n_vals,avg_line,'Color',figcolors(i*incr,:),'LineWidth',2)
    errorbar (n_vals,avg_line,2*sd_line/sqrt(numiter),'Color',figcolors(i*incr,:),'LineWidth',2)
    hold on
end

xlabel('$m$','Interpreter','latex','FontSize',20)
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)
plot(n_vals,1./(N-n_vals),main_colors{length(main_colors)},'LineWidth',2)

%avg_line=mean(fc_unseed(:,:,rho_int),2);
%sd_line = std(fc_unseed(:,:,rho_int),1,2);
rho = 0:0.1:0.5;
%errorbar (n_vals,avg_line,2*sd_line/sqrt(numiter),'Color',figcolors(rho_int*incr,:), ...
%    'LineStyle','-.','LineWidth',1.5)
qvals= num2str(rho');


legend(qvals)   
title('Simulation','FontSize',20)    
xlim([-10 max(n_vals)+2])
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
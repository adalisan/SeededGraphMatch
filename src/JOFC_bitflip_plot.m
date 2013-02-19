N= 100;

rho_len= length(rho);
numiter = 20;
n_vals =n_vals_bitflip;
figcolors= colormap(jet);
[num_colors,~]=size(figcolors);
incr=floor(num_colors/rho_len);
main_colors = { 'r' 'g' 'b'  'm'   'y' 'c' 'k'};


for i= 1:rho_len
 
    avg_line = JOFC_corr_bitflip(:,i,4)./(N-n_vals);
    sd_line = JOFC_corr_bitflip_sd(:,i,4)./(N-n_vals);
   % plot (n_vals,avg_line,'Color',figcolors(i*incr,:),'LineWidth',2)
    errorbar (n_vals,avg_line,2*sd_line,'Color',figcolors(i*incr,:),'LineWidth',2)
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
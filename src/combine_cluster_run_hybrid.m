%s = what; %look in current directory
s=what('./cache/cluster_run/hybrid/') %change dir for your directory name 
matfiles=s.mat;
fc_rqap_agg=[];
fc_hybrid_agg=[];
fc_slp_agg=[];

for a=1:numel(matfiles)

if a==1
    load(char(matfiles(a)))
    fc_rqap_agg=fc_rqap;
    fc_hybrid_agg=fc;
    fc_slp_agg=fc_slp;
end
    try
        load(char(matfiles(a)),'fc','fc_rqap','fc_slp')
        
    catch
        'unable to load'
        char(matfiles(a))
        fc=[];
        fc_slp=[];
        fc_rqap=[];
    end
fc_hybrid_agg=cat(2,fc_hybrid_agg,fc);
fc_rqap_agg=cat(2,fc_rqap_agg,fc_rqap);
fc_slp_agg=cat(2,fc_slp_agg,fc_slp);
end

[~,numiter,~]= size(fc_agg)


figure

figcolors= colormap(hsv);
[num_colors,~]=size(figcolors);
incr=floor(num_colors/q_len);
for i= 1:q_len
    q_i=q(i);
avg_line=mean(fc_hybrid_agg(:,:,i),2);
sd_line = std(fc_hybrid_agg(:,:,i),1,2);
    plot (n_vals,avg_line,'Color',figcolors(i*incr,:),'LineWidth',2)
    errorbar (n_vals,avg_line,2*sd_line/sqrt(numiter), ...
        'Color',figcolors(i*incr,:),'LineWidth',2,'LineStyle','-')
    hold on
    
    avg_line_2=mean(fc_rqap_agg(:,:,i),2);
sd_line_2 = std(fc_rqap_agg(:,:,i),1,2);
errorbar (n_vals,avg_line_2,2*sd_line_2/sqrt(numiter), ...
    'Color',figcolors(i*incr,:),'LineWidth',2,'LineStyle',':')
end




main_colors = { 'r-' 'g-' 'b-'  'm-'   'y-' 'c-' 'k-.'};


xlabel('$m$','Interpreter','latex','FontSize',20)
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)
plot(n_vals,1./(N-n_vals),main_colors{length(main_colors)},'LineWidth',2)

%avg_line=mean(fc_unseed(:,:,q_int),2);
%sd_line = std(fc_unseed(:,:,q_int),1,2);

%errorbar (n_vals,avg_line,2*sd_line/sqrt(numiter),'Color',figcolors(q_int*incr,:), ...
%    'LineStyle','-.','LineWidth',1.5)
qvals= num2str(q');


 
title('Simulation','FontSize',20)    
xlim([-5 max(n_vals)+20])
ylim([-0.1 1.1])
legend(qvals)  
legend('hybrid-FAQ','FAQ' )



figure
for i= 1:q_len
    q_i=q(i);
avg_line=mean(fc_agg(:,:,i),2);
sd_line = std(fc_agg(:,:,i),1,2);
    plot (n_vals,avg_line,'Color',figcolors(i*incr,:),'LineWidth',2)
    hold on
    
    avg_line_2=mean(fc_rqap_agg(:,:,i),2);
sd_line_2 = std(fc_rqap_agg(:,:,i),1,2);
plot (n_vals,avg_line_2, ...
    'Color',figcolors(i*incr,:),'LineWidth',2,'LineStyle',':')
end





main_colors = { 'r-' 'g-' 'b-'  'm-'   'y-' 'c-' 'k-.'};


xlabel('$m$','Interpreter','latex','FontSize',20)
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)
plot(n_vals,1./(N-n_vals),main_colors{length(main_colors)},'LineWidth',2)

%avg_line=mean(fc_unseed(:,:,q_int),2);
%sd_line = std(fc_unseed(:,:,q_int),1,2);

%errorbar (n_vals,avg_line,2*sd_line/sqrt(numiter),'Color',figcolors(q_int*incr,:), ...
%    'LineStyle','-.','LineWidth',1.5)
q_r = kron(q,ones(1,2))
qvals= num2str(q_r');


 
title('Simulation','FontSize',20)    
xlim([-5 max(n_vals)+20])
ylim([-0.1 1.1])
legend(qvals)  






n_vals_s = n_vals(n_vals<30 )
fc_s=fc(1:length(n_vals_s),:,:);


figure

figcolors= colormap(hsv);
[num_colors,~]=size(figcolors);
incr=floor(num_colors/q_len);
for i= 1:q_len
    q_i=q(i);
avg_line=mean(fc_s(:,:,i),2);
sd_line = std(fc_s(:,:,i),1,2);
   % plot (n_vals,avg_line,'Color',figcolors(i*incr,:),'LineWidth',2)
    errorbar (n_vals_s,avg_line,2*sd_line/sqrt(numiter),'Color',figcolors(i*incr,:),'LineWidth',2)
    hold on
end



xlabel('$m$','Interpreter','latex','FontSize',20)
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)
plot(n_vals_s,1./(N-n_vals_s),main_colors{length(main_colors)},'LineWidth',2)

%avg_line=mean(fc_unseed(:,:,q_int),2);
%sd_line = std(fc_unseed(:,:,q_int),1,2);

%errorbar (n_vals,avg_line,2*sd_line/sqrt(numiter),'Color',figcolors(q_int*incr,:), ...
%    'LineStyle','-.','LineWidth',1.5)
qvals= num2str(q');


legend(qvals)   
title('Simulation','FontSize',20)    
xlim([-1 max(n_vals_s)+2])
ylim([-0.1 1.1])







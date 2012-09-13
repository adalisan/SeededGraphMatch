%s = what; %look in current directory
s = what('./cache/cluster_run/wiki/') %change dir for your directory name 
matfiles=s.mat;
fc_agg=[];

for a=1:numel(matfiles)

if a==1
    load(char(matfiles(a)))
    fc_agg = fc;
else
    try
        load(char(matfiles(a)),'corr_match')
    catch
        'unable to load'
        char(matfiles(a))
        fc=[];
    end
    fc_agg=cat(2,fc_agg,fc);
end
end


select=size(fc_agg,2);
fc_agg_s= fc_agg(:,1:select,:);

fc = fc_agg_s;
[~,numiter,~] = size(fc_agg_s) 
%fc= pc./repmat((N-n_vals'),[1 numiter length(q)]);

figure

figcolors= colormap(hsv);
main_colors = { 'r-' 'g-' 'b-'  'm-'   'y-' 'c-' 'k-.'};

[num_colors,~]=size(figcolors);


    avg_line=mean(fc(:,:),2);
    sd_line = std(fc(:,:),1,2);
    % plot (n_vals,avg_line,'Color',figcolors(i*incr,:),'LineWidth',2)
    errorbar (n_vals,avg_line,2*sd_line/sqrt(numiter),main_colors{1},'LineWidth',2)
    hold on





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







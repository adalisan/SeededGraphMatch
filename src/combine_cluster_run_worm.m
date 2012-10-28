%s = what; %look in current directory
s=what('./cache/cluster_run/worm/') %change dir for your directory name 
matfiles=s.mat;
fc_rqap_agg=[];
fc_unwt_agg=[];
fc_seed_agg=[];

for a=1:numel(matfiles)

if a==1
    load(char(matfiles(a)))
    fc_rqap_agg=fc;
    fc_unwt_agg=fc_unwt;
    fc_seed_agg=fc_seed;
end
    try
        load(char(matfiles(a)),'fc','fc_unwt','fc_seed')
        
    catch
        'unable to load'
        char(matfiles(a))
        fc=[];
        fc_seed=[];
        fc_unwt=[];
    end
%     size(fc)
%     size(fc_unwt)
%     size(fc_seed)
%     size(fc_unwt_agg)
%      size(fc_rqap_agg)
%       size(fc_seed_agg)
fc_unwt_agg=cat(2,fc_unwt_agg,fc_unwt);
fc_rqap_agg=cat(2,fc_rqap_agg,fc);
fc_seed_agg=cat(2,fc_seed_agg,fc_seed);
end

[~,numiter,~]= size(fc_unwt_agg)

select = numiter
N=253
figure

figcolors= colormap(hsv);
main_colors = { 'r' 'g' 'b'  'm'   'y' 'c' 'k'};

[num_colors,~]=size(figcolors);


    avg_line=mean(fc_rqap_agg(:,:),2);
    sd_line = std(fc_rqap_agg(:,:),1,2);
    % plot (n_vals,avg_line,'Color',figcolors(i*incr,:),'LineWidth',2)
    errorbar (n_vals,avg_line,2*sd_line/sqrt(select),main_colors{1},'LineWidth',2)
    hold on


    avg_line=mean(fc_unwt_agg(:,:),2);
    sd_line = std(fc_unwt_agg(:,:),1,2);
    % plot (n_vals,avg_line,'Color',figcolors(i*incr,:),'LineWidth',2)
    errorbar (n_vals,avg_line,2*sd_line/sqrt(select),main_colors{2},'LineWidth',2)
    hold on



xlabel('$m$','Interpreter','latex','FontSize',20)
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)
plot(n_vals,1./(N-n_vals),main_colors{length(main_colors)},'LineWidth',2, ...
'LineStyle','-.')

avg_line=mean(fc_seed_agg,2);
sd_line = std(fc_seed_agg,1,2);

errorbar (n_vals',avg_line,2*sd_line/sqrt(numiter),'Color',main_colors{3}, ...
    'LineStyle','-','LineWidth',1.5)


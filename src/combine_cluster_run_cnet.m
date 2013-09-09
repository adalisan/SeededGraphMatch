%s = what; %look in current directory
result_dir='./cache/cluster_run/cnet/'
s = what(result_dir) %change dir for your directory name 
matfiles=s.mat;
fc_agg=[];
fc_seed_agg=[];
N=1382;
for a=1:numel(matfiles)

if a==1
    load(strcat(result_dir,char(matfiles(a))))
    fc_agg = fc;
%     fc_seed_agg = fc_seed;
else
    try
        load(strcat(result_dir,char(matfiles(a))),'fc','fc_seed')
    catch
        'unable to load'
        char(matfiles(a))
        fc=[];
    end
    fc_agg=cat(2,fc_agg,fc);
%     fc_seed_agg=cat(2,fc_seed_agg,fc_seed);
end
end


select=size(fc_agg,2);
fc_agg_s= fc_agg(:,1:select,:);


fc = fc_agg_s;
fc_unseed= fc_seed_agg;

[~,numiter,~] = size(fc_agg_s) 
%fc= pc./repmat((N-n_vals'),[1 numiter length(q)]);

figure

figcolors= colormap(hsv);
main_colors = { 'r' 'g' 'b'  'm'   'y' 'c' 'k'};

[num_colors,~]=size(figcolors);


    avg_line=mean(fc(:,:),2);
    sd_line = std(fc(:,:),1,2);
    % plot (n_vals,avg_line,'Color',figcolors(i*incr,:),'LineWidth',2)
    errorbar (n_vals,avg_line,2*sd_line/sqrt(select),main_colors{1},'LineWidth',2)
    hold on





xlabel('$m$','Interpreter','latex','FontSize',20)
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)
plot(n_vals,1./(N-n_vals),main_colors{length(main_colors)},'LineWidth',2, ...
'LineStyle','-.',)

avg_line=mean(fc_unseed,2);
sd_line = std(fc_unseed,1,2);

errorbar (n_vals',avg_line,2*sd_line/sqrt(numiter),'Color',main_colors{2}, ...
    'LineStyle','-','LineWidth',1.5)


 



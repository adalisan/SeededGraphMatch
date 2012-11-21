%s = what; %look in current directory
%s=what('./cache/cluster_run/worm/questionable (unnormalized) weighted') %change dir for your directory name 
s=what('./cache/cluster_run/worm/')
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

%s=what('./cache/cluster_run/worm/questionable (unnormalized) weighted/directed/') %change dir for your directory name 

matfiles=s.mat;
fc_dir_agg=[];
fc_unwt_dir_agg=[];
fc_unseed_dir_agg=[];


for a=1:numel(matfiles)

if a==1
    load(char(matfiles(a)))
    fc_dir_agg=fc;
    fc_unwt_dir_agg=fc_unwt;
    fc_unseed_dir_agg=fc_seed;
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
fc_unwt_dir_agg=cat(2,fc_unwt_dir_agg,fc_unwt);
fc_dir_agg=cat(2,fc_dir_agg,fc);
fc_unseed_dir_agg=cat(2,fc_unseed_dir_agg,fc_seed);
end






for a=1:numel(matfiles)

if a==1
    load(char(matfiles(a)))
    fc_dir_agg=fc_dir;
    fc_unwt_dir_agg=fc_unwt_dir;
    fc_unseed_dir_agg=fc_unseed_dir;
end
    try
        load(char(matfiles(a)),'fc_dir','fc_unwt_dir','fc_unseed_dir')
        
    catch
        'unable to load'
        char(matfiles(a))
        fc_dir=[];
        fc_unseed_dir=[];
        fc_unwt_dir=[];
    end
%     size(fc)
%     size(fc_unwt)
%     size(fc_seed)
%     size(fc_unwt_agg)
%      size(fc_rqap_agg)
%       size(fc_seed_agg)
fc_unwt_dir_agg=cat(2,fc_unwt_dir_agg,fc_unwt_dir);
fc_dir_agg=cat(2,fc_dir_agg,fc_dir);
fc_unseed_dir_agg=cat(2,fc_unseed_dir_agg,fc_unseed_dir);
end



figure

figcolors= colormap(hsv);
main_colors = { 'r' 'g' 'b'  'm'   'y' 'c' 'k'};

[num_colors,~]=size(figcolors);


    avg_line=mean(fc_rqap_agg(:,:),2);
    sd_line = std(fc_rqap_agg(:,:),1,2);
     d1=plot (n_vals,avg_line,strcat(main_colors{1},'--'),'LineWidth',2.5)
        hold on
    errorbar (n_vals,avg_line,2*sd_line/sqrt(select),strcat(main_colors{1},'.'),'LineWidth',1);
 


    avg_line=mean(fc_unwt_agg(:,:),2);
    sd_line = std(fc_unwt_agg(:,:),1,2);
    d2= plot (n_vals,avg_line,strcat(main_colors{1},'-'),'LineWidth',2.5)
    errorbar (n_vals,avg_line,2*sd_line/sqrt(select),strcat(main_colors{1},'.'),'LineWidth',1);
    hold on



xlabel('$m$','Interpreter','latex','FontSize',20)
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)
d3=plot(n_vals,1./(N-n_vals),main_colors{length(main_colors)},'LineWidth',1.5, ...
'LineStyle','-.');

avg_line=mean(fc_seed_agg,2);
sd_line = std(fc_seed_agg,1,2);
 d4=plot (n_vals,avg_line,strcat(main_colors{1},'-'),'LineWidth',1);
errorbar (n_vals',avg_line,2*sd_line/sqrt(numiter), ...
  strcat(main_colors{1},'.'),'LineWidth',1);


set(d1,'DisplayName','FAQ -  undirected weighted')
    
set(d2,'DisplayName','FAQ -   undirected unweighted');

set(d3,'DisplayName',' chance')
    
set(d4,'DisplayName','FAQ-  undirected no seed');
set(d4,'LineWidth',1);



add_JOFC=1;


if (add_JOFC)
    n_vals_worm=  [20  40  60  80 100 125 150 175 200];
    dim_JOFC = size(JOFC_corr_worm_wt_dice);
     num_it_JOFC= dim_JOFC(2);
    hold on

    avg_line=mean(JOFC_corr_worm_wt_dice,2);
    sd_line = std(JOFC_corr_worm_wt_dice,1,2);
    
    % plot (n_vals,avg_line,'Color',figcolors(i*incr,:),'LineWidth',2)
      e1=plot (n_vals_worm,avg_line,strcat(main_colors{2},'--'),'LineWidth',2.5);
    errorbar (n_vals_worm,avg_line,2*sd_line/sqrt(num_it_JOFC),strcat(main_colors{2},'.'),'LineWidth',1);
    
    hold on


    avg_line=mean(JOFC_corr_worm_wt_dice_unwt(:,:),2);
    sd_line = std(JOFC_corr_worm_wt_dice_unwt(:,:),1,2);
     e2=plot (n_vals_worm,avg_line,strcat(main_colors{2},'-'),'LineWidth',2.5);
    errorbar (n_vals_worm,avg_line,2*sd_line/sqrt(num_it_JOFC),strcat(main_colors{2},'.'),'LineWidth',1);
 
    hold on


    avg_line=mean(JOFC_corr_worm_wt_dice_directed(:,:),2);
    sd_line = std(JOFC_corr_worm_wt_dice_directed(:,:),1,2);
    e3=plot (n_vals_worm,avg_line,strcat(main_colors{2},':'),'LineWidth',2.5);
    errorbar (n_vals_worm,avg_line,2*sd_line/sqrt(num_it_JOFC),strcat(main_colors{2},'.'),'LineWidth',1);
     
    hold on


    avg_line=mean(JOFC_corr_worm_wt_dice_unwt_directed(:,:),2);
    sd_line = std(JOFC_corr_worm_wt_dice_unwt_directed(:,:),1,2);
    % plot (n_vals,avg_line,'Color',figcolors(i*incr,:),'LineWidth',2)
    e4=plot (n_vals_worm,avg_line,strcat(main_colors{2},'-.'),'LineWidth',2.5);

    errorbar (n_vals_worm,avg_line,2*sd_line/sqrt(num_it_JOFC),strcat(main_colors{2},'.'),'LineWidth',1);
 
    set(e1,'DisplayName','JOFC - undirected weighted')
    
set(e2,'DisplayName','JOFC - undirected unweighted');
%,'LineWidth',2
    set(e3,'DisplayName','JOFC - directed weighted')
    
set(e4,'DisplayName','JOFC - directed unweighted');

 
end

add_dir = 1;
if (add_dir)
    
    avg_line=mean(fc_dir_agg(:,:),2);
    sd_line = std(fc_dir_agg(:,:),1,2);
     dir1=plot (n_vals,avg_line,strcat(main_colors{1},':'),'LineWidth',2.5)
        hold on
    errorbar (n_vals,avg_line,2*sd_line/sqrt(select),strcat(main_colors{1},'.'),'LineWidth',1);
 


    avg_line=mean(fc_unwt_dir_agg(:,:),2);
    sd_line = std(fc_unwt_dir_agg(:,:),1,2);
    dir2= plot (n_vals,avg_line,strcat(main_colors{1},'-.'),'LineWidth',2.5)
    errorbar (n_vals,avg_line,2*sd_line/sqrt(select),strcat(main_colors{1},'.'),'LineWidth',1);
    hold on



xlabel('$m$','Interpreter','latex','FontSize',20)
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20)

avg_line=mean(fc_unseed_dir_agg,2);
sd_line = std(fc_unseed_dir_agg,1,2);
 %dir4=plot (n_vals,avg_line,strcat(main_colors{3},'-'),'LineWidth',1);
%errorbar (n_vals',avg_line,2*sd_line/sqrt(numiter), ...
  %strcat(main_colors{3},'.'),'LineWidth',1);


set(dir1,'DisplayName','FAQ -directed weighted ')
    
set(dir2,'DisplayName','FAQ -  directed unweighted ');

%set(dir4,'DisplayName','FAQ-directed no seed ');
%set(dir4,'LineWidth',1);

end


 xlim([-10 253]);


% Create xlabel
xlabel('$m$','Interpreter','latex','FontSize',20);

% Create ylabel
ylabel('$\delta^{(m)}$','Interpreter','latex','FontSize',20);

% Create title
title({'C. elegans'},'FontSize',20);

legend([d1,d2,d4,dir1,dir2,e1,e2,e3,e4,d3,])
% Create legend
legend SHOW;




addpath('~/projects/gurobi/gurobi562/linux64/matlab');
if  ~exist('p') 
    p=0.2;
end
if  ~exist('corr') 
    corr = 0.7;
end

if  ~exist('m') 
    m=10;
end
if  ~exist('n') 
    n=45;
    
end
mplusn = m+n;


[G_1,G_2,true_perm] = generate(p*ones(mplusn),corr,m);

if ~exist('tag') 
    tag='_';
end

model_fname_id = strcat(int2str(mplusn),'verts_edgeprob_pt',int2str(100*p), '_corr_coef_pt',int2str(100*corr),'_',int2str(m),'seeds',tag)
fname_1= GM_exact_model_writer  (G_1,G_2,m,model_fname_id)
fname_2=GM_exact_model_writer_infeasible  (G_1,G_2,m,model_fname_id)

[G_1,G_2,true_perm] = generate(p*ones(mplusn),corr,m);

if ~exist('tag') 
    tag='_';
end

model_fname_id = strcat(int2str(mplusn),'verts_edgeprob_',p, '_corr_coef_pt',num2str(corr),'_',m,'seeds',tag)
fname_1= GM_exact_model_writer  (G_1,G_2,m,model_fname_id)
fname_2=GM_exact_model_writer_infeasible  (G_1,G_2,m,model_fname_id)
load('elegansGraph.mat')

GF=Agap;
GF= GF-triu(GF); GF= (GF>0); GF=GF+GF';
%GE= GE-triu(GE); GE= (GE>0); GE=GE+GE';
N_dims=size(GF)
N_init= N_dims(1)

N = N_init


%GF=Agap;




%GE=GE(1:N,1:N);
GF=GF(1:N,1:N);
%rowsum_E=sum(GE,2);
%unconnected_verts_G1=find(rowsum_E==0)
rowsum_F=sum(GF,2);
unconnected_verts_G2=find(rowsum_F==0)

n_vals=[0 1 5 10 20 50 75 100 150 200 300 350 400 450];
n_vals = n_vals(find(n_vals<N))
num_iter = 1;
corr_match=zeros(length(n_vals),num_iter);
for n_i = 1:length(n_vals)
    
    shuffle= randperm(N-n_vals(n_i));
    shuffle= [1:n_vals(n_i) , shuffle+n_vals(n_i)]
    GE= GF(shuffle,shuffle);
    
    
    for i=1:num_iter
    i
    ordering=randperm(N);
    matching=ConVogHard_rQAP_order(GE,GF,n_vals(n_i),ordering);
    corr_match(n_i,i) =  sum(matching((n_vals(n_i)+1):N)== ...
        shuffle(ordering((n_vals(n_i)+1):N)))
    end
end


pc=mean(corr_match,2)
fc=pc./(N-n_vals')
sd_pc = std(corr_match,0,2)
sd_fc= sd_pc./(N-n_vals')


'Connectome Chem_Chem Finished'
random_chance= 1/(N-n_vals');
plot(n_vals,fc,'r-')
hold on

plot(n_vals,random_chance,'b-.')




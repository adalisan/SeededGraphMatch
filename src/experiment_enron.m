

evalR('setwd("./data")')
evalR('load("AAA-187As-184x184.Rbin")')
evalR('sink("test_R_int.txt")')
evalR('print(str(AAA[[1]]))')
evalR('sink()')

time_stamp_G1 = 131;
time_stamp_G2 = 132;

putRdata('t1', time_stamp_G1)
putRdata('t2', time_stamp_G2)
evalR('GE=AAA[[t1]]')
evalR('GF=AAA[[t2]]')


GE=getRdata('GE');
GF=getRdata('GF');

N_dims=size(GE)
N_init= N_dims(1)

N= N_init%floor(N_init/4)
GE=GE(1:N,1:N);
GF=GF(1:N,1:N);
rowsum_E=sum(GE,2);
find(rowsum_E==0)
rowsum_F=sum(GF,2);
find(rowsum_F==0)
%n_vals=[1 5 10 20 50 100];
n_vals = [0 1 5 10 20 50 60 90 100 140]
%n_vals=[0 90  ];
mismatched_verts_G1=[];
mismatched_verts_G2=[];
seeds=zeros(num_iter,100);
num_iter = 60;
corr_match=zeros(length(n_vals),num_iter);
for n_i = 1:length(n_vals)
    for i=1:num_iter
    i
    ordering = randperm(N);
    ordering(1:n_vals(n_i)) = sort(ordering(1:n_vals(n_i)));
    ordering(n_vals(n_i)+1:N) = sort(ordering(n_vals(n_i)+1:N));
    matching=ConVogHard_rQAP_order(GE,GF,n_vals(n_i),ordering);
    corr_match(n_i,i) =  sum(matching(n_vals(n_i)+1:N)==ordering(n_vals(n_i)+1:N));
    if (n_vals(n_i)==100)
        mismatched_ind=find(matching(n_vals(n_i)+1:N)~=ordering(n_vals(n_i)+1:N));
        mismatched_verts_G1=[mismatched_verts_G1 NaN matching(n_vals(n_i)+mismatched_ind)]
      
        mismatched_verts_G2=[mismatched_verts_G2 NaN ordering(n_vals(n_i)+mismatched_ind)]
        seeds(i,:) = ordering(1:n_vals(n_i));
    end
    end
end

pc=mean(corr_match,2)
sd_pc = std(corr_match,0,2)

fc= pc./(N-n_vals')
sd_fc= sd_pc./(N-n_vals')

'Enron Finished'
random_chance= 1./(N-n_vals');
figure
hold on

plot(n_vals,random_chance,'b-.')




if time_stamp_G1==130 && time_stamp_G2==131
corr_match_run0=corr_match;
pc_run0 = pc;
sd_pc_run0 = sd_pc;
fc_run0 = fc;
sd_fc_run0 = sd_fc;
mismatched_verts_G1_run0 = mismatched_verts_G1;
mismatched_verts_G2_run0 = mismatched_verts_G2;
elseif time_stamp_G1==130 && time_stamp_G2==131

corr_match_run1=corr_match;
pc_run1 = pc;
sd_pc_run1 = sd_pc;
fc_run1 = fc;
sd_fc_run1 = sd_fc;
mismatched_verts_G1_run1 = mismatched_verts_G1;
mismatched_verts_G2_run1 = mismatched_verts_G2;


else %(time_G1=130 and time_G2=132
corr_match_run2=corr_match;
pc_run2 = pc;
sd_pc_run2 = sd_pc;
fc_run2 = fc;
sd_fc_run2 = sd_fc;
mismatched_verts_G1_run2 = mismatched_verts_G1;
mismatched_verts_G2_run2 = mismatched_verts_G2;

end

figure
errorbar(n_vals,fc_run0,2*sd_fc_run0/sqrt(num_iter),'r-')
hold on 
errorbar(n_vals,fc_run1,2*sd_fc_run1/sqrt(num_iter),'b-')
errorbar(n_vals,fc_run2,2*sd_fc_run2/sqrt(num_iter),'g-')

plot(n_vals,random_chance,'g:')

title('Enron graph matching, graphs at t=130, 131, 132')
legend('130&131','131&132','130&132')





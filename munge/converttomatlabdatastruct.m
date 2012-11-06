
evalR('load("./data/wiki.Rdata")')


GE=getRdata('GE');
GF=getRdata('GF');

TE=getRdata('TE');
TF=getRdata('TF');





GM=getRdata('GM');
TM=getRdata('TM');
label= getRdata('label');


evalR('load("./data/AAA-187As-184x184.Rbin")')

AAA=zeros([184 184 187]);
for i=1:187
    
time_stamp_G1 = i;

putRdata('time_stamp_G1',time_stamp_G1);

evalR('Gi=AAA[[time_stamp_G1]]')

Gi=getRdata('Gi');
AAA(:,:,i) =Gi;
end


MAP=zeros(6,5,5);
kkk=99;
EVA='MultiEva';
load('BreaH100.mat')

len=[16,32,48,64,128,256];

I_tr=double(I_tr);
I_re=double(I_re);
I_te=double(I_te);

L_tr=double(L_tr);
L_re=double(L_re);
L_te=double(L_te);

for i=1:6
   i
fs=len(i);

MAP(i,1,:)=LBSE(kkk,fs,EVA,I_tr,I_te,I_re,L_tr,L_te,L_re);
MAP(i,2,:)=LBSE(kkk,fs,EVA,I_tr,I_te,I_re,L_tr,L_te,L_re);
MAP(i,3,:)=LBSE(kkk,fs,EVA,I_tr,I_te,I_re,L_tr,L_te,L_re);
MAP(i,4,:)=LBSE(kkk,fs,EVA,I_tr,I_te,I_re,L_tr,L_te,L_re);
MAP(i,5,:)=LBSE(kkk,fs,EVA,I_tr,I_te,I_re,L_tr,L_te,L_re);

end
save('MultiEVA_BreaH100.mat','MAP');

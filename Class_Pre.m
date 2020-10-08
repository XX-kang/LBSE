function[nb] = Class_Pre(kkk,B_te,B_re,L_query,L_re)

Y_te=L_query;
Y_re=L_re;

[dig,c]=size(L_query);

B=B_te;
H=B_re;

comspace= H;
comspace1=B;


nb=zeros(dig,c);
%     tic;
for i=1:dig
%     i
    a=foundHammingCF(comspace1(i,:),comspace,Y_te,Y_re);
    a=a(1:kkk,:);

    sum_class_ori=sum(a,1)/sum(sum(a));
%     save('sum_class_ori.mat');
    topK=sort(sum_class_ori,'descend');topK=topK(1);
    sum_class_all=(sum_class_ori>=topK);
    nb(i,:)=sum_class_all;
end

end



function[num] = foundHammingCF(tbf,comspace,Y_te,Y_re)
% load('cifar10.mat')
cat=Y_te;
cat2=Y_re;
[~,c]=size(Y_te);
[numtrain,~]=size(comspace);
% numtrain
tbf=repmat(tbf,numtrain,1);
% size(tbf)
% size(comspace)
dis=sum(bitxor(tbf,comspace),2);
% size(dis)
% size(cat2)
dis2=horzcat(dis,cat2);
dis=sortrows(dis2,1);

num=dis(:,2:c+1);
end
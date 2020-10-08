
  function[mAP] =LBSE(kkk,fs,EVA,I_tr,I_query,I_re,L_tr,L_query,L_re)
%
I_tr=I_tr';
L_tr=L_tr';
Y=L_tr;

len = fs;

[c,n] = size(Y);

[doI,~] = size(I_tr);

lamda =1e-4;
alpha=0.1;
beta =20;
gamma=1;
W=rand(len,c)-0.5;

B=sign(rand(len,n)-0.5);
Z=rand(len,n)-0.5;


dl=10; 
S=L_tr'*L_tr>0;
S=2*S-1;

P_temp=pinv(I_tr*I_tr'+lamda*diag(ones(doI,1)))*(I_tr);
Pv=P_temp*B';
loss=zeros(dl,1);
tic;
for loops=1:dl
    %     loops
    Pv=pinv(I_tr*I_tr'+beta*diag(ones(doI,1)))*(I_tr*B');
    B=sign(1e-0*S*Z'+lamda*I_tr'*Pv+alpha*Z'+gamma*Y'*W');
    B=B';
    
    if fs>=c
        [U, ~, V] = svd(Y*B');
        W =V(:,1:c)*U;
    elseif fs<c
        [U, ~, V] = svd(B*Y');
        W =U*V(:,1:fs)';
    end


    Z = S*B'+beta*B';
    Temp = Z'*Z-1/n*(Z'*ones(n,1)*(ones(1,n)*Z));
    [~,Lmd,QQ] = svd(Temp); 
    idx = (diag(Lmd)>1e-6);
    Q = QQ(:,idx); Q_ = orth(QQ(:,~idx));
    P = (Z-1/n*ones(n,1)*(ones(1,n)*Z)) *  (Q / (sqrt(Lmd(idx,idx))));
    P_ = orth(randn(n,fs-length(find(idx==1))));
    Z = sqrt(n)*[P P_]*[Q Q_]';
    Z=Z';
end

mAP=toc;

B_te=I_query*Pv>0;
B_re=I_re*Pv>0;


if  strcmp(EVA,'MultiEva')

     mAP=Class_Pre(kkk,B_te,B_re,L_query,L_re);
    [num,~]=size(L_query);
    EM=mAP-L_query;
 
    TP=zeros(1,1);
    TN=zeros(1,1);
    FP=zeros(1,1);
    FN=zeros(1,1);
    
    for i=1:num
        if sum(abs(EM(i,:)))==0&&L_query(i,5)==0&&L_query(i,6)==0&&L_query(i,7)==0&&L_query(i,8)==0
            TP=TP+1;
        end
        if sum(abs(EM(i,:)))==0&&L_query(i,1)==0&&L_query(i,2)==0&&L_query(i,3)==0&&L_query(i,4)==0
            TN=TN+1;
        end
        if sum(abs(EM(i,:)))==2&&L_query(i,5)==0&&L_query(i,6)==0&&L_query(i,7)==0&&L_query(i,8)==0
            FP=FP+1;
        end
        
        if sum(abs(EM(i,:)))==2&&L_query(i,1)==0&&L_query(i,2)==0&&L_query(i,3)==0&&L_query(i,4)==0
            FN=FN+1;
        end
    end

    EvaMetric=zeros(5,1);
    EvaMetric(1,1)=(TP+TN)/(TP+TN+FP+FN);
    EvaMetric(2,1)=(TP)/(TP+FP);
    EvaMetric(3,1)=(TP)/(TP+FN);
    EvaMetric(4,1)=2*(EvaMetric(2,1)*EvaMetric(3,1))/(EvaMetric(2,1)+EvaMetric(3,1));
    EvaMetric(5,1)=(TP*TN-FP*FN)/sqrt((TP+FP+1e-8)*(TP+FN+1e-8)*(TN+FP+1e-8)*(TN+FN+1e-8));
    mAP=EvaMetric;
end
end
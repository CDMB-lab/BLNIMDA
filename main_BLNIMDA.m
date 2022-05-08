clc
clear
tic

dSemanticArray_1=textread('./Data/disease semantic similarity1.txt');
dSemanticArray_2=textread('./Data/disease semantic similarity2.txt');
rFunctionalArray=textread('./Data/miRNA functional similarity.txt');
dWeightArray=textread('./Data/disease semantic similarity weight.txt');
rWeightArray=textread('./Data/miRNA functional similarity weight.txt');
KnownAssociation_ori=xlsread('./Data/miRNA-disease.xlsx');

dSemanticArray=(dSemanticArray_1+dSemanticArray_2)/2;

m=max(KnownAssociation_ori(:,1)); %495
n=max(KnownAssociation_ori(:,2)); %383
Y=zeros(m,n);    %创建m*n的全零数组

for i=1:length(KnownAssociation_ori)
    Y(KnownAssociation_ori(i,1),KnownAssociation_ori(i,2))=1; %初始关系邻接矩阵
end
%经过以上处理得到了miRNA-disease关联矩阵；疾病语义相似性平均矩阵；miRNA功能相似性矩阵

[score_ori]=BLNIMDA(rFunctionalArray,dWeightArray,rWeightArray,dSemanticArray,Y);
pp=length(KnownAssociation_ori);
L=Y;
index=find(1==Y);
score_0=score_ori;
score_0(index(:))=0;
[x,y]=size(score_0);
score_1=zeros(x,y);


% %% global loocv
% 
% for i=1:pp
%     i
%     L(index(i))=0;
%     
%     [finT]=WBNPMD(rFunctionalArray,dWeightArray,rWeightArray,dSemanticArray,L);
%     score_1(index(i))=finT(index(i));
% 
%     L=Y;
% end
% 
% score_final=score_0+score_1;
% auc=roc_1(score_final(:),Y(:),'green');

% auc(k)=roc_1(score_final(:),Y(:),'red');



%% 5-fold cross validation
auc=zeros(1,50);
for j=1:50
    j
    score_this=zeros(x,y);
    L=Y;
    indices = crossvalind('Kfold', pp, 5);
    for i=1:5
    i
    index1=find(i==indices);
    L(index(index1))=0;
    
    score=BLNIMDA(rFunctionalArray,dWeightArray,rWeightArray,dSemanticArray,L);
    score_this(index(index1))=score(index(index1));
    
    L=Y;
    end
    
    score_1=score_1+score_this;
    score_thisfin=score_0+score_this;
    auc=roc_1(score_thisfin(:),Y(:),'red');
    auc1=roc_1(score_thisfin(:),Y(:),'red'); 

%     score_final=score_0+score_1./100;
%     auc_final=roc_1(score_final(:),Y(:),'red');
end

toc
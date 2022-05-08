function [S_F]=BLNIMDA(rFunctionalArray,dWeightArray,rWeightArray,dSemanticArray,interaction)

%% 初始化miRNA-疾病关系矩阵
m=495;
n=383;
S_r=zeros(m,m);
S_d=zeros(n,n);

%% 高斯核相互作用
[GaussR,GaussD] = gaussiansimilarity(interaction,m,n);

for i=1:n
    for j=1:n
        if (dWeightArray(i,j) == 1) S_d(i,j)=dSemanticArray(i,j);
        else if (dWeightArray(i,j) == 0) S_d(i,j)=GaussD(i,j);
            end
        end
    end
end


for i=1:m
    for j=1:m
        if (rWeightArray(i,j) == 1) S_r(i,j)=rFunctionalArray(i,j);
        else if (rWeightArray(i,j) == 0) S_r(i,j)=GaussR(i,j);
            end
        end
    end
end

%% 构建初始信息
H_r=zeros(m,m);
H_d=zeros(n,n);
 for i=1:m
     for j=1:m
         if S_r(i,j)>0.02
             H_r(i,j)=S_r(i,j);
         else
             H_r(i,j)=0;
         end
     end
 end
 
 for i=1:n
     for j=1:n
         if S_d(i,j)>0.02
             H_d(i,j)=S_d(i,j);
         else
             H_d(i,j)=0;
         end
     end
 end
% y1=S_r*interaction./sum(S_r,2);%.*interaction; %miRNA权重计算
 y1=H_r*interaction./sum(S_r*interaction,2);  %miRNA-disease初始信息矩阵
% y2=S_d*test./sum(S_d,2); %疾病权重计算
y2=H_d*interaction'./sum(S_d*interaction',2); %disease-miRNA初始信息矩阵

S_rr=zeros(m,m);
S_dd=zeros(n,n);
for i=1:m
    for j=1:m
        if S_r(i,j)~=1
            S_rr(i,j)=S_r(i,j);
        end
    end
end
for i=1:n
    for j=1:n
        if S_d(i,j)~=1
            S_dd(i,j)=S_d(i,j);
        end
    end
end

% function[max_m,L]=max_1(matrix)
  
%     [o,p]=size(matrix);
%     max_m=zeros(o,p);
%     L=zeros(o,p);
%     for a=1:o
%       for b=1:p
%         if matrix(a,b)>max(max_m)
%             max_m(a)=matrix(a,b);
%             L(a)=b;
%         end
%       end
%     end
% end
nzero1=sum(H_r~=0,2);
nzero2=sum(H_d~=0,2);
%% 第一次 miRNA
W1=zeros(m,n);
[max_m,L]=max(S_rr,[],2);
for i=1:m
    for j=1:n
        if interaction(i,j)~=0
            W1(i,j)=exp(interaction(i,j));     
        elseif interaction(i,j)==0 && interaction(L(i),j)~=0
            W1(i,j)=exp(max_m(i));
        else
            W1(i,j)=exp(sum(H_r(i,[]))/383);
        end
    end
end


sumrow1=zeros(m,1);
sumrow1=sum(W1,2);
Wfin1=W1./repmat(sumrow1,1,size(W1,2));%第i行每个元素除以第i行元素之和
S_D=y1.*Wfin1;    



%% 第二次 疾病 （邻接矩阵383 by 495)

W2=zeros(m,n);
[max_n,M]=max(S_dd,[],2);
for i=1:m
    for j=1:n
        if interaction(i,j)~=0
            W2(i,j)=exp(1);
        elseif interaction(i,j)==0 && interaction(i,M(j))~=0
            W2(i,j)=exp(max_n(j));
        else
            W2(i,j)=exp(sum(H_d([],j)/495));
        end
    end
end
sumrow2=zeros(n,1);
sumrow2=sum(W2,2);
Wfin2=W2./repmat(sumrow2,1,size(W2,2));
S_M=y2.*Wfin2';   

%% 计算最终得分
S_F=(S_D+S_M')/2;

% S_F=S_M;
% S_F=S_D';

toc
end

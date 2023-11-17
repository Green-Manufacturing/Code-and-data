%Dagum基尼系数分解主函数文件
function [G_Total,G_W,G_sub,G_nb,G_jh,G_t,G_test]=Dagum_Gini(varargin)
%% 输入变量和输出变量说明
%
% 输入变量：假定有子群：A，B,C，...，首先，按各子群均值由大到小排序，
%         如:B的均值>A的均值>C的均值>...,则输入命令为：Dagum_Gini(B,A,C,...)
%         注：所有子群变量均以列向量形式输入。
%
% 输出变量：
%         G_Total：   总体基尼系数
%         G_W：       子群内差异贡献
%         G_sub：     子群内基尼系数
%         G_nb：      子群间差异贡献
%         G_jh：      子群间基尼系数
%         G_t：       超变密度贡献
%         G_test：    检验G_W+G_nb+G_t是否等于G_Total

%% 检查变量输入格式
n_sub_group=nargin;                     %输入变量数（子群数）
all_sub_group=varargin;                 %输入的所有子群（以元胞数组形式存储）

for i_group=1:n_sub_group
    [~,check_size]=size(all_sub_group{i_group});
    if check_size~=1
        error('输入格式错误：输入的各变量必须为列向量！')
    end
end

for i_group=1:n_sub_group-1
    if mean(all_sub_group{i_group})>=mean(all_sub_group{i_group+1})
        %disp('变量输入格式正确！')
    else
        error('输入格式错误：请按子群均值从大到小依次输入变量！')
    end
end

%% 计算整体和子群的样本数、均值，子群样本数占整体比重等
Total_group=cell2mat(all_sub_group');    %将子群合并成总体
n_T=length(Total_group);                 %计算总体样本数
m_T=mean(Total_group);                   %计算总体均值
        
% 预分配内存
n_sample_sub_group=zeros(n_sub_group,1); %存储各子群样本数量
m_sub_group=n_sample_sub_group;          %存储各子群均值
P=n_sample_sub_group;                    %存储各子群样本数占总体样本数比重
S=n_sample_sub_group;                    %存储各子群元素之和占总体元素之和的比重

for i_group=1:n_sub_group
    n_sample_sub_group(i_group)=length(all_sub_group{i_group});%b{i}是矩阵（列向量）
    m_sub_group(i_group)=mean(all_sub_group{i_group});
    P(i_group)=length(all_sub_group{i_group})/n_T;
    S(i_group)=length(all_sub_group{i_group})*mean(all_sub_group{i_group})/(n_T*m_T);
end

%% 载入基尼系数求解子函数
my_GINI=GINI_COMPUTE;

%% 计算各子群基尼系数
G_sub=zeros(n_sub_group,1);%预分配内存
for i_group=1:n_sub_group
    %循环计算各子群的基尼系数，并将结果存储在G_sub矩阵中
    G_sub(i_group)=my_GINI.GINI_1(all_sub_group{i_group},all_sub_group{i_group});
end

%% 使用for循环+矩阵上三角化计算G_jh、D_jh和P_S
%预分配内存
G_jh_T=zeros(n_sub_group,n_sub_group);
D_jh_T=zeros(n_sub_group,n_sub_group);
P_S_T=zeros(n_sub_group,n_sub_group);

for i_group=1:n_sub_group
    for j_group=1:n_sub_group
        %循环计算子群间基尼系数矩阵
        G_jh_T(i_group,j_group)=my_GINI.GINI_1(all_sub_group{i_group},all_sub_group{j_group});
        %循环计算子群间的相对影响矩阵
        D_jh_T(i_group,j_group)=my_GINI.D_jh(all_sub_group{i_group},all_sub_group{j_group});
        %循环计算权重矩阵
        P_S_T(i_group,j_group)=P(i_group)*S(j_group);
    end
end

% G_jh、D_jh和P_S的矩阵上三角化
G_jh=triu(G_jh_T,1);
D_jh=triu(D_jh_T,1); 
D_jh_1=1-D_jh;
D_jh_1=triu(D_jh_1,1);
P_S=triu(P_S_T,1)+triu(P_S_T',1); %（重要的一步）原矩阵与转置矩阵相加

%% 计算基尼系数分解结果
G_Total=my_GINI.GINI_1(Total_group,Total_group);   %计算总体基尼系数
G_nb=sum(sum(P_S.*D_jh.*G_jh));        %计算子群间差异贡献额
G_t=sum(sum(P_S.*D_jh_1.*G_jh));       %计算超变密度  
G_W=sum(P.*S.*G_sub);                  %计算组内基尼系数
G_test=G_W+G_nb+G_t;                   %检验总体基尼系数（G_Total）是否等于G_W+G_nb+G_t

%% 显示计算结果
disp('%%%%%%%%%%%%%%%%Dagum基尼系数分解结果%%%%%%%%%%%%%%%')
disp(' ')
disp(['总体基尼系数（G_T）: ',num2str(G_Total)])
disp(' ')
disp(['子群内差异贡献（G_W）: ',num2str(G_W)])
disp('------------各子群基尼系数------------')
for i_group=1:n_sub_group
    disp(['子群内基尼系数','（G_sub(',num2str(i_group),')):', num2str(G_sub(i_group))])
end
disp('-----------------------------------')
disp(' ')
disp([char('子群间差异贡献（G_nb）: ') num2str(G_nb)])

disp('------------子群间基尼系数------------')
for i_group=1:n_sub_group
    for j_group=1:n_sub_group
        if j_group>i_group
        disp(['子群间基尼系数','（G_jh_',num2str(i_group),'_',num2str(j_group),'):', num2str(G_jh(i_group,j_group))])
        end
    end
end 

disp('-----------------------------------')
disp(' ')
disp([char('超变密度贡献（G_t）: ') num2str(G_t)])
disp(' ')

disp('===========各部分贡献率(%)===========')
disp(['子群内差异贡献率:',num2str(G_W/G_Total*100)])
disp(['子群间差异贡献率:',num2str(G_nb/G_Total*100)])
disp(['超变密度贡献率:',num2str(G_t/G_Total*100)])
disp('===================================')
disp(' ')

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

end
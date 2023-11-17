%% 计算整体、子群内、子群间基尼系数和子群间相对影响

function GINI_MAIN=GINI_COMPUTE

%将一组功能独立的函数写在同一个函数文件中
%主函数GINI_MAIN(输出)=GINI_COMPUTE(输入)
%注1：运行该函数时可使用：my_GINI=GINI_COMPUTE
%注2：对应的子函数输入形式为：my_GINI.GINI_1

GINI_MAIN.GINI_1=@GINI_1;%计算整体、子群内和子群间间基尼系数
GINI_MAIN.D_jh=@D_jh;    %计算子群间相对影响（D_ij）
end

%% 子函数1：计算整体、子群内和子群间基尼系数

function GINI1=GINI_1(xx,yy)
n1=length(xx);
n2=length(yy);
m1=mean(xx);
m2=mean(yy);
D_ij_2=zeros(n1,n2);

for ii=1:n1
    for jj=1:n2
        D_ij_2(ii,jj)=xx(ii)-yy(jj);%列向量元素两两相减的差值矩阵
    end
end

sum_d2=sum(sum(abs(D_ij_2))); %对矩阵所有元素先取绝对值，然后求和 
GINI1=sum_d2/(n1*n2*(m1+m2)); %计算总体、组内和组间基尼系数
end

%% 子函数2:计算子群j和子群h变量的相对影响D_jh

function D_jh1=D_jh(xx,yy)
n1=length(xx);
n2=length(yy);
D_ij_3=zeros(n1,n2);

for ii=1:n1
    for jj=1:n2
        D_ij_3(ii,jj)=xx(ii)-yy(jj);
    end
end
%从矩阵中找出数值大于0的元素 
E11=sort(D_ij_3(find(D_ij_3>0)));
%子群j和子群h间差值(算术平均)
d_jh=sum(E11)/(n1*n2);   

%从矩阵中找出数值小于0的元素并对他们取绝对值
E12=abs(sort(D_ij_3(find(D_ij_3<0))));
%超变一阶矩（算术平均）
p_jh=sum(E12)/(n1*n2);

%j区域和h区域变量的相对影响
D_jh1=(d_jh-p_jh)/(d_jh+p_jh);
%注1：若子群j的样本均值和子群h的样本均值相同，则D_jh1=0
%注2：若子群j的所有元素均大于子群子群h的元素，则D_jh1=1
%注3：其他情况下，0<D_jh1<1
end

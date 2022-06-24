clc,clear
x=[1, 3, 5, 9, 6]'; % 将3等效到9上
[R,y]=givens(x,4,2) % 注意3的下标为2，9的下标为4
function [R,y]=givens(x,i,j)
% 求解标准正交的Given变换矩阵R，使用Rx=y，其中y(j)=0，y(i)=sqrt(x(i)^2+x(j)^2)
% 参数说明
% x：需要进行Givens变换的列向量
% i：变为sqrt(x(i)^2+x(j)^2)的元素下标
% j：变为0的元素的下标
% R：Givens变换矩阵
% y：Givens变换结果
% 实例说明
% x=[1 3 5 9 6]'; % 将3等效到9上
% [R,y]=givens(x,4,2) % 注意3的下标为2，9的下标为4
% R*x-y % 验证Rx=y
% 关于Givens变换说明
% 1.Givens矩阵是标准正交矩阵，也叫平面旋转矩阵，它是通过坐标旋转的原理将元素j的数值等效到元素i上
% 2.Givens变换每次只能将一个元素变为0，而Householder变换则一次可以将任意个元素变为0
% 3.Givens变换常用于将矩阵A变为对角阵
% see also http://www.matlabsky
    xi=x(i);
    xj=x(j);
    r=sqrt(xi^2+xj^2);
    cost=xi/r;
    sint=xj/r;
    R=eye(length(x));
    R(i,i)=cost;
    R(i,j)=sint;
    R(j,i)=-sint;
    R(j,j)=cost;
    y=x(:);
    y([i,j])=[r,0];
end
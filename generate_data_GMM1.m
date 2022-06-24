clear, clc;
model_class = 2;
dim = 3;
% 期望值
m = [0, 0, 0;
    0, 0, 0];
% 协方差阵
s(:, :, 1) = [0.02, 0, 0;
              0, 0.02, 0;
              0, 0, 0.02];
s(:, :, 2) = [0.02, 0, 0;
              0, 0.02, 0;
              0, 0, 0.02];
num = [10, 3];
data = generate_data_GMM(dim, model_class, m, s, num);

% 提出每组数据
data1 = data((data(:,4) == 1), (1:3));
data2 = data((data(:,4) == 2), (1:3));

% 绘图分析（以第二组数据为例）
figure(1)
plot(data2(:, 1), data2(:, 2), '.k');
axis equal;
figure(2)
plot(data2(:, 2), data2(:, 3), '.k');
axis equal;
figure(3)
plot(data2(:, 1), data2(:, 3), '.k');
axis equal;


function data = generate_data_GMM(dim, model_class, m, s, num)
    %   生成多组多维高斯分布数据
    %   dim为高斯分布的维数
    %   model_class为生成的数据组数
    %   m为高斯分布的期望值，大小为model_class*dim
    %   s为高斯分布的协方差阵，大小为dim*dim*model_class
    %   num为各组高斯分布的数据量，其大小为model_class*1
    %   返回值data为生成的数据，其大小为num*(dim+1)
    %   前dim列为高斯分布数据，第(dim+1)列为组的类别编号
    %   该函数用于生成多组高斯分布数据，可为聚类算法提供数据
    data = [];
    for i = 1 : model_class
        data1 = ones(num(i), dim + 1);
        data1(:, (1 : dim)) = mvnrnd(m(i, :), s(:, :, i), num(i));
        for j = 1 : num(i)
            data1(j, dim + 1) = i;
        end
        data = [data; data1];
    end
end
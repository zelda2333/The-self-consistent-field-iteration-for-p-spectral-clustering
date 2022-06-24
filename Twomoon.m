clc,clear

% model_class = 1;
% dim = 10;
% % 期望值
% mu = zeros(model_class,dim);
% % 协方差阵
% sigma(:, :, 1) = 0.02*diag(ones(dim,1));
% 
% num = [20];
% data = generate_data_GMM(dim, model_class, mu, sigma, num);
% 
% % 提出每组数据
% data1 = data((data(:,dim+1) == 1), (1:dim));
% for i = 1:dim-1
%         x_((i-1)*num+1:i*num) = data1(:, i);
%         y_((i-1)*num+1:i*num) = data1(:, i+1);
% end


aplha=0:pi/200:pi;
r=1;
x=r*cos(aplha);
y=r*sin(aplha);
x_=randn(size(x))./10+x;
y_=randn(size(y))./10+y;


aplha2=pi:pi/200:2*pi;
x2=1+r*cos(aplha2);
y2=0.5+r*sin(aplha2);

x2_=randn(size(x2))./10+x2;
y2_=randn(size(y2))./10+y2;



plot(x,y,'ko');
axis equal
hold on

plot(x_,y_,'ko');
hold on
plot(x2,y2,'ko');
hold on

plot(x2_,y2_,'ko');


% x=0:0.05:2*pi;
% 
% c=cos(x);
% 
% s=sin(x);
% 
% ss=s(find(s.*(x>=0&x<=3)));
% 
% cc=c(find(c.*(x>=1.5&x<=4.5)));
% 
% xc=x(find(x.*(x>=1.5&x<=4.5)));
% 
% xs=x(find(x.*(x>=0&x<=3)));
% 
% % s1=randn(size(ss))./10+ss;
% % 
% % c1=randn(size(cc))./10+cc;
% 
% plot(xc,cc,'ro');
% 
% hold on
% 
% plot(xs,ss,'bd');
% % axis equal
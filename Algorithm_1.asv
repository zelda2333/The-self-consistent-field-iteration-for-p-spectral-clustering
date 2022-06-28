%  REFERENCES: Upadhyaya P, Jarlbering E, Tudisco F. The self-consistent 
% field iteration for p-spectral clustering.2021 
clc,clear
load('Twomoons.mat'); data = twomoons;
num_cluster = 2;
[n, ~] = size(data);
x = data;
%% Construct a similarity matrix based on k-nearest neighbor.
d = zeros(n); K = zeros(n); S = zeros(n);
for i = 1:n
    for j = i+1:n
        % Euclidean distance of each points.
        d(i,j) = norm(x(i,:) - x(j,:)); d(j,i) = d(i,j);
    end
    % find undirected k-nearest neighbors for each of the points.
    [sort_d, index] = sort(d(i,:));
    index(index == i) = [];
    % Marks the k-nearest neighbor of each element.
    K(i,index(1:10)) = 1;
    K(index(1:10),i) = 1;
    % measure similarity based on k-nearest neighbors.
    for j = i+1:n
        if K(i,j) == 1
            % Gaussian kernel function
            S(i,j) = exp(-2*norm(x(i,:) - x(j,:))^2 / (sort_d(2)^2));
            S(j,i) = S(i,j);
        end
    end 
end

% total edge m
edge = sum(K,'all') / 2;
% relationship between edges and corresponding nodes encoded in the graph incidence matrix B
B = zeros([edge,n]);
% D_w  is a diagonal matrixthe kth diagonal entry denoting the weight of edge k
D_w =zeros(edge);
m = 1;
for i = 1:n
    for j = i+1:n      
       if K(i,j) == 1
        B(m,i) = -1;
        B(m,j) = 1;        
        D_w(m,m) = S(i,j);
        m = m+1;   
       end       
    end
end

%% Compute the unnormalized 2-Laplacian the secound smallest veter to get v0.
D = sum(S,2); D = diag(D);
L = D - S;
% [v0, ~] = eigs(L, 2, 'smallestabs');  
% \L is singular matrix,
[v0, va] = eig(L);
%% Iteratively solve eigenpair  SCF
V = v0(:,2);
res = 1;
k = 1;
a = 100;
p = 1.9;

% [V,k] = Algorithm1(B,D_w,V,a,res,p,k);
% [V,k] = Algorithm2(B,D_w,V,a,res,p,k);
[V,k] = Algorithm3(B,D_w,V,a,res,p,k);


%% Cluster the points with k-means algorithm.
C = kmeans(V(:,k), num_cluster);

%% Compute ratio cut (RCut) and ratio Cheeger cut (RCC) 
num_class_1 = size(data(C==1),1);
class_1_point = find(C==1);

num_class_2 = n-num_class_1;
class_2_point = find(C==2);

% Compute the weight of the cut graph between class_1 and class_2
w = 0;
for i = 1:num_class_1
    for j = 1:num_class_2
        a = class_1_point(i);
        b = class_2_point(j);
        w = w + S(a,b);
    end
end

RCut = (1/num_class_1 + 1/num_class_2) * w
RCC = w / min(num_class_1,num_class_2)

%% Plot
figure(1)
clf
plot(data(C==1,1),data(C==1,2),'r.','MarkerSize',10,'linewidth',5); 
hold on
plot(data(C==2,1),data(C==2,2),'b.','MarkerSize',10,'linewidth',5);
title('Spectral clustering');
set(gca,'FontSize',18)
set(gca,'FontName','times')
set(gcf,'color',[1,1,1]);
% end



function [V,k] = Algorithm1(B,D_w,V,a,res,p,k)
    while res > 1/a
        N = B.'*D_w*(diag(abs(B*V(:,k))))^(p-2)*B; 
        R = (diag(abs(V(:,k))))^(p-2);
        M = R\N;
    %     [vecter, value] = eigs(M, 2, 'smallestabs');  M is singular matrix,
        [vecter, value] = eig(M);
        V(:,k+1) = vecter(:,2); 
        lambda = value(2,2);
        k = k+1; 
        res = norm((N - lambda*R)*V(:,k));
    end
end

function [V,k] = Algorithm2(B,D_w,V,a,res,p,k)
    while res > 1/a
        N = B.'*D_w*(diag(abs(B*V(:,k))))^(p-2)*B;
        R = (diag(abs(V(:,k))))^(p-2);
        M = R^(-1/2)*N*R^(-1/2);
    %     [vecter, value] = eigs(M, 2, 'smallestabs');  M is singular matrix,
        [vecter, value] = eig(M);
        V(:,k+1) = R^(-1/2)*vecter(:,2);
        lambda = value(2,2);    
        k = k+1;  
        res = norm((N - lambda*R)*V(:,k));
    end
end

function [V,k] = Algorithm3(B,D_w,V,a,res,p,k)
    while res > 1/a
        sfBx = 2*log(1+exp(-a*B*V(:,k)))/a;
        sfx = 2*log(1+exp(-a*V(:,k)))/a;
        N = B.'*D_w*(diag(abs(sfBx)))^(p-2)*B;
        R = (diag(abs(sfx)))^(p-2);
        M = R^(-1/2)*N*R^(-1/2);
    %     [vecter, value] = eigs(M, 2, 'smallestabs');  M is singular matrix,
        [vecter, value] = eig(M);
        V(:,k+1) = R^(-1/2)*vecter(:,2);
        lambda = value(2,2);    
        k = k+1;  
        res = norm((N - lambda*R)*V(:,k));
    end
end
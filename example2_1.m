clear; clc; close all

x1 = [0,0]; x2 = [0,10]; x3 = [10,0]; x4 = [10,10];
X = [x1;x2;x3;x4]'; x = [2,3]'; L = size(X,2); 
d = (sqrt(sum((x*ones(1,L)-X).^2,1))).';

dB = 30; 
sigma2 = d.^2/10^(dB/10);  
r = d + randn(L,1).*sqrt(sigma2);
iter = 30;

% NLS Newton-Raphson  %

x=[3,2]'; %initial guess value
for i = 1:iter
    H = hessian_nls(X,x,r);
    g = grad_nls(X,x,r);
    x = x-inv(H)*g;
    x_nr(i,:) = x;
end

% NLS Gauss-Newton algorithm  %
x=[3,2]';
for i = 1:iter
    G = jacob(X,x);
    f_TOA = sqrt(sum((ones(L,1)*x'-X').^2,2));
    x = x+inv(G'*G)*G'*(r-f_TOA);
    x_gn(i,:) = x;
end

% NLS steepest descent algorithm  %
x=[3,2]';
mu= 0.1;
for i = 1:iter
    g = grad_nls(X,x,r);
    x = x - mu*g;
    x_sd(i,:) = x;
end

% == displaying results == %
iter_no = 1:iter;
figure
plot(iter_no, x_nr(:,1), 'k.', iter_no, x_gn(:,1), 'ko', iter_no, x_sd(:,1), 'k+');
legend('Newton-Raphson','Gauss-Netwon','steepest descent');
xlabel('Number of Iteration'); ylabel('Estimate of x'); 

figure
plot(iter_no, x_nr(:,2), 'k.', iter_no, x_gn(:,2), 'ko', iter_no, x_sd(:,2), 'k+');
legend('Newton-Raphson','Gauss-Netwon','steepest descent');
xlabel('Number of Iteration'); ylabel('Estimate of x'); 




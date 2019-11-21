function g = grad_nls(X,x,r)
L = size(X,2); t1 = 0; t2 = 0;
ds = sum((x*ones(1,L)-X).^2,1);
ds = ds';
for i=1:L
    t1 = t1 + (r(i)-ds(i)^(0.5))*(x(1)-X(1,i))/ds(i)^(0.5);
    t2 = t2 + (r(i)-ds(i)^(0.5))*(x(2)-X(2,i))/ds(i)^(0.5);
end
g = -2*[t1; t2];
end
function G = jacob(X,x)
[dim,L] = size(X); 
f_TOA = sqrt(sum((ones(L,1)*x'-X').^2,2));
G = (ones(L,1)*x' - X')./(f_TOA*ones(1,dim));
end
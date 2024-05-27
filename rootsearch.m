function [x1,x2] = rootsearch(func,a,b,dx)
% source: Kiusalaas, J., 2010. Numerical Methods in Engineering with MATLAB, Cambridge university press
% Incremental search for a root of f(x).
% USAGE: [x1,x2] = rootsearch(func,a,d,dx)
% INPUT:
% func = handle of function that returns f(x).
% a,b = limits of search.
% dx = search increment.
% OUTPUT:
% x1,x2 = bounds on the smallest root in (a,b);
% set to NaN if no root was detected
x1 = a; f1 = feval(func,x1);
x2 = a + dx; f2 = feval(func,x2);
while f1*f2 > 0.0
    if x1 >= b
        x1 = NaN; x2 = NaN; return
    end
    x1 = x2; f1 = f2;
    x2 = x1 + dx; f2 = feval(func,x2);
end
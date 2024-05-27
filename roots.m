function [root, nroots] = roots( mu )
% source: Kiusalaas, J., 2010. Numerical Methods in Engineering with MATLAB, Cambridge university press
func = @(x) (x * tan(x)-mu);
a = 0.0; b = 4000.0; dx = 0.001;
nroots = 0;
while 1
    [x1,x2] = rootsearch(func,a,b,dx);
    if isnan(x1)
        break
    else
        a = x2;
        x = bisect(func,x1,x2,1);
        if ~isnan(x)
            nroots = nroots + 1;
            root(nroots) = x;
        end
    end
end


end


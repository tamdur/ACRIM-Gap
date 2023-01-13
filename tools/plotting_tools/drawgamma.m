function y = drawgamma(m,b,a,n)
    y = zeros(n,1);
    for ii = 1:n
        z = m + randn(b,1);
        y(ii) = a./(z'*z); 
    end
end
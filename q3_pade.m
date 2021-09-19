format long e
syms x
range_L = [1,3,5,10];  % max order of numerator
f = @(x) (1+x)^(1/2);   % equation of f_1(x)
x = linspace(1,100,10000);
for i = 1:numel(range_L)
    L = range_L(i);
    % R stores the pade approximations at different x
    R = zeros(numel(x),1);
    % cache the coefficients
    [p,q] = CoefPade(f,0,L,L);
    for j = 1: numel(x)
        R(j) = polyval(p,x(j))/polyval(q,x(j));   % find the approximation at x(j)
    end
    plot(x,R,'DisplayName',sprintf('L=%g',L))
    hold on
end
plot(x,(1+x).^(1/2),'DisplayName', 'f(x)')
legend()

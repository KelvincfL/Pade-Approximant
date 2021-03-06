format long
syms x
f = @(x) (1+x+x^2)^(1/2);   % fill in equation of f_i(x), i=1,3,4,5,6
range_L = [8,9];      % range of L we want to look at
for i = 1:numel(range_L)    
    L = range_L(i);
    [p,q] = CoefPade(f,0,L,L); % coefficients of numerator and denominator
    zeros = complex(RootFinder(p));  % zeros are roots of the numerator poly
    poles = complex(RootFinder(q));  % poles are roots of the denominator poly
    figure()
    plot(zeros,'o','DisplayName','zeros')
    hold on
    plot(poles,'x','DisplayName','pole')
    legend()
end


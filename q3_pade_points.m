format long e
syms x
range_L = [1:30];  % max order of numerator
f = @(x) (1+x)^(1/2);   % equation of f_1(x)
x = [100];
for j = 1:numel(x)
        % R stores the pade approximations at different L
        R = zeros(numel(range_L),1);
        error = zeros(numel(range_L),1);
    for i = 1:numel(range_L)
        L = range_L(i);
        % calculating R,error at diferent levels of L
        R(i) = CoefSolver(f,0,L,L,x(j));
        error(i) = abs(R(i)-(1+x(j))^(1/2));
    end
        % plotting the error at x(j) against range of L
        plot(range_L,log(error),'DisplayName',sprintf('x=%g',x(j)))
        coefficients = polyfit((range_L(1:15)), log(error(1:15)), 1);
        slope = coefficients(1);
        disp(slope)
        legend()
        figure()
end

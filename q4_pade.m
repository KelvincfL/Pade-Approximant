format long e
x = linspace(0,20,200);
range_L = [15];  % order of L we wish to plot at
for l = 1:numel(range_L)
    L = range_L(l);
    c = zeros(1,2*L+1);     % stores from c_0 to c_L+M+1
    R = zeros(numel(x),1);
    % find the coefficients directly
    for i = 1:2*L+1
        c(i)= ((-1)^(i-1))*factorial(i-1);
    end
    for j = 1:numel(x)
        R(j) = CoefSolverPS(c,L,L,x(j));    % find pade estimation at x(j)
    end
    plot(x,R,'DisplayName',sprintf('L=%g',L))
    hold on
end
% the numerical integration result from project sheet
x_num = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];
y_num = [0.91563334, 0.85211088, 0.80118628, 0.75881459, 0.72265723, 0.69122594, 0.66351027...
    0.63879110, 0.61653779, 0.59634736, 0.46145532, 0.38560201, 0.33522136, 0.29866975,...
    0.27063301, 0.24828135, 0.22994778, 0.21457710, 0.20146425, 0.19011779,0.18018332,...
    0.17139800, 0.16356229, 0.15652164, 0.15015426, 0.14436271, 0.13906806, 0.13420555, 0.12972152];
plot(x_num,y_num,'DisplayName', 'f2(x)')
legend()

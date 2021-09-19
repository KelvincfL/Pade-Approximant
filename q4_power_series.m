% truncation power
range_N = [4,7];      % range of N we want to test
% plot against x from 0 to 20 with increment of 0.1
x = linspace(0,1,50); 
for i = 1:numel(range_N)
    result = zeros(numel(x),1);
    N = range_N(i);     % N is the number of terms we wish to add
    c = zeros(N,1);   % c stores the value of c_k for k=1,...N
    % find the coefficients directly
    for j = 1:N
        c(j)= ((-1)^(j))*factorial(j);
    end
    % reverse the order of c_k to descending
    c = c(N:-1:1);
    % c_0 = 1
    c(N+1) = 1;
    % evaluate the power series at points of x
    for l = 1:numel(x)
        result(l) = polyval(c,x(l));
    end
    plot(x,result,'DisplayName', sprintf('N=%g',N))
    hold on
end
% compare with f_1(x)
x_num = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
y_num = [0.91563334, 0.85211088, 0.80118628, 0.75881459, 0.72265723, 0.69122594, 0.66351027...
    0.63879110, 0.61653779, 0.59634736];
plot(x_num,y_num,'DisplayName', 'f2(x)')
legend()
ylim([0,1])
xlim([0,1])
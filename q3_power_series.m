% power series estimate for f_1(x), expansion about x=0
range_N = [4];      % range of N we want to test
% plot against x from 1 to 100 with increment of 0.1
x = linspace(1,100,10000); 
for i = 1:numel(range_N)
    result = zeros(numel(x),1);
    N = range_N(i);     % N is the number of terms we wish to add
    c = zeros(N,1);   % c stores the value of c_k for k=1,...N
    c(1) = 1/2;  % value of c_1
    for k = 2:N
        % temp is (2k-3)*(2k-1)*...*3*1
        temp = 1;
        for j = 1:2:(2*k-3)
            temp = temp*j;
        end
        % calculate c_k
        c(k) = ((-1)^(k-1))*temp/...
            (factorial(k)*2^(k));
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
plot(x,(1+x).^(1/2),'DisplayName','f(x)')
legend()

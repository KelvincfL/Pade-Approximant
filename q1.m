format long
range_N = [10:10:150]';     % the range of N we wish to test with
NumTestN = numel(range_N);
partial_sum = zeros(NumTestN,1);   % to store the results for partial sums
errors = zeros(NumTestN,1);        % to store the results for error term
for i = 1:numel(range_N)
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
    temp_sum = 1+sum(c);     % need to add c_0 to the partial sum  
    % compare values between sqrt(2) and the partial sum
    temp_error = abs(sqrt(2)-temp_sum);  
    partial_sum(i) = temp_sum;
    errors(i) = temp_error;
end
disp([range_N,partial_sum,errors])


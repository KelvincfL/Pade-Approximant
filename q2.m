format long e
syms x
range_L = 1:20;  % max order of numerator
f = @(x) (1+x)^(1/2);   % equation of f_1(x)
% R stores the pade approximations at different values of L
R = zeros(numel(range_L),1);
for i = 1:numel(range_L)
    L = range_L(i);
    R(i) = CoefSolver(f,0,L,L,1);       % find the approximation at x=1
end
disp([range_L',R,abs(sqrt(2)-R)])
% plot log-log of L against error
plot(log(range_L), log(abs(sqrt(2)-R)))
xlabel('log(L)')
ylabel('log(error)')
% plot log(error) against L
figure()
plot(range_L, log(abs(sqrt(2)-R)))
xlabel('L')
ylabel('log(error)')

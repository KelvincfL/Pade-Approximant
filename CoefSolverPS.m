function [R] = CoefSolverPS(c,L,M,y)
% parameters:
% R is the pade approximant at y
% c is the power series truncated to L+M+1 terms
% L is the max order term of the numerator
% M is the max order term of the denominator
% y is the value we want to evaluate the approximant at
% symbolic calculation to avoid numerical error
c = sym(c);
% A is the matrix of simultaneous equation in c
A = zeros(M);
% B is vector of (-c_L+1,...,-C_L+M) 
B = [-c(L+2:L+M+1)]';    
% c(L+2)=c_(L+1) and c(L+M+1)=c_(L+M)
% since indexing start from 1 and c(1)=c_0

for r = L+1:L+M
    upper_limit = min(r,M);
    % update the matrix with coef c_(r-k)
    % if index larger than min(r,M), then A(index)=0
    for k = 1 : upper_limit
        A(r-L,k) = c(r-k+1);
    end
end
% symbolic calculation to avoid numerical error
A = sym(A);
% q is the coefficients of the denominator start from q1
% q is the solution to Ax=B, row vector
q = mldivide(A,B)';

% p is the coefficients of the numerator start from p0
% now we find the coefficients p
p = zeros(1,L+1);
p(1) = c(1);
% equation 5
for k = 1:L
    upper_limit = min(k,M);
    % temporary result of sum of q_s*c_(k-s)
    temp = 0;
    for s = 1 : upper_limit
        temp = temp + q(s)*c(k-s+1);
    end
    p(k+1) = c(k+1)+ temp;
end

% reverese order to descending powers
p = p(L+1:-1:1);
q = q(M:-1:1);
% set the q_0 to be 1 for easier conversion to polynomials
q(M+1) = 1;
q = double(q);
% evaluate R at point y
R = polyval(p,y)/polyval(q,y);
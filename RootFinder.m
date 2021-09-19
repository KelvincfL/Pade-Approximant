function [r] = RootFinder(f)
% f is the coefficients of the polynomial we want to solve
% f should be in descending order
% r is the complex roots of the polynomial

r = roots(f);

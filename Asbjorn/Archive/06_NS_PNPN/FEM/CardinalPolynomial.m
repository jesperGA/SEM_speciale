function [ l ] = CardinalPolynomial(knots,i,t)
%kardinalpoly: Calculates the function values of i'th cardinal 
%              polynomial in all t-values based on the knots 
%
% Input:
%    i      the index (between 1 and n+1) of the cardinal polynomial
%    knots: the knots [x_1 ... x_(n+1)]
%    t      a row [t1 ... tm] of all the values that the 
%           i'th cardinal polynomial is to be evaluated in  
% Output:
%    l      a row with the m function values of the
%           i'th cardinal polynomial

l=ones(size(t)); % Accepting both a row and a column
m=length(knots);
for j=1:m
    if j~=i
        l=l.*(t-knots(j))./(knots(i)-knots(j));
    end
end


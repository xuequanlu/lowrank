function [ rowNum_patch, colNum_patch, Patchsize ] = decomposeInteger( K, allFactors )
%DECOMPOSEINTEGER Summary of this function goes here
%   Detailed explanation goes here

%global allFactors;

% decompose a number into two integers which are the closest
factors = factor(K);% from small to big
n = numel(factors);
% check if tow factors, and the smallest one is <=3
while n==2 && factors(1)==3 && abs(factors(1)-factors(2))>=6
    K = K-3;
    factors = factor(K);
    n = numel(factors);
end

Patchsize = K/3;


rowNum_patch = allFactors(K/3,1);
colNum_patch = allFactors(K/3,2);


end


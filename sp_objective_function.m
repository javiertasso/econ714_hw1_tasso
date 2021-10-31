function [outputArg1] = sp_objective_function(X, omega, alpha, lambda)
% This is the objective function of the social planner in problem 4
%   The function takes the consumption of each good and each agent (m*n)
%   and gives us the social welfare 
%   The input is a matrix X m*n where m is the number of goods and n is the
%   number of people
%   Also I need to specify the matrix of parameters omega and vectors alpha
%   and lambda 

[m,n] = size(X);
M = zeros(m,n); 
 


for ii = 1:m
    
   for jj = 1:n
       
       M(ii,jj) = alpha(jj) * X(ii,jj)^(1+omega(ii,jj))/(1+omega(ii,jj)); 
       
   end
   
end

outputArg1 = sum(M,1)*lambda; 

clear M  
end


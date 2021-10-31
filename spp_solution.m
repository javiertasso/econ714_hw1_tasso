function [outputArg1, outputArg2] = spp_solution(endow, var, omega, alpha, lambda)
% spp_solution: This function computes the solution of the social planner
% problem
%   I need a matrix of initial endowments
%   The variable "var" is the standard deviation we use in the algorithm 
%   Omega, alpha and lambda are the parameters associated with the
%   inidividuals 

% Random walk Metropolis-Hastings - Modified 

% Set up the algorithm 
[m,n] = size(endow); 
x_prev = ones(m,n) .* sum(endow,2) / m ;  
diff = 1; 
max_n_it = 10000000;
tol = 10^(-6);
it = 0;
x_star = x_prev; 
f_values = zeros(max_n_it,1);
f_values(1) = sp_objective_function(x_prev, omega, alpha, lambda); 

% The initial guess is allocating the total endowment of each good equally
% between all the people in this economy. This is indeed the solution when
% the problem is symmetric, but this will not be the solution once we start
% including heterogeneities 

% As long as the difference is higher than the tolerance level, continue  
 
while diff > tol && it < max_n_it 
   
   
   % We go good by good. Start with good 1, then good 2 and all the way up to good m  
   % For each good, we will pick two people at random and reallocate a
   % random amount of the good 
   
    
   for jj = 1:m 
        
        it = it + 1;
        permut = randperm(n); 
        reassign_from_to = permut(:,1:2); % Ramdomly choose the reallocation 
        epsilon = normrnd(0,var); 
        x_star(jj,reassign_from_to(1)) = max(x_prev(jj,reassign_from_to(1)) - epsilon,0);
        x_star(jj,reassign_from_to(2)) = max(x_prev(jj,reassign_from_to(2)) + epsilon,0);
        x_star(jj, :) = x_star(jj,:) ./ sum(x_star(jj,:)) .* sum(endow(jj,:),2); 
        
        p = min(1,sp_objective_function(x_star,omega, alpha, lambda) / sp_objective_function(x_prev, omega, alpha, lambda)); 
        bern_draw = binornd(1,p,1); 
        x =  bern_draw * x_star + (1-bern_draw) * x_prev; % Ramdomly draw from x_star or the previous value 
        
       if sp_objective_function(x, omega, alpha, lambda) >= sp_objective_function(x_prev, omega, alpha, lambda)

            x_prev = x;
            f_values(it+1,1) = sp_objective_function(x, omega, alpha, lambda);
            % If there was an improvement, then keep this new draw 

       else

            f_values(it+1,1) = sp_objective_function(x_prev, omega, alpha, lambda);
            % Do nothing if there was no improvement 

       end 
       
       % We will stop if we see that the objective function has not changed
       % after the following number of iterations
       % This could be modified to address the fact that as we have more
       % goods we will be performing more iterations
       
       if it > 100000 
        
            diff = abs(f_values(it-100000,1) - f_values(it,1));
        
       end
        
        
   end
   
   
end

% The function returns the solution and the number of iterations 

outputArg1 = x_prev;
outputArg2 = it; 

end


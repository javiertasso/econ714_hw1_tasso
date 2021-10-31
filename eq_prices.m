function [outputArg1] = eq_prices(endow, alpha, omega)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

tol = 10^(-6); 
tol_norm = 10^(-3);
max_it = 10000000; 
[m,n] = size(endow);
% omega = -2 * ones(m,n);
% alpha = transpose(1:1:m); 
alpha = transpose(alpha); 
p = ones(1,m); 
X = zeros(m,n); 
excess_demand_norm = 1; 
iterations = 0; 

while excess_demand_norm > tol_norm && iterations < max_it 


    for ii = 1:n 

        % Start one person 

        % Calculate her income and an initial guess
        income = p * endow(:,ii); 
        x = transpose(1./(m * p(1,:)) * income); 



        % Set up the while
        n_it = 0; 
        norm_of_F = 1; 

        while norm_of_F > tol && n_it <= max_it

            n_it = n_it+1; 

            % Set up the system we want to solve
            F = ones(m,1); 
            F(1,1) = p*x - p*endow(:,ii); 
            F(2:m,1) = (alpha(2:m,1) / alpha(1,1)) .* (x(2:m,1).^(omega(2:m,ii)) / x(1,1)^(omega(1,ii))) ...
                    - transpose(p(1,2:m)) / p(1,1); 

            % Set up the Jacobian matrix 
            J(1,:) = p; 
            J(2:m,1) = - (omega(1,ii) / x(1,1)) * (F(2:m, 1) + transpose(p(1,2:m)) / p(1,1)); 
            J(2:m,2:m) = diag((omega(2:m,ii) ./ x(2:m,1)) .* (F(2:m, 1) + transpose(p(1,2:m)) / p(1,1)));

            % Norm of F
            norm_of_F = norm(F); 
            
                        % I do not want to move too fast
            
            step = -J\F;
            
            % Here I change the step. I do not want to move more than 5% of
            % everything I could spend in that good 
            
            for jj = 1:m
                
               if step(jj,1) > income / (20 * p(1,jj))
                   
                   step(jj,1) = income / (20 * p(1,jj));
                   
               end
               
               if step(jj,1) < -income / (20 * p(1,jj))
                   
                   step(jj,1) = -income / (20 * p(1,jj));
                   
               end
                
            end

            % New point 
            x = max(x + step, 10^(-6));


        end

        X(:,ii) = x; 

    end

    % Now check if the prices are the correct ones 
    aggregate_demand = sum(X,2); 
    aggregate_endow = sum(endow, 2); 

    % Adjust the price 
    for jj = 2:m 

       if aggregate_demand(jj,1) - aggregate_endow(jj,1) > tol  

           p(1,jj) = p(1,jj) + abs(normrnd(0,(aggregate_demand(jj,1) - aggregate_endow(jj,1))/...
               (1 * aggregate_endow(jj,1))));

       end

       if aggregate_demand(jj,1) - aggregate_endow(jj,1) < -tol

          p(1,jj) = max(p(1,jj) - abs(normrnd(0,(aggregate_endow(jj,1) - aggregate_demand(jj,1))/...
               (1 * aggregate_endow(jj,1)))),10^(-6)); 

       end

    end
    
    % Excess demand norm 
    
    excess_demand_norm = norm(aggregate_demand - aggregate_endow);  
    
    iterations = iterations + 1;   
    
    

end 


outputArg1 = p;


end


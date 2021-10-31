% -------------------------------------------------------------------------
% Econ 714
% Homework 1 
% Fall 2021 
% Javier Tasso
% -------------------------------------------------------------------------
clearvars
clc 
cd 'C:\Users\Javier\Dropbox\Aplicaciones\Overleaf\Econ714-PS1'

% -------------------------------------------------------------------------
% Part 2
% -------------------------------------------------------------------------

% Set up parameter values
T = 100;
rho = 0.04;
lambda = 0.02;
n = [100, 1000, 10000, 100000];
integral_midpoint = zeros(length(n),1);
integral_midpoint_t = zeros(length(n),1);
integral_simpson = zeros(length(n),1);
integral_simpson_t = zeros(length(n),1);
integral_trapezoid = zeros(length(n),1);
integral_trapezoid_t = zeros(length(n),1);
integral_montecarlo_crude = zeros(length(n),1);
integral_montecarlo_crude_t = zeros(length(n),1);
integral_montecarlo_importance = zeros(length(n),1);
integral_montecarlo_importance_t = zeros(length(n),1);

% Define the function we want to integrate
F = @(x) exp(-rho*x)*(-exp(-(1-exp(-lambda*x))));

% Plot the function to see what is going on 
xx = 1:1:100;
f_values = zeros(100,1);

for ii = 1:100
    
    f_values(ii,1) = F(xx(ii));
    
end

figure(1)
plot(xx, f_values, 'linewidth', 2);
hold on
h = area(xx, f_values, 0);
set(h,'FaceColor',[91, 207, 244] / 255)
xlabel('x')
ylabel('y')
title('Plot of the area we want to calculate')
hold off

saveas(gcf,'econ714_homework1_question2_plot.png');

% Remove things 
close(figure(1)) 
clear f_plot f_values h ii xx

% This is the main loop that generates all the objects we need
for jj = 1:length(n)

    % Midpoint rule -----------------------------------------------------------

    tic

    % Set up the grid
    grid = linspace(0, T, n(jj));
    step_size = (T)/(n(jj)-1);

    % Generate a vector to store values of F
    f_values = zeros(n(jj)-1, 1); 

    % Store values of F 
    for ii = 1:(n(jj)-1)

        f_values(ii,1) = F((grid(1,ii)+grid(1,ii+1))/2);

    end

    % Calculate the integral
    integral_midpoint(jj) = step_size * sum(f_values); 

    % Clean things
    clear f_values grid ii step_size 
    integral_midpoint_t(jj) = toc; 

    % -------------------------------------------------------------------------

    % Trapezoid rule ----------------------------------------------------------

    tic

    % Set up parameters
    a = 0;
    b = T; 

    % Store values of f to sum 
    f_values = zeros(n(jj)-1, 1); 

    % Fill f_values
    for ii = 1:(n(jj)-1)

        f_values(ii,1) = F(a + ii * (b-a) / (n(jj)-1));

    end

    % Calculate the integral
    integral_trapezoid(jj) = (b-a) / n(jj) * (F(a)/2 + sum(f_values) + F(b)/2);

    % Clean things
    clear f_values ii a b 
    integral_trapezoid_t(jj) = toc; 

    % -------------------------------------------------------------------------

    % Simpson rule ------------------------------------------------------------

    tic

    % Set up parameters 
    a = 0;
    b = T;
    h = (b-a) / n(jj);

    % Store values of f to sum 
    f_values = zeros(n(jj)-1, 1); 

    % Fill f_values
    for ii = 1:(n(jj)-1)

        f_values(ii,1) = F(a + ii * h);

    end

    f_even = f_values(2:2:end);

    % Calculate the integral 
    integral_simpson(jj) = h/3 * (F(a) + 2 * sum(f_values) + 2 * sum(f_even) + F(b));

    % Remove things
    clear a b f_even f_values h ii 
    integral_simpson_t(jj) = toc; 

    % -------------------------------------------------------------------------

    % Monte carlo - Crude method ----------------------------------------------

    tic

    % Set up parameters
    a = 0;
    b = T;
    x_rand = (b-a) * rand(n(jj), 1);
    f_values = zeros(n(jj), 1);

    for ii = 1:n(jj)

        f_values(ii,1) = F(x_rand(ii));

    end

    % Calculate the integral
    integral_montecarlo_crude(jj) = ((b-a)/n(jj)) * sum(f_values);

    % Clean things
    clear f_values a b ii x_rand
    integral_montecarlo_crude_t(jj) = toc; 

    % -------------------------------------------------------------------------

    % Monte carlo - Importance sampling ---------------------------------------

    tic

    % Set up parameters
    a = 0;
    b = T;
    x_values = a:((b-a)/(n(jj)-1)):b;
    f_values_abs = zeros(n(jj),1); 

    for ii = 1:n(jj) 

        f_values_abs(ii,1) = abs(F(x_values(ii)));

    end

    % Generate the cdf
    total_sum = sum(f_values_abs); 
    cdf = cumsum(f_values_abs) / total_sum;

    % Clean some things
    clear total_sum f_values_abs x_values ii

    % Generate random numbers using the cdf and evaluate f/p
    random_numbers = rand(n(jj),1); 
    f_values = zeros(n(jj),1);

    for ii = 1:n(jj) 

        [~, index] = min(abs(cdf - random_numbers(ii,1)));
        element = index / n(jj) * (b-a); 

        if index == 1 

            f_values(ii,1) = F(element) / (cdf(index));

        else

            f_values(ii,1) = F(element) / (cdf(index)-cdf(index-1));

        end

    end

    % Calculate the integral 
    integral_montecarlo_importance(jj) = ((b-a)/(n(jj)^2)) * sum(f_values);

    % Remove things
    clear a b cdf index element f_values ii random_numbers 
    integral_montecarlo_importance_t(jj) = toc; 

    % -------------------------------------------------------------------------

end

% Clean
clear jj lambda rho T F 

% Create a matrix to export output 
format short 
values_l = [integral_midpoint, integral_midpoint_t, integral_trapezoid, integral_trapezoid_t, ...
    integral_simpson, integral_simpson_t, integral_montecarlo_crude, integral_montecarlo_crude_t, ...
    integral_montecarlo_importance, integral_montecarlo_importance_t]';
 
values = round(values_l, 4);

row_names = {'Midpoint'; 'Time (in seconds)'; 'Trapezoid'; 'Time (in seconds)'; 'Simpson'; ...
    'Time (in seconds)'; 'Montecarlo Crude'; 'Time (in seconds)'; 'Montecarlo Importance'; ...
    'Time (in seconds)'};
  
vv = [row_names, num2cell(values)];
vvv = cell2table(vv, 'VariableNames', {'Method', 'N100', 'N1000', 'N10000', 'N100000'});
 

% Export output 
writetable(vvv,'econ714_homework1_question2_output.csv')

% Clean 
clearvars

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
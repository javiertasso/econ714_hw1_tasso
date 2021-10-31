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
% Part 3
% -------------------------------------------------------------------------

% Define function, gradient and hessian
F = @(x,y) 100 * (y-x^2)^2 + (1-x)^2;
Fx = @(x,y) -400 * x * (y-x^2) - 2 * (1-x);
Fy = @(x,y) 200 * (y-x^2);
Fxx = @(x,y) -400 * y + 1200 * x^2 + 2;
Fyy = @(x,y) 200;
Fxy = @(x,y) -400 * x;
 
% Tolerance level
tol = 10^(-6); 
max_it = 500000;

% Set up the loop
n_initial_guesses = 4;
initial_guess = [2,2,-2,-2;3,-3,3,-2];
bfgs_n_it = zeros(n_initial_guesses,1);
bfgs_time = zeros(n_initial_guesses,1);
bfgs_solution = zeros(2,n_initial_guesses); 
cong_n_it = zeros(n_initial_guesses,1);
cong_time = zeros(n_initial_guesses,1); 
cong_solution = zeros(2,n_initial_guesses); 
nw_n_it = zeros(n_initial_guesses,1);
nw_time = zeros(n_initial_guesses,1); 
nw_solution = zeros(2,n_initial_guesses); 
steepest_d_n_it = zeros(n_initial_guesses,1);
steepest_d_time = zeros(n_initial_guesses,1);
steepest_d_solution = zeros(2,n_initial_guesses); 


for ii = 1:n_initial_guesses

  % Initial guess 

  initial_guess_1 = initial_guess(:,ii); 
  initial_guess_2 = 0.8 * initial_guess_1;
  f_value_initial_guess = F(initial_guess_1(1), initial_guess_1(2));

  % Newton Raphson Method ---------------------------------------------------

  tic 

  % Define variables 
  n_it_nw = 0; 
  x_value_nw = initial_guess_1;
  diff = 1;
  f_value = f_value_initial_guess;

  while diff > tol && n_it_nw < max_it
     n_it_nw = n_it_nw + 1;
     x_value_nw = x_value_nw - [Fxx(x_value_nw(1), x_value_nw(2)), Fxy(x_value_nw(1), x_value_nw(2)); ...
         Fxy(x_value_nw(1), x_value_nw(2)), Fyy(x_value_nw(1), x_value_nw(2))] ... 
         \ [Fx(x_value_nw(1), x_value_nw(2));Fy(x_value_nw(1), x_value_nw(2))];
     diff = abs(f_value - F(x_value_nw(1), x_value_nw(2)));
     f_value = F(x_value_nw(1), x_value_nw(2));
        
  end

  nw_time(ii) = toc;
  nw_n_it(ii) = n_it_nw;
  nw_solution(:,ii) = x_value_nw;


  % Clean 
  clear diff f_value n_it_nw
   

  % -------------------------------------------------------------------------

  % BFGS --------------------------------------------------------------------

  tic 

  % Define variables 
  n_it_bfgs = 0; 

  % Need two guesses here so I can start approximating the inv(Hessian)
  x_value_bfgs_1 = initial_guess_1;
  x_value_bfgs_2 = initial_guess_2; 
  diff = 1;
  f_value_2 = F(initial_guess_2(1), initial_guess_2(2));

  % Use a something that resembles the true hessian 
  H_inv = inv([Fxx(1,1),0;0,200]);
      

  while diff > tol && n_it_bfgs < max_it
      
     n_it_bfgs = n_it_bfgs + 1;
     s = x_value_bfgs_2 - x_value_bfgs_1; 
     y = [Fx(x_value_bfgs_2(1),x_value_bfgs_2(2));Fy(x_value_bfgs_2(1),x_value_bfgs_2(2))] - ...
         [Fx(x_value_bfgs_1(1),x_value_bfgs_1(2));Fy(x_value_bfgs_1(1),x_value_bfgs_1(2))];
     H_inv = (eye(2) - s*transpose(y) / (transpose(y)*s)) * H_inv * (eye(2) - y*transpose(s) ... 
         / (transpose(y)*s)) + s*transpose(s) / (transpose(y)*s);
     x_value_bfgs_1 = x_value_bfgs_2;
     x_value_bfgs_2 = x_value_bfgs_2 - H_inv * [Fx(x_value_bfgs_2(1),x_value_bfgs_2(2));...
         Fy(x_value_bfgs_2(1),x_value_bfgs_2(2))];
     f_value_1 = f_value_2;
     f_value_2 = F(x_value_bfgs_2(1), x_value_bfgs_2(2));
     diff = abs(f_value_2 - f_value_1);
      
  end

  bfgs_time(ii) = toc; 
  bfgs_n_it(ii) = n_it_bfgs;
  bfgs_solution(:,ii) = x_value_bfgs_2;

  % Clean
  clear y s H_inv f_value_1 f_value_2 diff x_value_bfgs_1 

  % -------------------------------------------------------------------------

  % Steepest descent method -------------------------------------------------

  tic 

  % Define variables
  x_value_steepest_descent = initial_guess_1;
  n_it_steepest_descent = 0;
  diff = 1;
  X_grid = zeros(2, max_it);
  X_grid(:,n_it_steepest_descent+1) = x_value_steepest_descent;


  while diff > tol && n_it_steepest_descent < max_it
      
     n_it_steepest_descent = n_it_steepest_descent + 1; 
     gradient = [Fx(x_value_steepest_descent(1),x_value_steepest_descent(2)); ...
         Fy(x_value_steepest_descent(1),x_value_steepest_descent(2))];
     H = [Fxx(x_value_steepest_descent(1),x_value_steepest_descent(2)), ... 
         Fxy(x_value_steepest_descent(1),x_value_steepest_descent(2)); ... 
         Fxy(x_value_steepest_descent(1),x_value_steepest_descent(2)), ...
         Fyy(x_value_steepest_descent(1),x_value_steepest_descent(2))];
     step_size = transpose(gradient) * gradient / (transpose(gradient) * H * gradient);
     x_value_steepest_descent = x_value_steepest_descent - step_size * gradient;
     diff = norm(gradient);
     X_grid(:,n_it_steepest_descent+1) = x_value_steepest_descent;

  end

  steepest_d_time(ii) = toc; 
  steepest_d_n_it(ii) = n_it_steepest_descent;
  steepest_d_solution(:,ii) = x_value_steepest_descent;

  % Plot the steps
  X_grid = X_grid(:,1:n_it_steepest_descent);


  figure(1)
  plot(X_grid(1,:),X_grid(2,:),'linewidth', 2)
  hold on
  plot(1,1,'r*')
  plot(initial_guess_1(1),initial_guess_1(2),'-o')
  xlabel('x')
  ylabel('y')
  title(['Number of iterations: ',num2str(n_it_steepest_descent), ...
      ' - Time (seconds): ', num2str(steepest_d_time(ii))])
  hold off
  saveas(gcf,sprintf('econ714_homework1_question3_plot_steepest_descent_%d.png',ii));
  close(figure(1))

  % Clean 
  clear diff f_value gradient H step_size X_grid   
   

  % -------------------------------------------------------------------------

  % Conjugate Gradient Method -----------------------------------------------

  tic

  max_eigen = max(eig([Fxx(initial_guess_1(1),initial_guess_1(2)), Fxy(initial_guess_1(1),initial_guess_1(2));...
      Fxy(initial_guess_1(1),initial_guess_1(2)), Fyy(initial_guess_1(1),initial_guess_1(2))]));

  step_size = 1/max_eigen; 
  n_it_cong = 0;
  diff = 1;
  X_values_cong = zeros(2, max_it);
  X_values_cong(:,1) = initial_guess_1;

  
  while diff > tol && n_it_cong < max_it 
      
     n_it_cong = n_it_cong + 1;
     r = - [Fx(X_values_cong(1,n_it_cong), X_values_cong(2,n_it_cong)); ...
         Fy(X_values_cong(1,n_it_cong), X_values_cong(2,n_it_cong))];
     X_values_cong(:, n_it_cong+1) = X_values_cong(:, n_it_cong) + step_size * r;
     gradient = [Fx(X_values_cong(1,n_it_cong+1), X_values_cong(2,n_it_cong+1)); ...
         Fy(X_values_cong(1,n_it_cong+1), X_values_cong(2,n_it_cong+1))];
     r_new = -gradient + (norm(gradient))^2 / (norm(r))^2 * r;
     diff = norm(r_new); 
     
      
  end

  cong_time(ii) = toc;
  cong_n_it(ii) = n_it_cong; 
  cong_solution(:,ii) = X_values_cong(:,n_it_cong);


  X_grid = X_values_cong(:,1:n_it_cong); 
  figure(2)
  plot(X_grid(1,:),X_grid(2,:),'linewidth', 2)
  hold on
  plot(1,1,'r*')
  plot(initial_guess_1(1),initial_guess_1(2),'-o')
  xlabel('x')
  ylabel('y')
  title(['Number of iterations: ',num2str(n_it_cong), ' - Time (seconds): ', num2str(cong_time(ii))])
  hold off
  saveas(gcf,sprintf('econ714_homework1_question3_plot_cong_%d.png',ii));
  close(figure(2))

  % Clean things
  clear diff gradient r r_new step_size X_grid X_values_cong 

end

% Clean things 
clear F Fx Fy Fxx Fxy Fyy ii initial_guess_1 initial_guess_2 max_eigen max_it ...
    n_initial_guesses n_it_bfgs n_it_cong n_it_steepest_descent tol x_value_bfgs_2 ...
    x_value_cong x_value_nw x_value_steepest_descent f_value_initial_guess f_values_cong

v = round([transpose(nw_n_it); transpose(nw_time); transpose(bfgs_n_it); transpose(bfgs_time)],4);
row_names = {'Newton Raphson'; 'Time (in seconds)'; 'BFGS'; 'Time (in seconds)'};
vv = [row_names, num2cell(v)];
vvv = cell2table(vv, 'VariableNames', {'Method'; 'GuessI'; 'GuessII'; 'GuessIII'; 'GuessIV'});

% Export output 
writetable(vvv,'econ714_homework1_question3_output.csv')

% Clean 
clearvars 






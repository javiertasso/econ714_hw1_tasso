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

% Initial guess 
initial_guess_1 = [2;3]; 
initial_guess_2 = [1.8;2.8];
f_value_initial_guess = F(initial_guess_1(1), initial_guess_1(2));

% Newton Raphson Method ---------------------------------------------------

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

display(x_value_nw)
display(n_it_nw)

% -------------------------------------------------------------------------

% BFGS --------------------------------------------------------------------

% Define variables 
n_it_bfgs = 0; 

% Need two guesses here so I can start approximating the inv(Hessian)
x_value_bfgs_1 = initial_guess_1;
x_value_bfgs_2 = initial_guess_2; 
diff = 1;
f_value_1 = f_value_initial_guess;
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

clear y s 
display(x_value_bfgs_2)
display(n_it_bfgs)

% -------------------------------------------------------------------------

% Steepest descent method -------------------------------------------------

% Define variables
x_value_steepest_descent = initial_guess_1;
n_it_steepest_descent = 0;
diff = 1;
f_value = f_value_initial_guess; 
X_grid = zeros(2, max_it);
X_grid(:,n_it_steepest_descent+1) = x_value_steepest_descent;
% H = [Fxx(0,0), Fxy(0,0); Fxy(0,0), Fyy(0,0)];


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
   f_value = F(x_value_steepest_descent(1),x_value_steepest_descent(2));
   X_grid(:,n_it_steepest_descent+1) = x_value_steepest_descent;

end

display(n_it_steepest_descent)
display(x_value_steepest_descent)

% Maybe plot the steps
X_grid = X_grid(:,1:n_it_steepest_descent);


figure(1)
plot(X_grid(1,:),X_grid(2,:),'linewidth', 2)
hold on
plot(1,1,'r*')
plot(initial_guess_1(1),initial_guess_1(2),'-o')
xlabel('x')
ylabel('y')
title('Steepest Descent Method - Steps')
hold off
saveas(gcf,'econ714_homework1_question3_plot_steepest_descent.png');
close(figure(1))

% -------------------------------------------------------------------------

% Conjugate Gradient Method -----------------------------------------------

% I have to think about this 

d0 = - [Fx(initial_guess_1(1),initial_guess_1(2));...
    Fy(initial_guess_1(1),initial_guess_1(2))];
 

% G = @(z) 100 * (d0(2) + z - (d0(1)+z)^2)^2 + (1-d0(1)-z)^2;
 
%alpha = fminbnd(G, 0, 1000);
%x = initial_guess_1 + alpha * d0;
%r1 = - [Fx(x(1),x(2)); Fy(x(1),x(2))];
%beta1 = transpose(r1)*r1 / (transpose(d0)*d0);
%d1 = r1 + beta1 * d0 

initial_guess_1 = [1.2;1.2];

H = [Fxx(initial_guess_1(1),initial_guess_1(2)), Fxy(initial_guess_1(1),initial_guess_1(2));...
    Fxy(initial_guess_1(1),initial_guess_1(2)), Fyy(initial_guess_1(1),initial_guess_1(2))];
alpha0 = transpose(d0) * d0 / (transpose(d0) * H * d0);
x1 = initial_guess_1 + alpha0*d0;
r1 = d0-alpha0*H*d0;
beta1 = transpose(r1)*r1 / (transpose(d0)*d0);
d1 = r1 + beta1 * d0;
H = [Fxx(x1(1),x1(2)), Fxy(x1(1),x1(2));...
    Fxy(x1(1),x1(2)), Fyy(x1(1),x1(2))];
alpha1 = transpose(r1)*r1 / (transpose(d1)* H * d1);
x2 = x1 + alpha1*d1



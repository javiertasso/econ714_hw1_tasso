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
% Part 4
% -------------------------------------------------------------------------

% Simple case - all identical 
    % rows are goods
    % columns are individuals       
    
% Set up the parameters 

n = 3; % Agents
m = 3; % Goods  
omega = -2 * ones(m,n); % They are all identical 
lambda = transpose(1:1:m); % But they are valued differently by the social planner
alpha = transpose(ones(m,1)); 
endow = eye(m,n) * 100; % Matrix of endowments 
var = max(sum(endow,2))/10000; 

% Initial try 
[cons_1, iter_1] = spp_solution(endow, var, omega, alpha, lambda); 
disp(cons_1)
disp(iter_1)  
omega_1 = omega; 
p_1 = eq_prices(endow, alpha, omega); 
disp(p_1)
 
% Try another set of parameters
omega = -2 * ones(m,n);
alpha = [3,2,1]; 
[cons_2, iter_2] = spp_solution(endow, var, omega, alpha, lambda); 
disp(cons_2)
disp(iter_2) 
omega_2 = omega;  
p_2 = eq_prices(endow, alpha, omega); 
 

% Now we include some heterogeneity in the omegas 
omega = -[2 * ones(1,n); 2.2 * ones(1,n); 1.8 * ones(1,n)];
alpha = transpose(ones(m,1)); 
[cons_3, iter_3] = spp_solution(endow, var, omega, alpha, lambda); 
disp(cons_3)
disp(iter_3) 
omega_3 = omega; 
p_3 = eq_prices(endow, alpha, omega); 
disp(p_3)
 

% More heterogeneity in the omegas 
omega = - (sqrt(chi2rnd(1,m,n))/4 + 1);
[cons_4, iter_4] = spp_solution(endow, var, omega, alpha, lambda);
disp(cons_4)
disp(iter_4)
omega_4 = omega; 
p_4 = eq_prices(endow, alpha, omega); 
 
 

% Let's see what happens with more individuals
n = 10;
m = 10;
omega = -2 * ones(m,n);
lambda = transpose(1:1:m);  
alpha = ones(1,n); 
endow = eye(m,n) * 100;
[cons_5, iter_5] = spp_solution(endow, var, omega, alpha, lambda); 
disp(cons_5)
disp(iter_5) 
omega_5 = omega;  
p_5 = eq_prices(endow, alpha, omega); 
 

% Let's try something crazy 
omega = - (sqrt(chi2rnd(1,m,n))/4 + 1);
[cons_6, iter_6] = spp_solution(endow, var, omega, alpha, lambda);
disp(cons_6)
disp(iter_6)
omega_6 = omega; 
p_6 = eq_prices(endow, alpha, omega); 
 

% Export things -----------------------------------------------------------

% For the case of 3x3
% Export some data
v = round(cons_2,2);
v = [v, round(transpose(p_2),4)]; 
row_names = {'GI'; 'GII'; 'GIII'};
vv = [row_names, num2cell(v)];
vvv = cell2table(vv, 'VariableNames', {'Goods';'AgentI'; 'AgentII'; 'AgentIII'; 'Price'});
writetable(vvv,'econ714_homework1_question4_output_cons_2.csv')

% Export some data
v = round(cons_4,2);
v = [v, round(transpose(p_4),4)];
row_names = {'GI'; 'GII'; 'GIII'};
vv = [row_names, num2cell(v)];
vvv = cell2table(vv, 'VariableNames', {'Goods';'AgentI'; 'AgentII'; 'AgentIII'; 'Price'});
writetable(vvv,'econ714_homework1_question4_output_cons_4.csv')

% Export some data
v = round(omega_4,4);
row_names = {'GI'; 'GII'; 'GIII'};
vv = [row_names, num2cell(v)];
vvv = cell2table(vv, 'VariableNames', {'Goods';'AgentI'; 'AgentII'; 'AgentIII'});
writetable(vvv,'econ714_homework1_question4_output_omega_4.csv')

% For the case of 10x10
% Export some data
v = round(cons_5,2);
v = [v, round(transpose(p_5),4)];
row_names = {'GI'; 'GII'; 'GIII'; 'GIV'; 'GV'; 'GVI'; 'GVII' ...
    ; 'GVIII'; 'GIX'; 'GX'};
vv = [row_names, num2cell(v)];
vvv = cell2table(vv, 'VariableNames', {'Goods'; 'AgntI'; 'AgntII'; 'AgntIII'; 'AgntIV'; 'AgntV'; 'AgntVI'; ...
    'AgntVII'; 'AgntVIII'; 'AgntXI';'AgntX'; 'Price'});
writetable(vvv,'econ714_homework1_question4_output_cons_5.csv')

% Export some data
v = round(cons_6,2);
v = [v, round(transpose(p_6),4)];
vv = [row_names, num2cell(v)];
vvv = cell2table(vv, 'VariableNames', {'Goods';'AgntI'; 'AgntII'; 'AgntIII';'AgntIV'; 'AgntV'; 'AgntVI'; ...
    'AgntVII'; 'AgntVIII'; 'AgntXI';'AgntX'; 'Price'});
writetable(vvv,'econ714_homework1_question4_output_cons_6.csv')

% Export some data
v = round(omega_6,4);
vv = [row_names, num2cell(v)];
vvv = cell2table(vv, 'VariableNames', {'Goods';'AgntI'; 'AgntII'; 'AgntIII';'AgntIV'; 'AgntV'; 'AgntVI'; ...
    'AgntVII'; 'AgntVIII'; 'AgntXI';'AgntX'});
writetable(vvv,'econ714_homework1_question4_output_omega_6.csv')

% Clean things
clear cons_1 cons_2 cons_3 cons_4 cons_5 cons_6 endow omega m n var lambda iter_1 iter_2 ... 
    iter_3 iter_4 iter_5 iter_6 alpha row_names v vv vvv p_1 p_2 p_3 p_4 p_5 p_6

%%

% -------------------------------------------------------------------------
% Part 5
% -------------------------------------------------------------------------

% Scenario 1
m = 3;
n = 3;
endow = eye(m,n) * 100;
omega = ones(m,n) * (-2);
alpha = [3,2,1]; 
p_1 = eq_prices(endow, alpha, omega);

% Scenario 2 
endow = ones(m,n) * 100 / n;
p_2 = eq_prices(endow, alpha, omega); 
    % Redistributing doesn't change anything here
    
% Scenario 3 
endow = eye(m,n) .* [90;100;120];
p_3 = eq_prices(endow, alpha, omega);

% Scenario 4 
endow = ones(m,n) .* [90;100;120] / n;
p_4 = eq_prices(endow, alpha, omega);
    
    
% Scenario 5
endow = eye(m,n) * 100;
omega = - ones(m,n) .* [2,2.05,2.1; 2,2.05,2.1; 2,2.05,2.1];
p_5 = eq_prices(endow, alpha, omega);
    % Changing omega this way does not change prices
    
    

    
   
% Scenario 6
endow = eye(m,n) * 100;
omega = - ones(m,n) .* transpose([2,2.05,2.1; 2,2.05,2.1; 2,2.05,2.1]);
p_6 = eq_prices(endow, alpha, omega);
    % Here it changes 
    
% Scenario 7
endow = eye(m,n) * 100;
omega = ones(m,n) * (-2);
alpha = [1,2,3]; 
p_7 = eq_prices(endow, alpha, omega);
  
% Export output
v = round([p_1; p_2; p_3; p_4; p_5; p_6; p_7], 4);
row_names = {'Scenario I'; 'Scenario II'; 'Scenario III'; 'Scenario VI'; 'Scenario V';...
    'Scenario VI';'Scenario VII'};
vv = [row_names, num2cell(v)];
vvv = cell2table(vv, 'VariableNames', {'Scenarios';'GoodI'; 'GoodII'; 'GoodIII'});
writetable(vvv,'econ714_homework1_question5_output.csv')


 
% Clean
clear v vv vvv 
clearvars 
disp("The end") 


   





 





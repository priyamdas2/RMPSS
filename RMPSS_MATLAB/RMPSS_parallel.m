clear all

%%%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 3; % Dimension of simplex space

% fun denotes the function to be MINIMIZED (Ackley's 2-D function, 
% TRUE global minimum at [0,0,1]).
fun = @(theta)(-20*exp(-0.2*sqrt(0.5*((theta(1))^2+(theta(2))^2)))...
    -exp(0.5*((cos(2*pi*theta(1)))+cos(2*pi*theta(2))))+exp(1)+20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

no_loops = 1000;
maximum_iteration = 50000;
theta_cut_off = 10^(-3);
epsilon_cut_off = 10^(-3);
rho_1 = 2;
rho_2 = 1.05;
tol_fun = 10^(-15);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function_count = 0;
size_here = M;


array_of_values = zeros(maximum_iteration,1);

rng('default')
starting_point_initial = zeros(M,1);
for kk = 1:size_here
    starting_point_initial(kk) = rand(1);
end
starting_point = starting_point_initial/sum(starting_point_initial);

tic;
theta_array = zeros(no_loops, M);

for iii = 1:no_loops
    epsilon = 1;
    epsilon_decreasing_factor = rho_2; % default is 2
    if(iii == 1)
        epsilon_decreasing_factor = rho_1; % default is 2
        theta = starting_point;
    else
        theta = transpose(theta_array((iii-1),:));
    end
    M = max(size(theta,1),size(theta,2));
    
    
    for i = 1:maximum_iteration
        current_lh = fun(theta);
        if(min(ge(theta,0)) == 0)
            stop('error')
        end
        
        temp_possibility_pos = zeros(M,1);
        temp_possibility_neg = zeros(M,1);

        total_lh = zeros(2*M,1);
        matrix_update_at_h = zeros(M,2*M);
        
        
        parfor parallel_no = 1:(2*M)
            if(mod(parallel_no,2) == 1)
                positive_change_loc = (parallel_no + 1)/2;
                possibility_pos = theta;
                temp_possibility_pos = theta;
                temp_possibility_pos(positive_change_loc) = 0; % To find all significant positions except the h-th
                significant_positions = find(gt(temp_possibility_pos, theta_cut_off*ones(M,1)));
                if(isempty(significant_positions) == 1)
                    possibility_pos = theta;
                else
                    possibility_pos(positive_change_loc) = theta(positive_change_loc) + epsilon;
                    possibility_pos(significant_positions) = possibility_pos(significant_positions) - epsilon/size(significant_positions,1);
                    epsilon_temp_pos = epsilon;
                    
                    if(min(ge(possibility_pos,0)) == 0 && epsilon_temp_pos > epsilon_cut_off)
                        epsilon_temp_pos = epsilon_temp_pos/epsilon_decreasing_factor;
                        possibility_pos = theta;
                        possibility_pos(positive_change_loc) = theta(positive_change_loc) + epsilon_temp_pos;
                        possibility_pos(significant_positions) = possibility_pos(significant_positions) - epsilon_temp_pos/size(significant_positions,1);
                    end
                end
                
                if(min(ge(possibility_pos,0)) == 0 || isequal(possibility_pos,theta) == 1)
                    possibility_pos = theta;
                    total_lh(parallel_no) = current_lh;
                else
                    function_count = function_count+1;
                    total_lh(parallel_no) = fun(possibility_pos);
                end
                matrix_update_at_h(:,parallel_no) = possibility_pos;
            else
                negative_change_loc = parallel_no/2;
                possibility_neg = theta;
                temp_possibility_neg = theta;
                temp_possibility_neg(negative_change_loc) = 0; % To find all significant positions except the h-th
                significant_positions = find(gt(temp_possibility_neg, theta_cut_off*ones(M,1)));
                if(isempty(significant_positions) == 1)
                    possibility_neg = theta;
                else
                    possibility_neg(negative_change_loc) = theta(negative_change_loc) - epsilon;
                    possibility_neg(significant_positions) = possibility_neg(significant_positions) + epsilon/size(significant_positions,1);
                    epsilon_temp_neg = epsilon;
                    
                    if(min(ge(possibility_neg,0)) == 0 && epsilon_temp_neg > epsilon_cut_off)
                        epsilon_temp_neg = epsilon_temp_neg/epsilon_decreasing_factor;
                        possibility_neg = theta;
                        possibility_neg(negative_change_loc) = theta(negative_change_loc) - epsilon_temp_neg;
                        possibility_neg(significant_positions) = possibility_neg(significant_positions) + epsilon_temp_neg/size(significant_positions,1);
                    end
                end
                
                if(min(ge(possibility_neg,0)) == 0 || isequal(possibility_neg,theta) == 1)
                    possibility_neg = theta;
                    total_lh(parallel_no) = current_lh;
                else
                    function_count = function_count+1;
                    total_lh(parallel_no) = fun(possibility_neg);
                end
                matrix_update_at_h(:,parallel_no) = possibility_neg;
            end
        end
        
        [Min_here,I] = min(total_lh);
        
        if(Min_here < current_lh)
            theta = matrix_update_at_h(:,I);
        end
        array_of_values(i) =  min(Min_here,current_lh);
        % [i,theta(1), theta(2), theta(3), current_lh, array_of_values(i), epsilon]
        
        % sparsity control
        sparsity_positions = lt(theta,theta_cut_off*ones(M,1));
        garbage = sum(theta(sparsity_positions));
        if(garbage > 0)
            theta(sparsity_positions) = 0;
            rich_positions = ge(theta,theta_cut_off*ones(M,1));
            theta(rich_positions) = theta(rich_positions)+garbage/nnz(rich_positions);
        end
        
        
        if(i > 1)
            if(abs(array_of_values(i) - array_of_values(i-1)) < tol_fun)
                if(epsilon > epsilon_decreasing_factor*epsilon_cut_off)
                    epsilon = epsilon/epsilon_decreasing_factor;
                else
                    break
                end
            end
        end
        
    end
    
    theta = round((10^8)*theta)/10^8; % adjustment for computation
    
    theta_array(iii,:) = transpose(theta);
    transpose(theta);
    if(iii > 1)
        if(theta_array(iii,:) == theta_array((iii-1),:))
            break
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time required

RMPSS_algo_time = toc

% starting point

starting_point

% value at starting point

fun(starting_point)

% solution

theta
% value at solution

fun(theta)


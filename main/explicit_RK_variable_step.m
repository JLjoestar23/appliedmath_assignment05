% This function computes the value of X at the next time step for any 
% arbitrary embedded RK method also computes the next step size to use, 
% and whether or not to accept/reject the projected value of X(t+h)
% (for variable time step methods)
% 
% INPUTS:
% rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
% t: the value of time at the current step
% XA: the value of X(t)
% h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
% BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
% p: how error scales with step size (error = k*hˆp)
% error_desired: the desired local truncation error at each step
% 
% OUTPUTS:
% XB: the approximate value for X(t+h)
% num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
% h_next: the time-step size at the next iteration
% redo: False if the estimated error was less than error_desired
% True if the estimated error was larger than error_desired
function [XB, num_evals, h_next, redo] = explicit_RK_variable_step(rate_func_in, t, XA, h, BT_struct, p, desired_error)
    % parse the struct
    A = BT_struct.A;
    B = BT_struct.B;
    C = BT_struct.C;
    
    % # of stages = order of the method
    stages = length(C);

    % state dimension
    n = length(XA);

    % pre-allocate vector of approximations
    % length of vector should match the order
    K = zeros(n, stages);
    
    % compute all K_i stages
    for i=1:stages
        % evaluate the sum of a_{i,j}*k_j terms as dot product
        sum_val = K*(A(i, :)');
        % evaluate the ith K approx
        K(:, i) = rate_func_in(t + C(i)*h, XA + h*(sum_val));
    end
    
    % compute both embedded next timestep estimates
    XB1 = XA + h * (K*B(1,:)');
    XB2 = XA + h * (K*B(2,:)');

    % various constants
    alpha = 6; % step size increase factor cap
    
    % approximate the local truncation error
    err = norm(XB1-XB2);
    
    % update timestep
    h_next = min(alpha, 0.9*(desired_error/err).^(1/p)) * h;
    
    if err > desired_error
        XB = XA; % don't update next step estimate
        redo = true; % trigger a recalculation
    else
        % return XB if h is sufficient
        XB = XA + h * (K*B(1,:)');
        redo = false;
    end
    
    % update num_evals
    num_evals = stages;
end
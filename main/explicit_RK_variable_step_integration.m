% Runs numerical integration arbitrary RK method using variable time steps
% 
% INPUTS:
% rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
% tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
% X0: the vector describing the initial conditions, X(t_start)
% BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
% p: how error scales with step size (error = k*hˆp)
% error_desired: the desired local truncation error at each step
%
% OUTPUTS:
% t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
% X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
% h_avg: the average step size
% num_evals: total number of calls made to rate_func_in during the integration
function [t_list, X_list, h_avg, num_evals, fail_fraction] = explicit_RK_variable_step_integration(rate_func_in, tspan, X0, BT_struct, p, desired_error)
    
    % intialize
    t = tspan(1);
    t_list(1) = t;
    t_end = tspan(2);
    XA = X0(:);
    X_list = XA;
    h = (t_end - t) / 500; % arbitrary value to iterate upon
    num_evals = 0;
    failed_steps = 0;
    attempted_steps = 0;
    
    % while time is within the timespan
    while t < t_end
        
        % in order to hit t_end exactly
        if t + h > t_end
            h = t_end - t;
        end

        % evaluate next step
        [XB, current_evals, h_next, redo] = explicit_RK_variable_step(rate_func_in, t, XA, h, BT_struct, p, desired_error);
        
        % if |XB1-XB2| > desired_error
        if redo
            % update timestep
            h = h_next;
            failed_steps = failed_steps + 1;
            attempted_steps = attempted_steps + 1;
            continue
        else
            t = t + h; % update time
            XA = XB; % update XA value
            X_list(:, end+1) = XB; % store result as row
            t_list(end+1) = t; % append new timestep
            
            % update timestep
            h = h_next;

            attempted_steps = attempted_steps + 1;
        end
        
        % Accumulate evaluations
        num_evals = num_evals + current_evals;
    end

    % compute average step size
    h_avg = mean(diff(t_list));
    
    % compute fail fraction
    fail_fraction = failed_steps / attempted_steps;
end
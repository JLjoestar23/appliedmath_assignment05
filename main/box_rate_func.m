% Rate function encoding the nonlinear equations of motion for the 
% vibrating box.
% 
% INPUTS:
% t: the current time
% V = [x;y;theta;dxdt;dydt;dthetadt]: state vector
% box_params: a struct containing the parameters that describe the system
%   box_params.m: mass of the box
%   box_params.I: moment of inertia w/respect to centroid
%   box_params.g: acceleration due to gravity
%   box_params.k_list: list of spring stiffnesses
%   box_params.l0_list: list of spring natural lengths
%   box_params.P_world: 2 x n list of static mounting points for the 
%                       spring (in the world frame)
%   box_params.P_box: 2 x n list of mounting points for the spring 
%                     (in the box frame)
%
% OUTPUTS
% dVdt = [dxdt;dydt;dthetadt;d2xdt2;d2ydt2;d2thetadt2]: the time 
% derivative of the state vector
function dVdt = box_rate_func(t, V, box_params)
    % seperate into individual variables for convenience
    x = V(1);
    y = V(2);
    theta = V(3);
    vx = V(4);
    vy = V(5);
    vtheta = V(6);
    
    % compute acceleration components
    [ax, ay, atheta] = compute_accel(x, y, theta, box_params);
    
    % structure state derivative
    dVdt = [vx; vy; vtheta; ax; ay; atheta];
end
% Computes the linear and angular acceleration of the box given its 
% current position and orientation
% 
% INPUTS:
% x: current x position of the box
% y: current y position of the box
% theta: current orientation of the box
% box_params: a struct containing the parameters that describe the system
%   box_params.m: mass of the box
%   box_params.I: moment of inertia w/respect to centroid
%   box_params.g: acceleration due to gravity
%   box_params.k_list: list of spring stiffnesses
%   box_params.l0_list: list of spring natural lengths
%   box_params.P_world: 2 x n list of static mounting
%   points for the spring (in the world frame)
%   box_params.P_box: 2 x n list of mounting points for the spring 
%                     (in the box frame)
%
% OUTPUTS
% ax: x acceleration of the box
% ay: y acceleration of the box
% atheta: angular acceleration of the box
function [ax, ay, atheta] = compute_accel(x, y, theta, box_params)
    % assign struct fields to variables for convenience
    m = box_params.m;
    I = box_params.I;
    g = box_params.g;
    k_list = box_params.k_list;
    l0_list = box_params.l0_list;
    P_world = box_params.P_world;
    P_box = box_params.P_box;
    
    F_g = [0; -(m*g)]; % defining gravitational force
    
    num_spr = size(l0_list, 2); % # of springs
    F_spr = [0; 0]; % define net spring force
    T_spr = 0; % define net spring torque
    
    % compute spring-box contact points in the world frame
    Plist_world = compute_rbt(x, y, theta, P_box);
    
    % iterate through each spring to sum up resulting forces/torques
    for i=1:num_spr
        % indexing through each parameter
        k_i = k_list(i);
        l0_i = l0_list(i);
        PA = P_world(:, i);
        PB = Plist_world(:, i);
        
        % calculate force from each spring
        F_i = compute_spring_force(k_i, l0_i, PA, PB);
        F_spr = F_spr + F_i;
        
        % calculate torque from each spring
        r = PB - [x; y]; % vector from spring-box contact point to CoM
        T_i = cross([r; 0], [F_i; 0]);
        T_spr = T_spr + T_i(3);
    end
    
    % calculate net force from gravity and springs
    F_net = F_g + F_spr;
    T_net = T_spr;
    
    % applying linear and angular momentum principles
    ax = F_net(1)/m;
    ay = F_net(2)/m;
    atheta = T_net/I;
end
% Computes the rigid body transformation that maps a set of points in the 
% box-frame to their corresponding world-frame coordinates
% 
% INPUTS:
% x: the x position of the centroid of the box
% y: the y position of the centroid of the box
% theta: the orientation of the box
% Plist_box: a 2 x n matrix of points in the box frame
% 
% OUTPUTS:
% Plist_world: a 2 x n matrix of points describing the world-frame 
% coordinates of the points in Plist_box
function Plist_world = compute_rbt(x, y, theta, Plist_box)
    % number of points on the box
    num_pts = size(Plist_box, 2);
    % pre-allocating 2 x n for world frame coordinates
    Plist_world = zeros(2, num_pts);
    % Define rotation and translation operations
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    t = [x; y];
    % iterate through each point
    for i=1:num_pts
        % convert from body frame to world frame
        Plist_world(:, i) = R * Plist_box(:, i) + t;
    end
end
% Computes the force exerted by the spring at one of its ends
% 
% INPUTS:
% PA: the position of the first end of the spring attached to the wall
% PB: the position of the second end of the spring attached to the box
% k: the stiffness of the spring
% l0: the natural length of the spring
% 
% OUTPUTS:
% F: the vector exerted by the spring at end B
function F = compute_spring_force(k, l0, PA, PB)
    % current length of the spring
    l = norm(PB - PA);

    % unit vector pointing from PA to PB
    e_s = (PB - PA) / l;

    % Force exerted by spring at point B
    F = -k*(l - l0) * e_s;
end
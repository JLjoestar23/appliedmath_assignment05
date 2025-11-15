%% Simulate Box

%define system parameters
num_springs = 4;
box_params = struct();
box_params.m = 1; % kg
box_params.I = 1/6 * box_params.m * 2; % m^4
box_params.g = 9.8; % m/s^2
box_params.k_list = 5 * ones(1, num_springs); % 8 uniform springs, N/m
box_params.l0_list = 0.2 * ones(1, num_springs); % m
box_params.P_world = 4*[-1, 1, 1, -1;
                      -1, -1, 1, 1]; % m
% spring-box contact points are at the 4 corners
box_params.P_box = [-1, 1, 1, -1; 
                    -1, -1, 1, 1]; % m

% initial conditions
x0 = 100;
y0 = 0;
theta0 = deg2rad(50);
vx0 = 0;
vy0 = 0;
vtheta0 = 0;
V0 = [x0; y0; theta0; vx0; vy0; vtheta0]; 

% integration timespan
tspan = [0 30];

simulate_box(V0, tspan, box_params);

%% Define linear model

% solve for equilibrium point
solver_params.max_iter = 1000;
V_eq = multivariate_newton_solver(@(V) box_rate_func(0, V, box_params), zeros(6, 1), solver_params);

% compute Jacobian approximation
J_approx = approximate_jacobian(@(V) box_rate_func(0, V, box_params), V_eq);

% define linearized rate function
lin_rate_func = @(t, V) J_approx * (V - V_eq);

% double check equilibrium
%simulate_box(V_eq, tspan, box_params)

%% Compare linear model with non-linear model

tspan = [0 5];

% integrating linear model
[tlist_lin, Vlist_lin] = ode45(lin_rate_func, tspan, V0);

% integrating non-linear model
nlin_rate_func = @(t, V) box_rate_func(t, V, box_params);
[tlist_nlin, Vlist_nlin] = ode45(nlin_rate_func, tspan, V0);

figure();
subplot(3, 1, 1);
plot(tlist_lin, Vlist_lin(:, 1), '.-', 'DisplayName', 'Linearized');
hold on;
plot(tlist_nlin, Vlist_nlin(:, 1), '.-', 'DisplayName', 'Non-linear');
title('X comparison');
xlabel('Time (s)');
ylabel('X (m)')
legend();
grid on;
hold off;

subplot(3, 1, 2);
plot(tlist_lin, Vlist_lin(:, 2), '.-', 'DisplayName', 'Linearized');
hold on;
plot(tlist_nlin, Vlist_nlin(:, 2), '.-', 'DisplayName', 'Non-linear');
title('Y comparison');
xlabel('Time (s)');
ylabel('Y (m)')
legend();
grid on;
hold off;

subplot(3, 1, 3);
plot(tlist_lin, Vlist_lin(:, 3), '.-', 'DisplayName', 'Linearized');
hold on;
plot(tlist_nlin, Vlist_nlin(:, 3), '.-', 'DisplayName', 'Non-linear');
title('Theta Comparison');
xlabel('Time (s)');
ylabel('Theta (rads)')
legend();
grid on;
hold off;

%% Modal analysis

% shape Q into 3x3 from the Jacobian approx
A = J_approx;
Q = A(4:6, 1:3);

% solve eigenvalue/eigenvector problem
[U_mode, w_n] = eigs(Q);

% small pertubation
epsilon = 0.1;

V0 = V_eq + epsilon*[U_mode; zeros(3,3)];

tspan = [0 2];

% integrating non-linear model
nlin_rate_func = @(t, V) box_rate_func(t, V, box_params);
[tlist_nlin_1, Vlist_nlin_1] = ode45(nlin_rate_func, tspan, V0(:,1));
[tlist_nlin_2, Vlist_nlin_2] = ode45(nlin_rate_func, tspan, V0(:,2));
[tlist_nlin_3, Vlist_nlin_3] = ode45(nlin_rate_func, tspan, V0(:,3));

%{
figure();
hold on;
plot(tlist_nlin_1, Vlist_nlin_1(:, 3), 'DisplayName', 'Mode1') % theta mode shape
plot(tlist_nlin_2, Vlist_nlin_2(:, 1), 'DisplayName', 'Mode2') % x mode shape
plot(tlist_nlin_3, Vlist_nlin_3(:, 2), 'DisplayName', 'Mode3') % y mode shape
legend();
hold off;
%}

figure();
subplot(3, 1, 1);
plot(tlist_nlin_2, Vlist_nlin_2(:, 1), '.', 'DisplayName', 'Non-linear Numerical');
hold on;
x_modal = V_eq(1) + epsilon * U_mode(1,2) * cos(sqrt(norm(w_n(2, 2))) * tlist_nlin_2);
plot(tlist_nlin_2, x_modal, '-', 'DisplayName', 'Linear Analytical')
title('X comparison');
xlabel('Time (s)');
ylabel('X (m)')
legend();
grid on;
hold off;

subplot(3, 1, 2);
plot(tlist_nlin_3, Vlist_nlin_3(:, 2), '.', 'DisplayName', 'Non-linear Numerical');
hold on;
y_modal = V_eq(2) + epsilon * U_mode(2,3) * cos(sqrt(norm(w_n(3, 3))) * tlist_nlin_3);
plot(tlist_nlin_3, y_modal, '-', 'DisplayName', 'Linear Analytical')
title('Y comparison');
xlabel('Time (s)');
ylabel('Y (m)')
legend();
grid on;
hold off;

subplot(3, 1, 3);
plot(tlist_nlin_1, Vlist_nlin_1(:, 3), '.', 'DisplayName', 'Non-linear Numerical');
hold on;
theta_modal = V_eq(3) + epsilon * U_mode(3,1) * cos(sqrt(norm(w_n(1, 1))) * tlist_nlin_1);
plot(tlist_nlin_1, theta_modal, '-', 'DisplayName', 'Linear Analytical')
title('Theta Comparison');
xlabel('Time (s)');
ylabel('Theta (rads)')
legend();
grid on;
hold off;
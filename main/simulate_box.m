function simulate_box()
    %define system parameters
    box_params = struct();
    box_params.m = 10; % kg
    box_params.I = 1/6 * box_params.m * 1; % m^4
    box_params.g = 9.8; % m/s^2
    box_params.k_list = 10 * ones(1,8); % 8 uniform springs, N/m
    box_params.l0_list = 1; % m
    box_params.P_world = [0, 1, 3, 4, 4, 3, 1, 0; 
                            1, 0, 0, 1, 3, 4, 4, 3]; % m
    % spring-box contact points are at the 4 corners
    box_params.P_box = 0.5*[-1, -1, 1, 1, 1, 1, -1 -1;
                            -1, -1, -1, -1, 1, 1, 1, 1]; % m

    % define rate function
    rate_func = @(t_in, V_in) box_rate_func(t_in, V_in, box_params);
    
    % define ICs
    x0 = 2;
    y0 = 2;
    theta0 = pi/3;
    vx0 = 0;
    vy0 = 0;
    vtheta0 = 0;
    V0 = [x0; y0; theta0; vx0; vy0; vtheta0]; 
    
    tspan = [0 60]; % 30s time span

    %run the integration
    [tlist, Vlist] = ode45(rate_func, tspan, V0);
    
    plot(tlist, Vlist(:, 1:3))
end


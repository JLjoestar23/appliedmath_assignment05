function simulate_box(x0, tspan, box_params, record_status, file_name)
    %{
    %define system parameters
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
    %}
    num_springs = box_params.num_springs;

    % define rate function
    rate_func = @(t_in, V_in) box_rate_func(t_in, V_in, box_params);
    
    % define ICs
    %{
    x0 = 2;
    y0 = 0;
    theta0 = deg2rad(10);
    vx0 = 0;
    vy0 = 6;
    vtheta0 = 0;
    V0 = [x0; y0; theta0; vx0; vy0; vtheta0]; 
    %}
    
    %tspan = [0 30]; % 30s time span

    %run the integration
    [tlist, Vlist] = ode45(rate_func, tspan, x0);

    % visualization
    figure();
    num_zigs = 5;
    w = 0.1;
    hold on;
    axis equal; axis square;
    axis([-5, 5, -5, 5]);
    xlabel('x (m)'); ylabel('y (m)');

    % initialize spring plotting structures for n springs
    for j = 1:num_springs
        spring_plot_struct(j) = initialize_spring_plot(num_zigs, w);
    end
    % initialize the box plot
    box_plot = plot(nan, nan, 'r-', 'LineWidth', 2);

    % adjustments in order to create a real time animation
    fps = 60; % desired frames-per-second
    % uniformly spaced animation times
    t_anim = (tspan(1) : 1/fps : tspan(2)); 
    V_anim = interp1(tlist, Vlist, t_anim, 'linear'); % linear interp
    dt_real = 1/fps; % seconds per frame

    % initialize video
    if record_status == true
        myVideo = VideoWriter(file_name); %open video file
        myVideo.FrameRate = 60;
        open(myVideo)
    end
    
    for i = 1:length(t_anim)
        title(sprintf('Box-Spring System Visualization (t = %.2f s)', t_anim(i)));

        % extract interpolated state
        x = V_anim(i, 1);
        y = V_anim(i, 2);
        theta = V_anim(i, 3);
    
        % compute geometry
        P1_all = box_params.P_world;
        P2_all = compute_rbt(x, y, theta, box_params.P_box);
    
        % update all springs
        for j = 1:num_springs
            update_spring_plot(spring_plot_struct(j), P1_all(:,j), P2_all(:,j));
        end

        % Update existing box line handle
        box_plot.XData = [P2_all(1, :), P2_all(1, 1)];
        box_plot.YData = [P2_all(2, :), P2_all(2, 1)];

        drawnow;

        if record_status == true
            frame = getframe(gcf); %get frame
            writeVideo(myVideo, frame);
        end

        pause(dt_real); % keeps playback consistent with physical time

    end

    if record_status == true
        close(myVideo);
    end
end

% helper function for visualizing springs
function spring_plotting_example()
    num_zigs = 5;
    w = .1;
    hold on;
    spring_plot_struct = initialize_spring_plot(num_zigs,w);
    axis equal; axis square;
    axis([-3,3,-3,3]);
    for theta=linspace(0,6*pi,1000)
        P1 = [.5;.5];
        P2 = 2*[cos(theta);sin(theta)];
        update_spring_plot(spring_plot_struct, P1, P2)
        drawnow;
    end
end

% updates spring plotting object so that spring is plotted with ends 
% located at points P1 and P2
function update_spring_plot(spring_plot_struct, P1, P2)
    dP = P2-P1;
    R = [dP(1),-dP(2)/norm(dP);dP(2),dP(1)/norm(dP)];
    plot_pts = R*spring_plot_struct.zig_zag;
    set(spring_plot_struct.line_plot,...
    'xdata',plot_pts(1,:)+P1(1),...
    'ydata',plot_pts(2,:)+P1(2));
    set(spring_plot_struct.point_plot,...
    'xdata',[P1(1),P2(1)],...
    'ydata',[P1(2),P2(2)]);
end

% create a struct containing plotting info for a single spring
% 
% INPUTS:
% num_zigs: number of zig zags in spring drawing
% w: width of the spring drawing
function spring_plot_struct = initialize_spring_plot(num_zigs,w)
    spring_plot_struct = struct();
    zig_ending = [.25,.75,1; ...
    -1,1,0];
    zig_zag = zeros(2,3+3*num_zigs);
    zig_zag(:,1) = [-.5;0];
    zig_zag(:,end) = [num_zigs+.5;0];
    for n = 0:(num_zigs-1)
        zig_zag(:,(3+3*n):2+3*(n+1)) = zig_ending + [n,n,n;0,0,0];
    end
    zig_zag(1,:)=(zig_zag(1,:)-zig_zag(1,1))/(zig_zag(1,end)-zig_zag(1,1));
    zig_zag(2,:)=zig_zag(2,:)*w;
    spring_plot_struct.zig_zag = zig_zag;
    spring_plot_struct.line_plot = plot(0,0,'k','linewidth',2);
    spring_plot_struct.point_plot = plot(0,0,'ro','markerfacecolor','r','markersize',7);
end
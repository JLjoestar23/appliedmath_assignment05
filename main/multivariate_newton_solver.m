function root_approx = multivariate_newton_solver(fun, x0, solver_params)
    % basic implementation of Newton's method for numerical root finding
    
    % unpacking solver parameters from struct
    % has default values in case they are not set
    max_iter = 1000;
    if isfield(solver_params, 'max_iter')
        max_iter = solver_params.max_iter;
    end

    dx_tol = 1e-14;
    if isfield(solver_params,'dx_tol')
        dx_tol = solver_params.dx_tol;
    end
    
    f_tol = 1e-14;
    if isfield(solver_params,'f_tol')
        f_tol = solver_params.f_tol;
    end
    
    dx_max = 1e8;
    if isfield(solver_params,'dx_max')
        dx_max = solver_params.dxmax;
    end

    approx_j = 1;
    if isfield(solver_params,'approx_j')
        approx_j = solver_params.approx_j;
    end
    
    status = 0; % convergence status
    x_n = x0(:); % initialize x_n for the first guess

    % loop until iterations reached the specified maximum number
    for i=1:max_iter

        if approx_j == 1
            % evaluate the function at the approximated root
            f_val = fun(x_n);
            % numerically compute the Jacobian at x_n
            J = approximate_jacobian(fun, x_n);
        else
            % evaluate the function at the approximated root
            [f_val, J] = fun(x_n);
        end
        
        % break if the determinant of the Jacobian is near 0, which will 
        % result in an invalid operation
        if rcond(J*J') < eps
            warning('Jacobian is near singular, method failed.');
            break
        end

        % calculate the root approximation for the next iteration
        x_next = x_n - J\f_val(:);
        
        % break if the update step is too large
        if norm(x_next - x_n) > dx_max
            warning('Updated step size is too large, method failed.');
            break
        end
        
        % check for convergence
        if norm(x_next - x_n) < dx_tol || norm(fun(x_n)) < f_tol
            status = 1; % set status to success
            break
        end

        x_n = x_next; % update x_n for next iteration

    end
    
    % if reached max number of iterations, status is failed
    if i == max_iter
        status = 0;
    end
    
    % if successful, return value of the approximated root
    if status == 1
        root_approx = x_n;
        %final_disp = strcat("Root Found, Number of Iterations: ", num2str(i));
        %disp(final_disp);
    else % warning flag if convergence failed
        warning("Convergence failed.");
        root_approx = NaN;
    end
end
%% prodQ_global
% Produce Q(z) using a global (ADMM) optimization paradigm.

%% Description
% 

function [P, q, state] = prodQ_global(z, opt_prob, state, varargin)

%% Input parameters

%% Output parameters

    %% Parse inputs

    % Short-hand variables from the values of the opt_prob parameter.
    fobj = [opt_prob.field_obj];
    pres = [opt_prob.phys_res];
    invA = {opt_prob.solve_A};
    invAd = {opt_prob.solve_A_dagger};

    N = length(fobj); % Number of fields/modes.

    % Determine the current state of the optimization.
    if isempty(state) % No previous state, use the default state.

        % Default value for x and u.
        for k = 1 : N
            mag = mean([fobj(k).alpha, fobj(k).beta], 2);
            C = fobj(k).C;
            x_default{k} = zeros(size(C,1), 1);
            u_default{k} = zeros(size(C,1), 1);
            for j = 1 : size(C, 2)
                x_default{k} = x_default{k} + ...
                                    mag(j) * C(:,j) ./ norm(C(:,j))^2;
            end
        end

        state = struct( 't', 1e3, ...
                        'rho', 10, ...
                        'newton_err_thresh', 1e-3, ...
                        'newton_max_steps', 10, ...
                        'line_search_err_thresh', 1e-6, ...
                        'vis_progress', @default_vis_progress, ...  
                        'x', {x_default}, ...
                        'u', {u_default}, ...
                        'update_u', false);
    end

    % Accept any forced parameters contained in varargin.
    for k = 2 : 2 : length(varargin)
        state = setfield(state, varargin{k-1}, varargin{k});
    end

    % Obtain state parameters.
    t = state.t;
    rho = state.rho;

    newton_max_steps = state.newton_max_steps;
    newton_err_thresh = state.newton_err_thresh; 
    line_search_err_thresh = state.line_search_err_thresh;
    vis_progress = state.vis_progress;
 
    x = state.x;
    u = state.u;
    update_u = state.update_u;

    % Check some parameters.
    if t < 1
        error('Values of t less than 1 not accepted.');
    end

    %% Update dual variables u

    % Calculate physics residual, useful for later as well, to test 
    % for improvement (and whether or not to reset u).
    for k = 1 : N
        physics_residuals{k} = pres(k).A(z) * x{k} - pres(k).b(z);
    end

    % Check if an update is needed.
    if update_u
        % Update u variables.
        for k = 1 : N 
            u{k} = u{k} + physics_residuals{k};
        end

    else 
        % Make sure to update u next time.
        % In order to force u to never be updated, update_u must be set to false
        % whenever prodQ_global is called.
        update_u = true; 
    end


    %% Compute A_dagger^-1 C values (tilde C)
    % The transformed values of C are used to efficiently calculate the Newton 
    % step.

    % Initiate A_dagger solves.
    cnt = 0;
    for k = 1 : N
        C = fobj(k).C;
        for l = 1 : size(C, 2)
            cnt = cnt + 1;
            cb{k}{l} = invAd{k}(z, C(:,l));
            done(cnt) = false;
        end
    end

    % Complete A_dagger solves.
    while ~all(done)
        cnt = 0;
        for k = 1 : N
            C = fobj(k).C;
            for l = 1 : size(C, 2)
                cnt = cnt + 1;
                [Ct_elem, done(cnt)] = cb{k}{l}(); 
                if done(cnt)
                    Ct{k}(:,l) = Ct_elem;
                end
            end
        end
    end

    % Check results.
    for k = 1 : N
        C = fobj(k).C;
        for l = 1 : size(C, 2)
            err(k, l) = norm(pres(k).A(z)' * Ct{k}(:,l) - C(:,l));
        end
    end


    %% Compose objective, gradient, and Hessian functional forms
    % This is the old way that we used to use to do it,
    % but it suffered from horrible overheads (~99%)
    % because of the use of anonymous functions.
    for k = 1 : N
        alpha = fobj(k).alpha;
        a0 = (alpha ~= 0); % Used to cancel lower barrier if alpha is 0.
        beta = fobj(k).beta;
        C = fobj(k).C;
        phi = angle(C' * x{k});

        A = pres(k).A(z);
        b = pres(k).b(z);
   
        % Scale factors used in gradient and Hessian.
        r{k} = @(x) a0 .* -exp(i*phi)./ (real(exp(-i*phi).*(C'*x)) - alpha) + ...
                    2 * (C'*x)./(beta.^2 - abs(C'*x).^2);
        s{k} = @(x) a0 .* 1./(real(exp(-i*phi).*(C'*x)) - alpha).^2 + ...
                    2./(beta.^2 - abs(C'*x).^2) + ...
                    4 * abs(C'*x).^2./(beta.^2 - abs(C'*x).^2).^2;

        % Function composition of the objective and its gradient and Hessian.
        f{k} = @(x) rho/2 * norm(A*x - b + u{k})^2 + ...
                    -1/t * sum(ln(real(exp(-i*phi).*(C'*x)) - alpha, a0) + ...
                                ln(beta.^2 - abs(C'*x).^2));
        grad{k} = @(x) rho * A' * (A*x - b + u{k}) + 1/t * C * r{k}(x); 
        multHess{k} = @(x, v) rho * A' * (A * v) + ...
                                1/t * (C * diag(s{k}(x)) * (C' * v));

%         % Alternative function composition of gradient and Hessian using C tilde.
%         grad{k} = @(x) rho * A' * (A*x - b + u{k} + 1/(rho*t) * Ct{k} * r{k}(x));
%         Hess{k} = @(x) rho * A' * ...
%                     (eye(n) + 1/(rho*t) * (Ct{k} * diag(s{k}(x)) * Ct{k}')) * A;

        % Forms needed for efficient calculation of the Newton step.
        abridged_grad{k} = @(x) (A*x - b + u{k} + 1/(rho*t) * Ct{k} * r{k}(x));
        multM{k} = @(x, v) v - ...
                    Ct{k} * (inv((diag(rho*t./s{k}(x))) + ...
                                Ct{k}'*Ct{k}) * (Ct{k}' * v));
    end

    % This is the faster way.

    % Compute the r coefficients, which are used in the gradient calculation.
    function [r] = compute_r(x, k)
        a0 = (fobj(k).alpha ~= 0); % Used to cancel lower barrier if alpha is 0.
        phi = angle(fobj(k).C' * x);
        r = a0 .* -exp(1i*phi)./ (real(exp(-1i*phi).*(fobj(k).C'*x)) - ...
                        fobj(k).alpha) + ...
                    2 * (fobj(k).C'*x)./(fobj(k).beta.^2 - abs(fobj(k).C'*x).^2);
    end

    % Compute the s coefficients, which are used in the Hessian calculation.
    function [s] = compute_s(x, k)
        a0 = (fobj(k).alpha ~= 0); 
        phi = angle(fobj(k).C' * x);
        s = a0.*1./(real(exp(-1i*phi).*(fobj(k).C'*x)) - fobj(k).alpha).^2 + ...
            2./(fobj(k).beta.^2 - abs(fobj(k).C'*x).^2) + ...
            4*abs(fobj(k).C'*x).^2./(fobj(k).beta.^2 - abs(fobj(k).C'*x).^2).^2;
    end

    % Compute the function value.
    function [f] = compute_f(x, k)
        a0 = (fobj(k).alpha ~= 0);
        phi = angle(fobj(k).C' * x);
        f = rho/2 * norm(pres(k).A(z)*x - pres(k).b(z) + u{k})^2 + ...
                    -1/t * sum(ln(real(exp(-1i*phi).*(fobj(k).C'*x)) - ...
                                    fobj(k).alpha, a0) + ...
                                ln(fobj(k).beta.^2 - abs(fobj(k).C'*x).^2));
    end

    % Precomputation to speed up the gradient calculation.
    for k = 1 : N
        grad_precompute_vector{k} = rho * pres(k).A(z)' * (-pres(k).b(z) + u{k});
        grad_precompute_matrix{k} = rho * pres(k).A(z)' * pres(k).A(z);
    end

    % Compute the gradient.
    function [grad] = compute_grad(x, k)
        a0 = (fobj(k).alpha ~= 0);
        phi = angle(fobj(k).C' * x);
        r = a0 .* -exp(1i*phi)./ ...
                    (real(exp(-1i*phi).*(fobj(k).C'*x)) - fobj(k).alpha) + ...
                    2*(fobj(k).C'*x)./(fobj(k).beta.^2 - abs(fobj(k).C'*x).^2);
        grad = grad_precompute_matrix{k} * x + grad_precompute_vector{k} + ...
                    1/t * fobj(k).C * r;
    end

    % Multiplication with the Hessian.
    function [result] = mult_Hess(x, v, k)
        result = rho * pres(k).A(z)' * (pres(k).A(z) * v) + ...
                    1/t * (fobj(k).C * diag(compute_s(x, k))*(fobj(k).C' * v));
    end

    % Computation of the abridged gradient, for efficient computation of 
    % the Newton step.
    function [ag] = compute_abridged_grad(x, k)
        ag = pres(k).A(z)*x - pres(k).b(z) + u{k} + ...
                1/(rho*t) * Ct{k} * compute_r(x, k);
    end

    % Multiplication with M, for efficiently computing the Newton step as well.
    function [result] = mult_M(x, v, k)
        result = v - Ct{k} * (inv((diag(rho*t./compute_s(x, k))) + ...
                                Ct{k}'*Ct{k}) * (Ct{k}' * v));
    end

    % Wrap everything up into cell arrays.
    for k = 1 : N
        f{k} = @(x) compute_f(x, k);
        grad{k} = @(x) compute_grad(x, k);
        multHess{k} = @(x, v) mult_Hess(x, v, k);
        abridged_grad{k} = @(x) compute_abridged_grad(x, k);
        multM{k} = @(x, v) mult_M(x, v, k);
    end


    %% Execute Newton's algorithm
    % Minimize using Newton's method.

    % Make sure we start with a feasible x.
    for k = 1 : N
        if isinf(f{k}(x{k}))
            error('Initial x is infeasible.');
        end
    end

    % Keeps track which modes have already been solved
    mode_done = false * ones(1, N);

    for j = 1 : newton_max_steps
        % Initialize the computation for the Newton step.
        unfinished_modes = find(~mode_done);
        for k = unfinished_modes
            cb{k} = invA{k}(z, -multM{k}(x{k}, abridged_grad{k}(x{k})));
        end

        % Complete the Newton step computation.
        done = mode_done; % Sets finished modes to true for computation.
        while ~all(done)
            for k = unfinished_modes
                [delta_x{k}, done(k)] = cb{k}();
            end
        end

        % Perform an iteration of the Newton algorithm.
        for k = unfinished_modes
            % Check the computation of the Newton step.
            newton_step_error = norm(multHess{k}(x{k}, delta_x{k}) + grad{k}(x{k}));

            % Calculate the Newton decrement.
            lambda = sqrt(abs(real(grad{k}(x{k})' * -delta_x{k})));

            % Perform a line search in the Newton step direction.
            step_size = line_search_convex(f{k}, grad{k}, delta_x{k}, x{k}, ...
                                            line_search_err_thresh);

            % Update x.
            x{k} = x{k} + step_size * delta_x{k};
            
            % Record the action and metrics of the step.
            progress{k}(j) = struct('fval', f{k}(x{k}), ...
                                    'grad_norm', norm(grad{k}(x{k})), ...
                                    'newton_dec', lambda^2/2, ...
                                    'step_size', step_size, ...
                                    'newton_step_err', newton_step_error);

            % Check termination condition.
            if lambda^2/2 <= newton_err_thresh
                mode_done(k) = true;
            end

        end % End of mode iteration loop.

        vis_progress(progress); % Visualize our progress.

        % Check overall termination condition.
        if all(mode_done)
            break
        end

    end % End of Newton iteration loop.

    if ~all(mode_done)
        newton_success = false;
        warning('Newton algorithm did not reach error threshold.');
    else
        newton_success = true;
    end

%     %% Check whether or not the scaled dual variables (u) need to be reset
%     % u_i is reset for every mode which has failed to decrease its physics 
%     % residual.
%     reset_u_flag = false * ones(N, 1);
%     for k = 1 : N
%         new_physics_residual{k} = pres(k).A(z) * x{k} - pres(k).b(z);
%         if norm(new_physics_residual{k}) > norm(physics_residuals{k})
%             u{k} = zeros(size(u{k})); % Reset this mode's u to zero.
%             reset_u_flag(k) = true;
%         end
%     end
% 
%     if any(reset_u_flag)
%         fprintf('[reset_u:');
%         for k = find(reset_u_flag)
%             fprintf(' %d', k);
%         end
%         fprintf(']');
%     end


    %% Form Q(z) = 1/2 || Pz - q ||^2
    P = [];
    q = [];
    for k = 1 : N
        P = [P; pres(k).B(x{k})];
        q = [q; (pres(k).d(x{k}) - u{k})];
    end
    
    % Scale by a factor of rho.
    P = sqrt(rho) * P;
    q = sqrt(rho) * q;


    %% Update state variable
    state.newton_success = newton_success;
    state.newton_progress = progress;

    state.x = x;
    state.u = u;
    state.update_u = update_u;

    % For computational efficiency, don't save this function handle.
    state = rmfield(state, 'vis_progress'); 

end % End of prodQ_global function.

%% Private functions
function [z] = ln(x, varargin)
% Custom log function which returns nan for non-zero imaginary part and 
% real part <= 0.

    if ~all(isreal(x))
        error('Cannot accept complex arguments!');
    end

    z = log(x);

    % Log gives complex results for negative x, eliminate these.
    z(find(real(x) <= 0)) = -Inf;

    % If needed, set some values to 0.
    if ~isempty(varargin)
        set_to_zero = varargin{1};
        z(find(set_to_zero == 0)) = 0;
    end

end % End of ln function.

function default_vis_progress(progress)
% Print progress.
    % Find the maximum number of iterations present over all the modes.
    for k = 1 : length(progress)
        len(k) = length(progress{k});
    end

    fprintf('%3d:', max(len)); % Print out iteration number.

    % Print out the Newton decrement for each mode.
    for k = 1 : length(progress) 
        if length(progress{k}) < max(len)
            fprintf(' [%1.1e]', progress{k}(end).newton_dec);
        else
            fprintf(' %1.3e', progress{k}(end).newton_dec);
        end
    end
    fprintf('\n');
end % End of default_vis_progress function.

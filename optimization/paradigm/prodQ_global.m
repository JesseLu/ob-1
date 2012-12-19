%% prodQ_global
% Produce Q(z) using a global (ADMM) optimization paradigm.

%% Description
%

%% How to test for Hessian?
% Hessian gives a quadratic estimation of the problem,
% one way is to simply look for quadratic convergence for the solve.
% Also, the Hessian can simply be calculated brute force.
function [P, q, status] = prodQ_global(z, opt_prob, state, varargin)

%% Input parameters

%% Output parameters

    %% Parse inputs

    % Short-hand variables from the values of the opt_prob parameter.
    fobj = [opt_prob.field_obj];
    pres = [opt_prob.phys_res];
    invA = {opt_prob.solve_A};
    invAd = {opt_prob.solve_A_dagger};

    N = length(fobj); % Number of fields/modes.
    n = length(z);

    % Determine the current state of the optimization.
    if isempty(state) % No previous state, use the default state.
        % Default value for t, used in the relaxed field objective.
        state.t = 1e-3; 
        state.rho = 10;

        state.newton_err_thresh = 1e-6;
        state.newton_max_steps = 100;
        state.line_search_err_thresh = 1e-9;
        state.vis_progress = @default_vis_progress;

        state.x = nan;
        for k = 1 : N
            state.u{k} = zeros(n, 1);
        end
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
    for k = 1 : N
        alpha = fobj(k).alpha;
        beta = fobj(k).beta;
        C = fobj(k).C;
        phi = angle(C' * x{k});

        A = pres(k).A(z);
        b = pres(k).b(z);
   
        % Scale factors used in gradient and Hessian.
        r{k} = @(x) -exp(i*phi)./ (real(exp(-i*phi).*(C'*x)) - alpha) + ...
                    2 * (C'*x)./(beta.^2 - abs(C'*x).^2);
        s{k} = @(x) 1./(real(exp(-i*phi).*(C'*x)) - alpha).^2 + ...
                    2./(beta.^2 - abs(C'*x).^2) + ...
                    4 * abs(C'*x).^2./(beta.^2 - abs(C'*x).^2).^2;

        % Function composition of the objective and its gradient and Hessian.
        f{k} = @(x) rho/2 * norm(A*x - b + u{k})^2 + ...
                    -1/t * sum(ln(real(exp(-i*phi).*(C'*x)) - alpha) + ...
                                ln(beta.^2 - abs(C'*x).^2));
        grad{k} = @(x) rho * A' * (A*x - b + u{k}) + 1/t * C * r{k}(x); 
        Hess{k} = @(x) rho * A' * A + 1/t * (C * diag(s{k}(x)) * C');

%         % Alternative function composition of gradient and Hessian using C tilde.
%         grad{k} = @(x) rho * A' * (A*x - b + u{k} + 1/(rho*t) * Ct{k} * r{k}(x));
%         Hess{k} = @(x) rho * A' * ...
%                     (eye(n) + 1/(rho*t) * (Ct{k} * diag(s{k}(x)) * Ct{k}')) * A;

        % Forms needed for efficient calculation of the Newton step.
        abridged_grad{k} = @(x) (A*x - b + u{k} + 1/(rho*t) * Ct{k} * r{k}(x));
        M{k} = @(x) eye(n) - ...
                    Ct{k} * inv((diag(rho*t./s{k}(x))) + Ct{k}'*Ct{k}) * Ct{k}';


    end

    %% Execute Newton's algorithm
    % Minimize using Newton's method.

    % Keeps track which modes have already been solved
    mode_done = false * ones(1, N);

    for j = 1 : newton_max_steps
        % Initialize the computation for the Newton step.
        unfinished_modes = find(~mode_done);
        for k = unfinished_modes
            cb{k} = invA{k}(z, -M{k}(x{k}) * abridged_grad{k}(x{k}));
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
            newton_step_error = norm(Hess{k}(x{k}) * delta_x{k} + grad{k}(x{k}));

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

        end

        vis_progress(progress); % Visualize our progress.

        % Check overall termination condition.
        if all(mode_done)
            break
        end

    end

end % End of prodQ_global function.

%% Private functions
function [z] = ln(x)
% Custom log function which returns nan for non-zero imaginary part and 
% real part <= 0.

    if ~all(isreal(x))
        error('Cannot accept complex arguments!');
    elseif any(real(x) <= 0)
        z = Inf;
    else
        z = log(x);
    end
end % End of ln function.

function default_vis_progress(progress)
% Print progress.
    % Find the maximum number of iterations present over all the modes.
    for k = 1 : length(progress)
        len(k) = length(progress{k});
    end

    fprintf('%d:', max(len)); % Print out iteration number.

    % Print out the Newton decrement for each mode.
    for k = 1 : length(progress) 
        fprintf(' %1.3e', progress{k}(end).newton_dec);
    end
    fprintf('\n');
end

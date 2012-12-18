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

        state.newton_err_thresh = 1e-6;
        state.newton_max_steps = 100;

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
    x_prev = state.x;
    u = state.u;
    newton_max_steps = state.newton_max_steps;
    newton_err_thresh = state.newton_err_thresh; 

    %% Check the gradient and Hessian of the relaxed field design objective
    % Note the special form that the Hessian must take in order to be tested.
% 
%     % Function handles for calculating the gradient and Hessian.
%     f = @(alpha, beta, C, t, phi, x) ...
%                 -sum(ln(real(exp(-i*phi).*(C'*x)) - alpha)) + ...
%                 -sum(ln(beta.^2 - abs(C'*x).^2));
%     grad_f = @(alpha, beta, C, t, phi, x) ...
%                 -sum((C * diag(exp(i*phi)./ ...
%                                 (real(exp(-i*phi).*(C'*x)) - alpha))), 2) + ...
%                 2 * sum(C * diag((C'*x)./(beta.^2 - abs(C'*x).^2)), 2);
%     Hess_f = @(alpha, beta, C, t, phi, x, dx) ...
%                 C * diag(exp(i*phi)./(real(exp(-i*phi).*(C'*x)) - alpha).^2) ...
%                     * real(exp(-i*phi) .* (C' * dx)) + ...
%                 C * diag(2./(beta.^2 - abs(C'*x).^2)) * C' * dx + ...
%                     4 * C * diag((C'*x)./(beta.^2 - abs(C'*x).^2).^2) * ...
%                     real(diag((x'*C))*C'*dx);
% 
%     % Test derivative for every mode.
%     for k = 1 : N
%         alpha = fobj(k).alpha;
%         beta = fobj(k).beta;
%         C = fobj(k).C;
%         phi = angle(C' * xv{k});
% 
%         fi = @(x) f(alpha, beta, C, t, phi, x);
%         grad_fi = @(x) grad_f(alpha, beta, C, t, phi, x);
%         Hess_fi = @(x, dx) Hess_f(alpha, beta, C, t, phi, x, dx);
% 
%         x = xv{k};
% 
%         derivative_tester(fi, grad_fi(x)', fi(x), x, 1e-6, @real);
%         hessian_tester(grad_fi, @(dx) Hess_fi(x, dx), grad_fi(x), x, 1e-6);
%     end

    %% Function handles for the gradients and Hessians.
    % Compose function handles that describe the form of the function,
    % its gradient, as well as its Hessian.

    % Barrier functions.
    f_lo = @(alpha, C, phi, x) -sum(ln(real(exp(-i*phi).*(C'*x)) - alpha));
    f_hi = @(beta, C, x) -sum(ln(beta.^2 - abs(C'*x).^2));

    % Gradient of barrier functions.
    grad_f_lo = @(alpha, C, phi, x) -sum(C * ...
                    diag(exp(i*phi)./ (real(exp(-i*phi).*(C'*x)) - alpha)), 2);
    grad_f_hi = @(beta, C, x) 2 * sum(C * ...
                    diag((C'*x)./(beta.^2 - abs(C'*x).^2)), 2);

    % Hessian of barrier functions.
    Hess_f_lo = @(alpha, C, phi, x, dx) ...
                C * diag(1./(real(exp(-i*phi).*(C'*x)) - alpha).^2) * C';
    Hess_f_hi = @(beta, C, x, dx) ...
                2 * C * diag(1./(beta.^2 - abs(C'*x).^2)) * C' + ...
                4 * C * diag(abs(C'*x).^2./(beta.^2 - abs(C'*x).^2).^2) * C';

    for k = 1 : N
        alpha = fobj(k).alpha;
        beta = fobj(k).beta;
        C = fobj(k).C;
        phi = angle(C' * x_prev{k});

        A = pres(k).A(z);
        b = pres(k).b(z);

        f{k} = @(x) 1/t * (f_lo(alpha, C, phi, x) + f_hi(beta, C, x));
        grad_f{k} = @(x) 1/t * (grad_f_lo(alpha, C, phi, x) + ...
                                grad_f_hi(beta, C, x));
        Hess_f{k} = @(x) 1/t * (Hess_f_lo(alpha, C, phi, x) + ...
                                Hess_f_hi(beta, C, x));

        
        f{k} = @(x) 1/t * (f_lo(alpha, C, phi, x) + f_hi(beta, C, x)) + ...
                    1/2 * norm(A*x - b + u{k})^2;
        grad_f{k} = @(x) 1/t * (grad_f_lo(alpha, C, phi, x) + ...
                                grad_f_hi(beta, C, x)) + ...
                        A' * (A*x - b + u{k});
        Hess_f{k} = @(x) 1/t * (Hess_f_lo(alpha, C, phi, x) + ...
                                Hess_f_hi(beta, C, x)) + ...
                        A' * A;
        
        f{k} = @(x) 1/2 * norm(A*x - b + u{k})^2 + ...
                    -1/t * sum(ln(real(exp(-i*phi).*(C'*x)) - alpha) + ...
                                ln(beta.^2 - abs(C'*x).^2));

        r{k} = @(x) -exp(i*phi)./ (real(exp(-i*phi).*(C'*x)) - alpha) + ...
                    2 * (C'*x)./(beta.^2 - abs(C'*x).^2);
        grad_f{k} = @(x) A' * (A*x - b + u{k}) + ...
                    1/t * sum(C * diag(r{k}(x)), 2);

        s{k} = @(x) 1./(real(exp(-i*phi).*(C'*x)) - alpha).^2 + ...
                    2./(beta.^2 - abs(C'*x).^2) + ...
                    4 * abs(C'*x).^2./(beta.^2 - abs(C'*x).^2).^2;
        Hess_f{k} = @(x) A' * A + ...
                    1/t * (C * diag(s{k}(x)) * C');
    end

    %% Form the problem to minimize.
    % Minimize using Newton's method.
    x = x_prev{1};
    k = 1;

    for j = 1 : newton_max_steps
        % Compute Newton step and decrement.
        delta_x = -(Hess_f{k}(x)) \ grad_f{k}(x);
        lambda = sqrt(abs(real(grad_f{k}(x)' * -delta_x)));

        fprintf('%d: %e, %e, %e\n', j, f{k}(x), norm(grad_f{k}(x)), lambda^2/2);

        % Check termination condition.
        if lambda^2/2 <= newton_err_thresh
            return
        end

        % Perform line search in Newton step direction.
        step_size = line_search_convex(f{k}, grad_f{k}, delta_x, x, 1e-9);

        % Update x.
        x = x + step_size * delta_x;
    end

end % End of prodQ_global function.

%% Other private functions
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

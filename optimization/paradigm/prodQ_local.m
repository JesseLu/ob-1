%% prodQ_local
% Produce Q(z) using a local (adjoint) optimization paradigm.

%% Description
%

function [P, q, state] = prodQ_local(z, opt_prob, state, varargin)

%% Input parameters

%% Output parameters
 
    %% Parse inputs.

    fobj = [opt_prob.field_obj];
    pres = [opt_prob.phys_res];
    invA = {opt_prob.solve_A};
    invAd = {opt_prob.solve_A_dagger};

    % Default values for the relaxed field objective.
    a = 1;
    p = 2;

    N = length(fobj); % Number of fields/modes.


    %% Solve for x_i
    % That is to say, obtain the updated field variables.

    % Initiate solves.
    for k = 1 : N
        cb{k} = invA{k}(z, pres(k).b(z));
    end

    % Complete solves.
    done = false * ones(N, 1);
    while ~all(done)
        for k = 1 : N
            [x{k}, done(k)] = cb{k}();
        end
    end

    % Compute error of x_i.
    for k = 1 : N
        err(k) = norm(pres(k).A(z) * x{k} - pres(k).b(z)) / norm(pres(k).b(z));
    end

    P = nan;
    q = nan;

    %% Compute df/dx for each mode.
    % We first use the form for the violation of the design objectives:
    %
    %
    % $$ f_v(x_i) = \mbox{neg}(|C_i^\dagger x_i - \alpha_i|) + 
    %               \mbox{neg}(|\beta_i - C_i^\dagger x_i|), $$
    %
    % where $\mbox{neg}(u) = u$ if $u < 0$ and $0$ otherwise.
    %
    % The derivative with respect to $x_i$ is then
    %
    % $$ \partial f_v(x_i) / \partial x_i = 
    %           \mbox{diag}(\tilde{a}) C_i^\dagger$$
    %
    % where the elements of vector $\tilde{a}$ are either
    % $(\pm c_{ij}^\dagger x_i)^\star / |c_{ij}^\dagger x_i|$, or 
    % $0$ respectively, 
    % depending on whether the $\alpha$, $\beta$, or neither limit is violated.
    %

    % Function handle to compute the degree that a mode currently violates
    % its design objectives.
    f_viol = @(alpha, beta, C, x) ...
                ((abs(C'*x) - alpha) < 0) .* (abs(C'*x) - alpha) + ...
                ((beta - abs(C'*x)) < 0) .* (beta - abs(C'*x));

    % Derivative of f_viol with respect to x.
    df_viol_dx = @(alpha, beta, C, x) diag( ...
                ((abs(C'*x) - alpha) < 0) .* conj(C'*x)./abs(C'*x) + ...
                ((beta - abs(C'*x)) < 0) .* conj(-C'*x)./abs(C'*x)) * C';
                
    
    % Calculate f and df/dx (for individual modes).
    for k = 1 : N
        alpha = fobj(k).alpha;
        beta = fobj(k).beta;
        C = fobj(k).C;

        f_i(k) = sum(f_viol(alpha, beta, C, x{k}).^p);
        df_i_dx{k} = sum(p * diag(f_viol(alpha, beta, C, x{k}).^(p-1)) ...
                            * df_viol_dx(alpha, beta, C, x{k}), 1);

%         % Use this to check derivative.
%         derivative_tester(@(x) sum(f_viol(alpha, beta, C, x).^p), ...
%                             df_i_dx, f_i, x{k}, 1e-6)
    end

    %% Compute grad_F

    %% Form Q(z)



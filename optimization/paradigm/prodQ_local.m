%% prodQ_local
% Produce Q(z) using a local (adjoint) optimization paradigm.

%% Description
%

function [P, q, state] = prodQ_local(z, opt_prob, state, varargin)

%% Input parameters

%% Output parameters
 
    %% Parse inputs

    % Short-hand variables from the values of the opt_prob parameter.
    fobj = [opt_prob.field_obj];
    pres = [opt_prob.phys_res];
    invA = {opt_prob.solve_A};
    invAd = {opt_prob.solve_A_dagger};

    % Determine the current state of the optimization.
    state_was_empty = false;
    if isempty(state) % No previous state, use the default state.
        state_was_empty = true;

        state.kappa = 1 / 1.1; % Default value for kappa.
        state.kappa_growth_rate = 1.1; % Percent increase for a successful step.
        state.kappa_shrink_rate = 0.5; % Percent decrease for a failed step.

        % Default values for the relaxed field objective.
        state.a = @max;
        state.p = 2;

        % Default visualization function.
        state.vis_progress = @default_vis_progress;
        
        state.f = nan;
        state.F = Inf;
        state.grad_F = nan;
        state.z = nan;
    end

    % Accept any forced parameters contained in varargin.
    for k = 2 : 2 : length(varargin)
        state = setfield(state, varargin{k-1}, varargin{k});
    end

    % Obtain state parameters.
    if state_was_empty % Allow for an initial_kappa setting.
        kappa = state.initial_kappa;
    else
        kappa = state.kappa;
    end

    kappa_growth_rate = state.kappa_growth_rate;
    kappa_shrink_rate = state.kappa_shrink_rate;
   
    a = state.a;
    p = state.p;

    prev_f = state.f;
    prev_F = state.F;
    prev_grad_F = state.grad_F;
    prev_z = state.z;

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


    %% Compute df/dx for each mode
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

        f{k} = 1 * f_viol(alpha, beta, C, x{k}).^p;
        df_dx{k} = sum((p) * diag(f_viol(alpha, beta, C, x{k}).^(p-1)) ...
                            * df_viol_dx(alpha, beta, C, x{k}), 1);
    end

    % Compute the overall value of the design objective.
    F = 0;
    for k = 1 : N
        F = F + sum(f{k});
    end
    
%     % Re-scale with respect to the a() function.
%     scale_factor = 1 / a(f);
%     for k = 1 : N
%         f{k} = scale_factor * f{k};
%         df_dx{k} = scale_factor * df_dx{k};
%     end
% 
    %% Compute grad_F
    % If previous step failed, simply use the previous value of grad_F.

    % First check if previous step succeeded or failed.

    if F <= prev_F % Previous step succeeded.
        prev_step_successful = true;

        % Update kappa and recompute grad_F.
        kappa = kappa * kappa_growth_rate;

        % Initiate A_dagger solves.
        for k = 1 : N
            cb{k} = invAd{k}(z, df_dx{k}');
        end

        % Complete A_dagger solves.
        done = false * ones(N, 1);
        while ~all(done)
            for k = 1 : N
                [d{k}, done(k)] = cb{k}(); 
            end
        end

        % Compute the gradient.
        for k = 1 : N
            df_dz(k,:) = -d{k}' * pres(k).B(x{k});
        end
        grad_F = sum(df_dz', 2);


    else % Previous step failed.
        prev_step_successful = false;

        % Update kappa, and use old value of grad_F and z.
        kappa = kappa * kappa_shrink_rate;
        f = prev_f;
        F = prev_F;
        grad_F = prev_grad_F;
        z = prev_z;

    end

%     % Use this to check grad_F.
%     function [f] = my_f(z)
%     % Private function to compute the value of f.
%         % Initiate solves.
%         for k = 1 : N
%             cb{k} = invA{k}(z, pres(k).b(z));
%         end
% 
%         % Complete solves.
%         done = false * ones(N, 1);
%         while ~all(done)
%             for k = 1 : N
%                 [x{k}, done(k)] = cb{k}();
%             end
%         end
%         
%         for k = 1 : N
%             f_indiv(k) = sum(1/a * f_viol(fobj(k).alpha, fobj(k).beta, ...
%                                     fobj(k).C, x{k}).^p);
%         end
%         f = sum(f_indiv);
%     end
%     derivative_tester(@my_f, grad_F', F, z, 1e-6, @real);



    %% Form Q(z)
    %
    % $$ Q(z) = \frac{1}{2}\|z - z_0\|^2 + \kappa \nabla_z F^\dagger (z - z_0) 
    %   = \frac{1}{2}\|z - (z_0 - \kappa \nabla_z F)\|^2 + \mbox{const.} $$

    if F ~= 0 
        fprintf(' %e', a(cell2mat(f)));
        my_vec = @(u) u(:);
        grad_F_scaled = 1/a(my_vec(cell2mat(f))) * grad_F;
        P = speye(length(grad_F_scaled));
        q = z - kappa * grad_F_scaled;
    else
        P = nan;
        q = nan;
    end

    default_vis_progress(kappa, prev_step_successful);
    
    %% Update state

    state = struct( 'kappa', kappa, ...
                    'kappa_growth_rate', kappa_growth_rate, ...
                    'kappa_shrink_rate', kappa_shrink_rate, ...
                    'a', a, ...
                    'p', p, ...
                    'F', F, ...
                    'grad_F', grad_F, ...
                    'f', {f}, ...
                    'x', {x}, ...
                    'z', z, ...
                    'prev_step_successful', prev_step_successful);

    

end % End of prodQ_local function.

function default_vis_progress(kappa, prev_step_successful)
% Print progress.
    fprintf(' [kappa: %1.2e]', kappa);
    if ~prev_step_successful
        fprintf(' [kappa_shrink]');
    end
end % End of default_vis_progress function.


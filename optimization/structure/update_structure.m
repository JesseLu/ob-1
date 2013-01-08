%% update_structure
% Updates z by solving argmin Q(z) + g(z).

%% Description
% Q(z) is described by matrix P and vector q, via Q(z) = 1/2 || Pz - q ||^2.
%
% Typically, the optimal value for z is found by solving the problem in terms
% of its parameterization variable p.
% For this reason, the problem is typically converted into "p-space".

function [z, p] = update_structure(P, q, g, p0)
    

    switch(g.scheme)
        %% Continuous case
        case 'continuous'
            p = p0;

            f_uncomp = @(p) 1/2 * norm(P * g.m(p) - q)^2 + g.w(p);

            my_compressor = @(p) ...
                (p <= g.p_range(:,1)) .* g.p_range(:,1) + ...
                (p >= g.p_range(:,2)) .* g.p_range(:,2) + ...
                ((p > g.p_range(:,1)) & p < g.p_range(:,2)) .* p;
            m = @(p) g.m(my_compressor(p));
            w = @(p) g.w(my_compressor(p));
            f = @(p) 1/2 * norm(P * m(p) - q)^2 + w(p);

            filter_grad = @(grad, p) ...
                ~(  ((p <= g.p_range(:,1)) & grad > 0) | ...
                    ((p >= g.p_range(:,2)) & grad < 0)) .* grad;

            f_prev = f(p);
            err = [];
            while true
                % Empirically find gradient.
                dp = get_gradient(@(p) f_uncomp(p), p)';

%                 % Plot error.
%                 err(end+1) = norm(filter_grad(dp, p));
%                 semilogy(err, '.-');

                % Perform line search
                optim_step = line_search_brute(f, -dp, p, 1e-6);

                if optim_step == 0 % If no step, we're done.
                    break
                end

                % Update p.
                p = my_compressor(p - optim_step * dp);
                
                % Make sure the residual never increases.
                f_curr = f(p);
                if f_curr > f_prev
                    f_curr - f_prev
                    error('Function value increased!');
                else
                    f_prev = f_curr;
                end
            end

            p = my_compressor(p);
            z = g.m(p);


        %% Continuous-linear case
        case 'continuous-linear'
            path(path, genpath(strrep(mfilename('fullpath'), ...
                                            'update_structure', 'cvx')));
            % Parameterize
            [A, b] = my_parameterize(P, q, g, p0);

            % Solve for p.
            cvx_quiet(true)
            cvx_begin
                variable dp(length(p0))
                minimize norm(A*dp - b)
                subject to
                    p0 + dp <= g.p_range(:,2)
                    p0 + dp >= g.p_range(:,1)
            cvx_end

            p = p0 + dp;
            z = g.m(p);

%             % Test
%             real(A'*(A * dp - b))
%             get_gradient(@(p) 1/2*norm(P*(g.m(p)) - q)^2 + g.w(p), p)'


        %% Discrete case
        case 'discrete'

            f = @(p) 1/2 * norm(P * g.m(p) - q)^2 + g.w(p);
            while true
                % Try all "one-off" combinations of p.
                f_min = Inf;
                for k = 1 : length(p0)
                    p = p0;
                    for l = 1 : length(g.p_range(k, :))
                        p(k) = g.p_range(k, l);
                        f_curr = f(p);
                        if f_curr < f_min
                            f_min = f_curr;
                            p_min = p;
                        end
                    end
                end

                % Test for termination condition.
                if f(p_min) < f(p0)
                    p0 = p_min;
                else
                    break
                end
            end
            p = p_min;
            z = g.m(p);

        %% Discrete-diagonal case
        case 'discrete-diagonal'
%             % Test if A is diagonal.
%             [A, b] = my_parameterize(P, q, g, p0); % Parameterize
%             A_diag = spdiags(spdiags(A, 0), 0, size(A,1), size(A,2));
%             if ~all(A == A_diag)
%                 warning('Matrix is not diagonal.');
%             end

            % Function handle to extract the individual components
            % of the structure design objective.
            my_diag = @(u) spdiags(u(:), 0, numel(u), numel(u));
            f = @(p) 1/2 * abs(P * g.m(p) - q).^2 + g.w(my_diag(p)).';

            % Try all combinations of p.
            for k = 1 : size(g.p_range, 2)
                p = g.p_range(:,k);
                res(:,k) = f(p);
            end

            % Find the optimal values of the residual.
            [~, ind] = min(res, [], 2);
            for k = 1 : length(ind)
                p(k, 1) = g.p_range(k, ind(k));
            end
            z = g.m(p);

        
        otherwise
            error('Invalid scheme for structure design objective.');
    end

end % End update_structure function.


%% Private parameterization function.
%
% We transform $Q(z) + g(z)$ into p-space via
% $$ 1/2\|P(A_m \Delta p - b_m)-q\|^2 + \mbox{real}(c_w^\dagger \Delta p) = 
%    1/2\|P A_m \Delta p - (q + P b_m - (P A_m)^{-\dagger} c_w)\|^2 + 
%       \mbox{const.}, $$
% where we have used the linearizations 
% $m(p_0 + \Delta p) \approx A_m \Delta p - b_m$ and 
% $w(p_0 + \Delta p) \approx c_w^\dagger \Delta p + \mbox{const.}$.
function [A, b] = my_parameterize(P, q, g, p0)

    % Linearization of m(p).
    A_m = get_gradient(g.m, p0);
    b_m = -g.m(p0);

    % Linearization of w(p).
    c_w = get_gradient(g.w, p0)';

    % Transform into p-space.
    A = P * A_m;
    b = q + P * b_m - (A' \ c_w);

end % End private function my_parameterize.
    

%% test_prodQ
% Test function for |prodQ_local| and |prodQ_global|.

%% Description
%

function test_prodQ(n_max, test_time)

    start_time = tic;
    N = 8;
    n = 100;
    p = 4;
    [opt_prob, x_valid] = my_create_test_problem(N, n, p);
    z = randn(n, 1) + i *randn(n, 1);
    prodQ_global(z, opt_prob, [], 'x', x_valid);
    return

    %% Simple tests
    fprintf('prodQ_local descent test:');
    while toc(start_time) < test_time
        % Create a random test problem.
        n = randi(n_max);
        N = randi(n);
        p = randi(n);
        [opt_prob] = my_create_test_problem(N, n, p);

        % Perform the local descent test from a random starting point.
        z = randn(n, 1) + i *randn(n, 1);
        test_result = test_prodQ_local_descent(z, opt_prob, ...
                                            abs(randn(1)), 1 + abs(randn(1)));

        % Check for success.
        if test_result == true
            fprintf('.');
        elseif isnan(test_result) % Test inconclusive (e.g. F = 0 initially).
            % Do nothing, just skip.
        else
            fprintf('f\n');
            warning('Test failure for parameters: N=%d, n=%d, p=%d', N, n, p);
        end
    end
    fprintf('\n');
end


%% Private function to test prodQ_local
% Tests to make sure that the design objective can actually be decreased.
% The only time that this is not possible is at the point of a local minimum with F > 0.
% With this caveat, we remark that it is highly unlikely for such a point 
% to be randomly chosen.
% For this reason, this function simply tests if F can be lowered at all
% from a random starting point.
function [success] = test_prodQ_local_descent(z, opt_prob, a, p)
    success = false;

    % Take an initial step, needed since first step is always a "success".
    [P, q, state] = prodQ_local(z, opt_prob, [], 'a', a, 'p', p);
    z = q; % Simple update.

    if state.F == 0 % If we are already at minimum, test is inconclusive.
        success = nan;
        return
    end

    state.kappa_shrink_rate = 0.1; % Make sure to shrink quickly
    for k = 1 : 12 
        [P, q, state] = prodQ_local(z, opt_prob, state); % Take a step.
        z = q; % Simple update.
        if state.prev_step_successful 
            success = true;
            break
        end
    end
end


%% Private function to create a test problem
function [opt_prob, x_valid] = my_create_test_problem(num_modes, n, p)
    z_sol = randn(n, 1) + i *randn(n, 1);
    for k = 1 : num_modes
        [fo, pr, sa, sd, x_valid{k}] = my_create_test_case(z_sol, n, p);
        opt_prob(k) = struct(   'field_obj', fo, ...
                                'phys_res', pr, ...
                                'solve_A', sa, ...
                                'solve_A_dagger', sd);
    end
end

%% Private function to create a test case
%
% A field design objective corresponding to
%
% $$ \alpha \le | C^\dagger x | \le \beta $$
%
% is created. 
% Critically, a random value for z is chosen, in order to ensure
% that a solution exists.
% In addition, a physics residual of the form
%
% $$ A_0 x - \mbox{diag}(z) x - b_0 = 0 $$
%
% where 
%
% * $A(z) = A_0 - \mbox{diag}(z)$,
% * $b(z) = b$,
% * $B(x) = -\mbox{diag}(x)$, and
% * $d(x) = A_0 x - b_0$
%
% is also created.
%
% Finally, function handles for solving A and its conjugate transpose are 
% also created.
%

function [field_obj, phys_res, solve_A, solve_A_dagger, x_valid] = ...
            my_create_test_case(z, n, p)

    % Create the physics residual.
    A0 = randn(n) + 1i * randn(n);
    b0 = randn(n, 1) + 1i * randn(n, 1);

    phys_res = struct(  'A', @(z) A0 - diag(z), ...
                        'b', @(z) b0, ...
                        'B', @(x) -diag(x), ...
                        'd', @(x) A0 * x - b0);

    % Create function handles for solving A.
    solve_A = @(z, b) get_callback(phys_res.A(z) \ b);
    solve_A_dagger = @(z, b) get_callback(phys_res.A(z)' \ b);

    % Create a random solution (x, z).
    x_valid = phys_res.A(z) \ phys_res.b(z);

    % Create the field objective.
    C = randn(n, p) + 1i * randn(n, p);
    % C = C * diag(exp(i*angle(C'*x_valid))); % Forces the initial angle to be 0.
    field_obj = struct( 'alpha', abs(C'*x_valid) * rand(1), ...
                        'beta',  abs(C'*x_valid) / rand(1), ...
                        'C', C);
end
           



%% Other private functions

function [cb] = get_callback(x)
% Return a simple (dummy) callback function.
    is_done = false;
    function [result, done] = my_simple_callback()
        if rand(1) < 0.5 || is_done
            is_done = true;
            result = x;
            done = true;  
        else
            result = nan;
            done = false;
        end
    end

    cb = @my_simple_callback;
end



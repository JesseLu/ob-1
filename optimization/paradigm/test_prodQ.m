%% test_prodQ
% Test function for |prodQ_local| and |prodQ_global|.

%% Description
%

function test_prodQ(n_max, test_time)

    start_time = tic;
%% Simple tests
    fprintf('prodQ_local descent test:');
    while toc(start_time) < test_time
        n = randi(n_max);
        N = randi(n);
        p = randi(n);
        [opt_prob] = my_create_test_problem(N, n, p);
        if test_prodQ_local_descent(randn(n, 1), opt_prob)
            fprintf('.');
        else
            fprintf('f');
            error('Test failure for parameters: N=%d, n=%d, p=%d', N, n, p);
        end
    end
    fprintf('\n');


%% Private function to test prodQ_local
% Tests to make sure that the design objective can actually be decreased.
% The only time that this is not possible is at the point of a local minimum with F > 0.
% With this caveat, we remark that it is highly unlikely for such a point 
% to be randomly chosen.
% For this reason, this function simply tests if F can be lowered at all
% from a random starting point.
function [success] = test_prodQ_local_descent(z, opt_prob)
    state = [];
    success = false;
    for k = 1 : 100
        [P, q, state] = prodQ_local(z, opt_prob, state);
        if state.prev_step_successful
            success = true;
            return
        end
    end




%% Private function to create a test problem
function [opt_prob] = my_create_test_problem(num_modes, n, p)
    z_sol = randn(n, 1);
    for k = 1 : num_modes
        [fo, pr, sa, sd] = my_create_test_case(z_sol, n, p);
        opt_prob(k) = struct(   'field_obj', fo, ...
                                'phys_res', pr, ...
                                'solve_A', sa, ...
                                'solve_A_dagger', sd);
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

function [field_obj, phys_res, solve_A, solve_A_dagger] = ...
            my_create_test_case(z, n, p)

    % Create the physics residual.
    A0 = randn(n) + 1i * randn(n);
    b0 = randn(n, 1) + 1i * randn(n, 1);

    phys_res = struct(  'A', @(z) A0 - diag(z), ...
                        'b', @(z) b0, ...
                        'B', @(x) -diag(x), ...
                        'd', @(x) A0 * x - b0);

    % Create function handles for solving A.
    solve_A = @(z, b) (@() my_callback(phys_res.A(z) \ b));
    solve_A_dagger = @(z, b) (@() my_callback(phys_res.A(z)' \ b));

    % Create a random solution (x, z).
    x = phys_res.A(z) \ phys_res.b(z);

    % Create the field objective.
    C = randn(n, p) + 1i * randn(n, p);
    C = C .* (rand(n,p) < 0.1);
    field_obj = struct( 'alpha', abs(C'*x) * rand(1), ...
                        'beta',  abs(C'*x) / rand(1), ...
                        'C', C);
           



%% Other private functions

function [x, done] = my_callback(x)
% Dummy callback function
    done = true;  

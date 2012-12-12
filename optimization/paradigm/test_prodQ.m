%% test_prodQ
% Test function for |prodQ_local| and |prodQ_global|.

%% Description
%

function [opt_prob] = test_prodQ()


%% Small test problem
    [z, opt_prob] = my_create_test_problem(2, 10, 3);
    [P, q, state] = prodQ_local(z, opt_prob, []);
    state



%% Private function to create a test problem
function [z, opt_prob] = my_create_test_problem(num_modes, n, p)
    z = randn(n, 1);
    for k = 1 : num_modes
        [fo, pr, sa, sd] = my_create_test_case(n, p);
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
% is created, as well as a physics residual of the form
%
% $$ A_0 x - \mbox{diag}(z) x - b_0 = 0 $$
%
% where 
%
% * $A(z) = A_0 - \mbox{diag}(z)$,
% * $b(z) = b$,
% * $B(x) = -\mbox{diag}(x)$, and
% * $d(x) = A_0 x - b_0$.
%
% Finally, function handles for solving A and its conjugate transpose are 
% also created.
%
function [field_obj, phys_res, solve_A, solve_A_dagger] = ...
            my_create_test_case(n, p)

    % Create the field objective.
    alpha = abs(randn(p, 1));
    field_obj = struct( 'alpha', alpha, ...
                        'beta',  alpha + abs(randn(p, 1)), ...
                        'C', randn(n, p) + 1i * randn(n, p));
    field_obj.C = field_obj.C / norm(field_obj.C);

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



%% Other private functions

function [x, done] = my_callback(x)
% Dummy callback function
    done = true;  

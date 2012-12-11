%% prodQ_local
% Produce Q(z) using a local (adjoint) optimization paradigm.

function [P, q, state] = prodQ_local(z, field_obj, phys_res, ...
                                    solve_A, solve_A_dagger, state, varargin)

%% Input parameters

%% Output parameters
 
    %% Parse inputs.

    N = length(field_obj); % Number of fields/modes.

    %% Solve for x_i
    % That is to say, obtain the updated field variables.

    % Initiate solves.
    for k = 1 : N
        cb{k} = solve_A{k}(z, phys_res(k).b(z));
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
        err(k) = norm(phys_res(k).A(z) * x{k} - phys_res(k).b(z)) ...
                    / norm(phys_res(k).b(z));
    end
    err

    P = nan;
    q = nan;

    %% Compute df_dx

    a = 1; % Normalization factor in relaxed field objective.
    p = 2; % Exponent in relaxed field objective.

    alpha = field_obj(k).alpha;
    C = field_obj(k).C;
    f_i = @(x) ((abs(C'*x) - alpha) < 0) .* (-abs(C'*x) + alpha);
    f = @(x) sum(f_i(x));

    df_dx = (- (p/a) * (-abs(C'*x{k}) + alpha).^(p-1) .* sign(C'*x{k}))' * C';

    fun = @(x) abs(C'*x);
    
    df_dx = ((C'*x{k})'./abs(C'*x{k})') * C';

    C = (C(:,:)); 
    fun = @(x) sum(abs(C'*x));
    df_dx = @(x) (conj(C'*x)./abs(C'*x)).' * C';

%     f_i = @(x) ((abs(C'*x) - alpha) < 0) .* (abs(C'*x) - alpha);
%     df_dx = @(x) diag(((abs(C'*x) - alpha) < 0) .* diag(
    fi = @(x) ((abs(C'*x) - alpha) < 0) .* (abs(C'*x) - alpha);
    dfi_dx = @(x) diag(((abs(C'*x) - alpha) < 0) .* conj(C'*x)./abs(C'*x)) * C';

    p = 1 + abs(randn(1));
    a = randn(1);
    fun = @(x) sum((1/a) * (-fi(x)).^p);
    df_dx = @(x) sum(diag((p/a) * (-fi(x)).^(p-1)) * (-dfi_dx(x)), 1);

    y = x{1};
    p
    a
    fi(y).^p
    fun(y)
    derivative_tester(fun, df_dx(y), fun(y), y, 1e-6)


    %% Compute grad_F

    %% Form Q(z)



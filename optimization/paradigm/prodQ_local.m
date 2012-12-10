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

    %% Compute grad_F

    %% Form Q(z)



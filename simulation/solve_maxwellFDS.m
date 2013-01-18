%% solve_maxwellFDS
% Solve the electromagnetic wave equation using maxwellFDS.

function [callback] = solve_maxwellFDS(omega, s_prim, s_dual, mu, epsilon, J)

    num_iters = 1e4; % Change to 1e5 later.
    err_thresh = 1e-6;
    dims = size(epsilon{1});
    E0 = {zeros(dims), zeros(dims), zeros(dims)};
    figure(5);
    cb = maxwellFDS(omega, s_prim, s_dual, mu, epsilon, E0, J, ...
                    num_iters, err_thresh);

    % Modify the callback.
    callback = @() cb(1);




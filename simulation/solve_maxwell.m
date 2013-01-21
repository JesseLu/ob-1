%% solve_maxwell
% Solve the electromagnetic wave equation using Maxwell.

function [callback] = solve_maxwell(omega, s_prim, s_dual, mu, epsilon, J)

    num_iters = 1e4; % Change to 1e5 later.
    err_thresh = 1e-6;
    dims = size(epsilon{1});
    E0 = {zeros(dims), zeros(dims), zeros(dims)};
    subplot(2, 6, 6+6);
    cb = maxwell(omega, s_prim, s_dual, mu, epsilon, E0, J, ...
                    num_iters, err_thresh);

    % Modify the callback.
    function [done, E, H] = my_callback(interval)
        start = tic;
        while toc(start) < interval
            [done, E, H] = cb();
        end
    end

    callback = @() my_callback(1);
end




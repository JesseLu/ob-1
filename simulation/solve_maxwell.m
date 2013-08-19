%% solve_maxwell
% Solve the electromagnetic wave equation using Maxwell.

function [callback] = solve_maxwell(omega, s_prim, s_dual, mu, epsilon, J)

    num_iters = 1e5;
    err_thresh = 1e-6;
    dims = size(epsilon{1});
    E0 = {zeros(dims), zeros(dims), zeros(dims)};
    subplot(2, 6, 6+6);

    if all([J{1}(:), J{2}(:), J{3}(:)] == 0) % We know the solution to this...
        cb = @() my_return_zeros(dims);
    else
        grid = struct(  'omega', omega, ...
                        's_prim', {s_prim}, ...
                        's_dual', {s_dual}, ...
                        'shape', {size(epsilon{1})}, ...
                        'origin', {[0 0 0]});
        cb = maxwell_solve_async(grid, [epsilon mu], J, ...
                                'E0', E0, ...
                                'max_iters', num_iters, ...
                                'err_thresh', err_thresh);
    end


    % Modify the callback.
    function [done, E, H] = my_callback(interval)
        start = tic;
        while toc(start) < interval
            [done, E, H] = cb();
        end
    end

    callback = @() my_callback(1);
end


function [done, E, H] = my_return_zeros(dims)
    done = true;
    E = {zeros(dims), zeros(dims), zeros(dims)};
    H = {zeros(dims), zeros(dims), zeros(dims)};
end

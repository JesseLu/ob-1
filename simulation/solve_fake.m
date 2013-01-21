%% solve_local
% Local solve of the electromagnetic wave equation.

function [callback] = solve_local(omega, s_prim, s_dual, mu, epsilon, J)

    dims = size(epsilon{1});
    E = {randn(dims), randn(dims), randn(dims)};
    H = {randn(dims), randn(dims), randn(dims)};

    %% Provide callback function
    % Use a random, artificially low "done" variable.
    done = false;
    function [done_local, E_local, H_local] = my_callback()
        if randn(1) < 1e-2
            done = true;
        end
        done_local = done;
        if done_local
            E_local = E;
            H_local = H;
        else
            E_local = {};
            H_local = {};
        end

    end

    callback = @my_callback;
end

    

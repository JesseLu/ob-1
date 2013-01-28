%% solve_local
% Local solve of the electromagnetic wave equation.

function [callback] = solve_local(omega, s_prim, s_dual, mu, epsilon, J)

    dims = size(epsilon{1});
    n = prod(dims);

    if n < 3 && ~any(dims == 1) % Don't accept 3D simulations...
        error('3D simulations not accepted.');
    end


    %% Form matrices and solve
    % We now form the necessary linear algebra components and solve the system 
    % using standard Matlab tools.

    % Get ingredient matrices and vectors.
    [A1, A2, m, e, b] = maxwell_matrices(omega, s_prim, s_dual, mu, epsilon, J); 

    % Form full matrix.
    my_diag = @(z) spdiags(z(:), 0, numel(z), numel(z));
    A = A1 * my_diag(m.^-1) * A2 - omega^2 * my_diag(e);

    % Solve
    x = A \ b;

    % Get the H-field.
    y = my_diag(1./(-i*omega*m)) * A2 * x;

    % Reshape solution.
    for k = 1 : 3
        E{k} = reshape(x((k-1)*n+1 : k*n), dims);
        H{k} = reshape(y((k-1)*n+1 : k*n), dims);
    end

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

    


function [err] = test_derivative(fun, df_dz, f0, z0, step_len)
% Check a derivative.
    
    % Produce a random direction.
    dz = randn(size(z0));
    dz = step_len * dz / norm(dz);

    % Evaluate delta in that direction empirically
    delta_empirical = fun(z0 + dz) - f0;
    delta_derivative = real(df_dz * dz);

    err = norm(delta_empirical - delta_derivative) / norm(delta_empirical);

    fprintf('Percent error in derivative: %e\n', err);

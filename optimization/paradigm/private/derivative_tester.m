
function [err] = derivative_tester(fun, df_dz, f0, z0, step_len, varargin)
% Check a derivative.
    
    if ~isempty(varargin)
        optional_fun = varargin{1};
    else
        optional_fun = @(z) z;
    end

    % Produce a random direction.
    dz = randn(size(z0)) + 1i * randn(size(z0));
    dz = step_len * dz / norm(dz);

    % Evaluate delta in that direction empirically
    delta_empirical = fun(z0 + dz) - f0;
    delta_derivative = optional_fun(df_dz * dz);

    err = norm(delta_empirical - delta_derivative) / norm(delta_empirical);

    fprintf('Percent error in derivative: %e', err);

    if norm(delta_empirical) == 0
        fprintf(' [no empirical change]');
    end

    fprintf('\n');

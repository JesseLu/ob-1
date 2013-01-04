%% translation_layer
% Takes the description of the modes and casts it into the language 
% of linear algebra which can then be optimized in the optimization module.

%% Description
%

function [opt_prob] = translation_layer(modes, solver)

%% Input parameters
% * |modes| is a structure array describing the optimization in physical
%   terms.
% * |solver| is a function handle which (asynchronously) produces the 
%   the result to the electromagnetic wave equation.

%% Output parameters
% * |opt_prob| is the the structure array describing the optimization
%   in linear algebra terms.

    for i = 1 : length(modes)
        opt_prob(i) = translate_mode(modes(i), solver);
    end

end % End of translation_layer function.

%% translate_mode private function
% Translates a single mode.
function [opt_prob] = translate_mode(mode, solver)

    % Allows us to stop writing mode.* everywhere...
    omega = mode.omega;
    s_prim = mode.s_prim;
    s_dual = mode.s_dual;
    mu = mode.mu;
    epsilon = mode.epsilon_const;
    S = mode.S;
    in = mode.in;
    outs = mode.out;

    % Useful helper functions.
    dims = size(epsilon{1});
    n = prod(dims);
    vec = @(f) [f{1}(:); f{2}(:); f{3}(:)];
    unvec = @(v) {reshape(v(1:n), dims), ...
                    reshape(v(n+1:2*n), dims), ...
                    reshape(v(2*n+1:3*n), dims)};
    my_diag = @(u) spdiags(u(:), 0, numel(u), numel(u));

    %% 
    % Compose the physics residual.

    % Get the excitation for the input mode.
    [beta, E, H, J] = solve_waveguide_mode(omega, s_prim, s_dual, ...
                                            mu, epsilon, ...
                                            in.pos, in.dir, in.mode_num);
    % Translate to linear algebra.
    [A1, A2, m, e, b] = maxwell_matrices(omega, s_prim, s_dual, ...
                                            mu, epsilon, J);

    b = sqrt(in.power) * b; % Scale the input excitation.

    phys_res = struct(  'A', @(z) A1 * my_diag(1./m) * A2 - ...
                                    omega^2 * my_diag(e + S * z), ...
                        'b', @(z) b, ...
                        'B', @(x) - omega^2 * my_diag(x) * S, ...
                        'd', @(x) b - (A1 * my_diag(1./m) * A2 - ...
                                    omega^2 * my_diag(e)) * x);


    %%
    % Compute the field design objective.
    
    for j = 1 : length(outs)
        out = outs(j);

        % Find the field pattern of the desired output mode.
        [~, E] = solve_waveguide_mode(omega, s_prim, s_dual, ...
                                       mu, epsilon, ...
                                       out.pos, out.dir, out.mode_num);
        C(:,j) = vec(E); % Vectorize
        alpha(j,1) = sqrt(out.power(1)) * norm(C(:,j))^2; % Scale.
        beta(j,1) = sqrt(out.power(2)) * norm(C(:,j))^2;
    end

    field_obj = struct( 'alpha', alpha, ...
                        'beta', beta, ... 
                        'C', C);

    %% 
    % Create function handles for solving for A and its conjugate transpose.

    solve_A = @(z, b) my_solve_A(solver, unvec, ...
                                omega, s_prim, s_dual, m, e + S * z, b);
    solve_A_dagger = @(z, b) my_solve_A_dagger(solver, unvec, ...
                                omega, s_prim, s_dual, m, e + S * z, b);


    %%
    % And now we can form the optimization problem structure.
    opt_prob = struct(  'field_obj', field_obj, ...
                        'phys_res', phys_res, ...
                        'solve_A', solve_A, ...
                        'solve_A_dagger', solve_A_dagger);



%     %% 
%     % This can be uncommented to test
%     % Test the physics residual.
%     n = size(S, 1);
%     l = size(S, 2);
%     zt = randn(l, 1) + 1i * randn(l, 1);
%     xt = randn(n, 1) + 1i * randn(n, 1);
%     pr1 = phys_res.A(zt) * xt - phys_res.b(zt);
%     pr2 = phys_res.B(xt) * zt - phys_res.d(xt);
%     norm(pr1 - pr2)
% 
% 
%     % Test solve.
%     cb = solve_A(zt, xt);
%     while ~cb()
%     end
%     [~, x] = cb();
%     norm(phys_res.A(zt) * x - xt)
% 
%     cb = solve_A_dagger(zt, xt);
%     while ~cb()
%     end
%     [~, x] = cb();
%     norm(phys_res.A(zt)' * x - xt)
% 
end % End translate_mode private function.

%% Private functions for solving A and A_dagger.
function [cb] = my_solve_A(solver, unvec, omega, s_prim, s_dual, m, e, b)
    b = b ./ (-1i * omega); % Transform into J.
    cb_orig = solver(omega, s_prim, s_dual, unvec(m), unvec(e), unvec(b));
    cb = @() my_solve_A_callback(cb_orig);
end 

function [x, done] = my_solve_A_callback(cb_orig)
    [done, E, H] = cb_orig();
    if done
        x = [E{1}(:); E{2}(:); E{3}(:)];
    else
        x = nan;
    end
end


function [cb] = my_solve_A_dagger(solver, unvec, ...
                                    omega, s_prim, s_dual, m, e, b)
%%
% This works via 
% $ A^\dagger y = A^\dagger S^\ast A^{-\ast} S^{-\ast} A^\dagger y $
% where the diagonal symmetrization matrix $S$ has the property 
% $ SA = A^T S$.
%

    % Form elements of diagonal symmetrization matrix S.
    [spx, spy, spz] = ndgrid(s_prim{1}, s_prim{2}, s_prim{3});
    [sdx, sdy, sdz] = ndgrid(s_dual{1}, s_dual{2}, s_dual{3});

    s = [sdx(:).*spy(:).*spz(:); ...
        spx(:).*sdy(:).*spz(:); ...
        spx(:).*spy(:).*sdz(:)];
     
    % Conjugate all the elements which compose the A matrix.
    for k = 1 : 3
        s_prim{k} = conj(s_prim{k});
        s_dual{k} = conj(s_dual{k});
    end

    omega = conj(omega);
    m = conj(m);
    e = conj(e);

    % Invert b by the conjugate of s.
    b = (b ./ (-1i * omega)) ./ conj(s);

    % Because of the conjugated terms, this actually solves A conjugate.
    cb_orig = solver(omega, s_prim, s_dual, unvec(m), unvec(e), unvec(b));

    % Special callback...
    cb = @() my_solve_A_dagger_callback(cb_orig, s);

end 

function [x, done] = my_solve_A_dagger_callback(cb_orig, s)
    [done, E, H] = cb_orig();
    if done
        x = conj(s) .* ([E{1}(:); E{2}(:); E{3}(:)]); % Re-transform by S.
    else
        x = nan;
    end
end

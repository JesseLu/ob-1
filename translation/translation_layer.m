function [opt_prob] = translation_layer(modes, solver)

    N = length(modes);

    for i = 1 : N
        opt_prob(i) = translate_mode(modes(i), solver);
    end

end % End of translation_layer function.

function [opt_prob] = translate_mode(mode, solver)

    omega = mode.omega;
    s_prim = mode.s_prim;
    s_dual = mode.s_dual;
    mu = mode.mu;
    epsilon = mode.epsilon_const;
    S = mode.S;
    in = mode.in;
    outs = mode.out;

    dims = size(epsilon{1});
    n = prod(dims);
    vec = @(f) [f{1}(:); f{2}(:); f{3}(:)];
    unvec = @(v) {reshape(v(1:n), dims), ...
                    reshape(v(n+1:2*n), dims), ...
                    reshape(v(2*n+1:3*n), dims)};
    my_diag = @(u) spdiags(u(:), 0, numel(u), numel(u));

    %% Compose the physics residual.
    [beta, E, H, J] = solve_waveguide_mode(omega, s_prim, s_dual, ...
                                            mu, epsilon, ...
                                            in.pos, in.dir, in.mode_num);
    [A1, A2, m, e, b] = maxwell_matrices(omega, s_prim, s_dual, ...
                                            mu, epsilon, J);
    b = sqrt(in.power) * b;

    phys_res = struct(  'A', @(z) A1 * my_diag(1./m) * A2 - ...
                                    omega^2 * my_diag(e + S * z), ...
                        'b', @(z) b, ...
                        'B', @(x) - omega^2 * my_diag(x) * S, ...
                        'd', @(x) b - (A1 * my_diag(1./m) * A2 - ...
                                    omega^2 * my_diag(e)) * x);

    % Test the physics residual.
    n = size(S, 1);
    l = size(S, 2);
    zt = randn(l, 1) + 1i * randn(l, 1);
    xt = randn(n, 1) + 1i * randn(n, 1);
    pr1 = phys_res.A(zt) * xt - phys_res.b(zt);
    pr2 = phys_res.B(xt) * zt - phys_res.d(xt);
    norm(pr1 - pr2)

    %% Compute the field design objective
    for j = 1 : length(outs)
        out = outs(j);
        [~, E] = solve_waveguide_mode(omega, s_prim, s_dual, ...
                                       mu, epsilon, ...
                                       out.pos, out.dir, out.mode_num);
        C(:,j) = vec(E);

        alpha(j,1) = sqrt(out.power(1)) * norm(C(:,j))^2;
        beta(j,1) = sqrt(out.power(2)) * norm(C(:,j))^2;
    end
    field_obj = struct( 'alpha', alpha, ...
                        'beta', beta, ... 
                        'C', C);

    %% Create function handles for solving A
    solve_A = @(z, b) solver(omega, s_prim, s_dual, mu, ...
                            unvec(e + S * z), ...
                            unvec(b./(-1i*omega)));
    solve_A = @(z, b) my_solve_A(solver, unvec, ...
                                omega, s_prim, s_dual, m, e + S * z, b);
    solve_A_dagger = @(z, b) my_solve_A_dagger(solver, unvec, ...
                                omega, s_prim, s_dual, m, e + S * z, b);

    % Test solve.
    cb = solve_A(zt, xt);
    while ~cb()
    end
    [~, x] = cb();
    norm(phys_res.A(zt) * x - xt)

    cb = solve_A_dagger(zt, xt);
    while ~cb()
    end
    [~, x] = cb();
    norm(phys_res.A(zt)' * x - xt)


    opt_prob = struct('field_obj', field_obj, ...
                        'phys_res', phys_res, ...
                        'solve_A', solve_A, ...
                        'solve_A_dagger', solve_A_dagger);

end % End translate_mode private function.

function [cb] = my_solve_A(solver, unvec, omega, s_prim, s_dual, m, e, b)
    b = b ./ (-1i * omega);
    cb_orig = solver(omega, s_prim, s_dual, unvec(m), unvec(e), unvec(b));
    cb = @() my_solve_A_callback(cb_orig);
end 

function [done, x] = my_solve_A_callback(cb_orig)
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
     
    for k = 1 : 3
        s_prim{k} = conj(s_prim{k});
        s_dual{k} = conj(s_dual{k});
    end

    omega = conj(omega);
    m = conj(m);
    e = conj(e);
    b = (b ./ (-1i * omega)) ./ conj(s);

    cb_orig = solver(omega, s_prim, s_dual, unvec(m), unvec(e), unvec(b));
    cb = @() my_solve_A_dagger_callback(cb_orig, s);

end 

function [done, x] = my_solve_A_dagger_callback(cb_orig, s)
    [done, E, H] = cb_orig();
    if done
        x = conj(s) .* ([E{1}(:); E{2}(:); E{3}(:)]);
    else
        x = nan;
    end
end

function [opt_prob] = translation_layer(modes, solver)

    N = length(modes);
    vec = @(f) [f{1}(:); f{2}(:); f{3}(:)];
    unvec = @(v, dims) {reshape(v(1 : prod(dims)), dims), ...
                        reshape(v(prod(dims)+1 : 2*prod(dims)), dims), ...
                        reshape(v(2*prod(dims)+1 : 3*prod(dims)), dims)};
    my_diag = @(u) spdiags(u(:), 0, numel(u), numel(u));

    for i = 1 : N
% 
%     % Create the physics residual.
%     A0 = randn(n) + 1i * randn(n);
%     b0 = randn(n, 1) + 1i * randn(n, 1);
% 
%     phys_res = struct(  'A', @(z) A0 - diag(z), ...
%                         'b', @(z) b0, ...
%                         'B', @(x) -diag(x), ...
%                         'd', @(x) A0 * x - b0);
% 
%     % Create function handles for solving A.
%     solve_A = @(z, b) get_callback(phys_res.A(z) \ b);
%     solve_A_dagger = @(z, b) get_callback(phys_res.A(z)' \ b);
% 
%     % Create a random solution (x, z).
%     x_valid = phys_res.A(z) \ phys_res.b(z);
% 
%     % Create the field objective.
%     C = randn(n, p) + 1i * randn(n, p);
%     % C = C * diag(exp(i*angle(C'*x_valid))); % Forces the initial angle to be 0.
%     field_obj = struct( 'alpha', abs(C'*x_valid) * rand(1), ...
%                         'beta',  abs(C'*x_valid) / rand(1), ...
%                         'C', C);
        m = modes(i);
        dims = size(m.epsilon_const{1});
        %% Compose the physics residual.
        [beta, E, H, J] = solve_waveguide_mode(m.omega, m.s_prim, m.s_dual, ...
                                                m.mu, m.epsilon_const, ...
                                                m.in.pos, m.in.dir, m.in.mode_num);
        [A1, A2, mu, eps, b] = maxwell_matrices(m.omega, m.s_prim, m.s_dual, ...
                                                m.mu, m.epsilon_const, J);
        b = sqrt(m.in.power) * b;

        phys_res = struct(  'A', @(z) A1 * my_diag(1./mu) * A2 - ...
                                        m.omega^2 * my_diag(eps + m.S * z), ...
                            'b', @(z) b, ...
                            'B', @(x) - m.omega^2 * my_diag(x) * m.S, ...
                            'd', @(x) b - (A1 * my_diag(1./mu) * A2 - ...
                                        m.omega^2 * my_diag(eps)) * x);

        % Test the physics residual.
        n = size(m.S, 1);
        l = size(m.S, 2);
        zt = randn(l, 1) + 1i * randn(l, 1);
        xt = randn(n, 1) + 1i * randn(n, 1);
        pr1 = phys_res.A(zt) * xt - phys_res.b(zt);
        pr2 = phys_res.B(xt) * zt - phys_res.d(xt);
        norm(pr1 - pr2)

        %% Compute the output excitations.
        for j = 1 : length(m.out)
            mo = m.out(j);
            [~, E] = solve_waveguide_mode(m.omega, m.s_prim, m.s_dual, ...
                                                    m.mu, m.epsilon_const, ...
                                                    mo.pos, mo.dir, mo.mode_num);
            C(:,j) = vec(E);

            alpha(j,1) = sqrt(mo.power(1)) * norm(C(:,j))^2;
            beta(j,1) = sqrt(mo.power(2)) * norm(C(:,j))^2;
        end
        field_obj = struct( 'alpha', alpha, ...
                            'beta', beta, ... 
                            'C', C);

        %% Create function handles for solving A
        solve_A = @(z, b) solver(m.omega, m.s_prim, m.s_dual, m.mu, ...
                                unvec(vec(m.epsilon_const) + m.S * z, dims), ...
                                unvec(b./(i*m.omega), dims));
        % solve_A_dagger = @(z, b) get_callback(phys_res.A(z)' \ b);

        % Test solve.
        cb = solve_A(zt, xt);
        while ~cb()
        end
        [~, E] = cb();
        x = vec(E);
        norm(phys_res.A(zt) * x - xt)

        opt_prob(i) = struct('field_obj', field_obj, ...
                            'phys_res', phys_res, ...
                            'solve_A', solve_A, ...
                            'solve_A_dagger', nan);


    end


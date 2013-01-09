function example()
    path(path, genpath('.'));


    %% Hard-coded constants.
    omega = 0.12;
    dims = [80 80 1];
    z_thickness = 10;
    z_center = dims(3)/2;
    eps_lo = 1.5;
    eps_hi = 13;


    %% Build up the base structure.

    mu = {ones(dims), ones(dims), ones(dims)};

    epsilon = {eps_lo*ones(dims), eps_lo*ones(dims), eps_lo*ones(dims)};

    my_shapes = {struct('type', 'rectangle', ...
                        'position', [0 0], ...
                        'size', [1e9 1e9], ...
                        'permittivity', eps_lo), ...
                struct('type', 'rectangle', ...
                        'position', [0 0], ...
                        'size', [1e9 12], ...
                        'permittivity', eps_hi)};

    epsilon_0 = add_planar(epsilon, z_center, z_thickness, my_shapes); 

    [s_prim, s_dual] = stretched_coordinates(omega, dims, [10 10 0]);

    % Build the selection matrix, and reset values of epsilon.
    [S, epsilon] = planar_selection_matrix('alternate', epsilon_0, ...
                                        {[21 21], [60 60]}, eps_lo, ...
                                        z_center, z_thickness);


    %% Specify structure design objective 
    % Otherwise known as the parameterization of z.
    struct_obj = struct('m', @(p) p, ...
                        'w', @(p) 0, ...
                        'p_range', ones(size(S,2), 1) * [0 1], ...
                        'scheme', 'discrete-diagonal');


    %% Specify modes

    % Allow for a simple way to specify a waveguide mode.
    wg = @(power, xpos, dir, mode_num) ...
                struct('type', 'wgmode', ...
                    'power', power, ...
                    'pos', {{[xpos 1 1], [xpos dims(2) dims(3)]}}, ...
                    'dir', dir, ...
                    'mode_num', mode_num);

    modes(1) = struct('omega', omega, ...
                    'in', wg(1, 15, 'x+', 1), ...
                    'out', [wg([0.99 1], 68, 'x+', 1), ...
                            wg([0 0.1], 68, 'x+', 3), ...
                            wg([0 0.1], 12, 'x-', 1)], ...
                    's_prim', {s_prim}, ...
                    's_dual', {s_dual}, ...
                    'mu', {mu}, ...
                    'epsilon_const', {epsilon}, ...
                    'S', (eps_hi - eps_lo) * S);


    %% Construct the optimization problem
    [opt_prob, J, E_out] = translation_layer(modes, @solve_local);
    % test_opt_prob(opt_prob, S); % Use to test opt_prob.


    %% Optimize
    p0 = struct_obj.p_range(:,2);
    [z, p, state] = run_optimization(opt_prob, struct_obj, p0, 'local');
    [z, p, state] = run_optimization(opt_prob, struct_obj, p0, 'global');

    %% Verify
    % modes = verification_layer(opt_prob, z);
    modes = verification_layer(opt_prob, z, state.x);

    fprintf('     alpha    result      beta\n');
    disp(modes(1).output_power)

    % Visualize.
    subplot 121; imagesc(modes(1).epsilon{2}'); axis equal tight;
    subplot 122; imagesc(abs(modes(1).E{3})'); axis equal tight;

%     %% Visualize and test.
%     cb = solve_local(omega, s_prim, s_dual, mu, epsilon_0, J{1})
%     while ~cb()
%     end
%     [~, E] = cb();
%     subplot 121; imagesc(real(E{3})');
%     subplot 122; imagesc(real(E_out{1}{2}{3})');
%     (abs(dot(E{3}(:), E_out{1}{1}{3}(:))) / norm(E_out{1}{1}{3}(:))^2)^2
%     (abs(dot(E{3}(:), E_out{1}{2}{3}(:))) / norm(E_out{1}{2}{3}(:))^2)^2
%     for k = 1 : 3
%         subplot(1, 3, k);
%         imagesc(abs(E{k})');
%         axis equal tight;
%     end


function test_opt_prob(opt_prob, S)

    % Test opt_prob.
    n = size(S, 1);
    l = size(S, 2);
    x = randn(n, 1) + 1i * randn(n, 1);
    z = randn(l, 1) + 1i * randn(l, 1);

    for i = 1 : length(opt_prob)
        pr = opt_prob(i).phys_res;
        phys_res_error(i) = norm(pr.A(z)*x-pr.b(z) - ...
                                    (pr.B(x)*z-pr.d(x)));  

        cb1 = opt_prob(i).solve_A(z, x);
        cb2 = opt_prob(i).solve_A_dagger(z, x);
        done = [false false];
        while ~all(done)
            [~, done(1)] = cb1();
            [~, done(2)] = cb2();
        end
        solve_A_error(i) = norm(pr.A(z) * cb1() - x);
        solve_A_dagger_error(i) = norm(pr.A(z)' * cb2() - x);
        fprintf('opt_prob_test [pr: %e, sA: %e, sAd: %e]\n', ...
                                            phys_res_error(i), ...
                                            solve_A_error(i), ...
                                            solve_A_dagger_error(i));
    end


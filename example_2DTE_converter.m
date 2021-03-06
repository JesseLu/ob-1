%% example_2DTE_converter
% A simple example of an optimization (2D).

%% Description
% Optimization, in 2D TE mode, of a waveguide mode converter.

function [vis_result] = example(paradigm, S_type, update_scheme, ...
                                        num_iters, err_thresh, varargin)
    help example_2DTE_converter

    %% Source code
    path(path, genpath('.'));
    %% Hard-coded constants.
    omega = 0.16;
    dims = [60 60 1];
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
                        'size', [1e9 13], ...
                        'permittivity', eps_hi)};

    epsilon_0 = add_planar(epsilon, z_center, z_thickness, my_shapes); 

    [s_prim, s_dual] = stretched_coordinates(omega, dims, [10 10 0]);

    % Build the selection matrix, and reset values of epsilon.
    [S, epsilon] = planar_selection_matrix(S_type, epsilon_0, ...
                                        {[12 24], [48 37]}, eps_lo, ...
                                        z_center, z_thickness);


    %% Specify structure design objective 
    % Otherwise known as the parameterization of z.
    struct_obj = struct('m', @(p) p, ...
                        'w', @(p) 0e-5 * sum(-p(:)), ...
                        'p_range', ones(size(S,2), 1) * [0 1], ...
                        'scheme', update_scheme);


    %% Specify modes

    % Allow for a simple way to specify a waveguide mode.
    wg = @(power, xpos, dir, mode_num) ...
                struct('type', 'wgmode', ...
                    'power', power, ...
                    'pos', {{[xpos 1 1], [xpos dims(2) dims(3)]}}, ...
                    'dir', dir, ...
                    'mode_num', mode_num);

    delta = 1e-2;
    modes(1) = struct('omega', omega, ...
                    'in', wg(1, 11, 'x+', 2), ...
                    'out', [wg([1-delta 2], 49, 'x+', 4), ...
                            wg([0 delta], 49, 'x+', 2)], ...
                    's_prim', {s_prim}, ...
                    's_dual', {s_dual}, ...
                    'mu', {mu}, ...
                    'epsilon_const', {epsilon}, ...
                    'S', (eps_hi - eps_lo) * S);

%     modes(2) = struct('omega', omega, ...
%                     'in', wg(1, 11, 'x+', 4), ...
%                     'out', [wg([0 delta], 69, 'x+', 4)], ...
%                     's_prim', {s_prim}, ...
%                     's_dual', {s_dual}, ...
%                     'mu', {mu}, ...
%                     'epsilon_const', {epsilon}, ...
%                     'S', (eps_hi - eps_lo) * S);


    %% Translate
    [opt_prob, J, E_out] = translation_layer(modes, @solve_local);
    % test_opt_prob(opt_prob, S); % Use to test opt_prob.

    %% Optimize
    p0 = struct_obj.p_range(:,2);
    mode_sel = 1 : length(modes);
    vis_layer = struct( 'component', 2, ...
                        'slice_dir', 'z', ...
                        'slice_index', round(dims(3)/2));

    options = struct(   'paradigm', paradigm, ...
                        'starting_iter', 1, ...
                        'num_iters', num_iters, ...
                        'err_thresh', err_thresh, ...
                        'paradigm_args', {{'t', 1e6, ...
                            'newton_err_thresh', 1e-3, ...
                            'vis_progress', ...
                            @(progress) vis_newton_progress(progress, mode_sel)}}, ...
                        'structure_args', {{}}, ...
                        'state_file', 'ex2D_state.mat', ...
                        'history_file', 'ex2D_hist.h5', ...
                        'vis_progress', @(k, x, z, p, progress) ...
                            track_progress(opt_prob, struct_obj, vis_layer, ...
                                            mode_sel, k, x, z, p, progress) ...
                        );

    if ~isempty(varargin)
        opt_state = load(varargin{1});
        opt_prob = opt_state.opt_prob;
        struct_obj = opt_state.g;
        p0 = opt_state.p;
        options = opt_state.options;
        options.starting_iter = opt_state.k + 1;
        options.num_iters = num_iters;
        state = opt_state.state;
        progress = opt_state.progress;
    else
        state = [];
        progress = [];
    end

    [z, p, state] = run_optimization(opt_prob, struct_obj, p0, options, ...
                                        state, progress);


    %% Visualize.
    vis_result = @() track_progress(opt_prob, struct_obj, vis_layer, ...
                                    mode_sel, [], [], z, p);

end % End example function.


function test_opt_prob(opt_prob, S)
% Tests the opt_prob structure.

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
end % End  test_opt_prob private function.

function [opt_prob] = translation_layer(z_param, modes, solver)

    opt_prob = nan;
    N = length(modes)

    for i = 1 : N
        m = modes(i);

        %% Compute the input excitation.
        m.in
        [beta, E, H, J] = solve_waveguide_mode(m.omega, m.s_prim, m.s_dual, ...
                                                m.mu, m.epsilon_const, ...
                                                m.in.pos, m.in.dir, m.in.mode_num);

        cb{i} = solver(m.omega, m.s_prim, m.s_dual, m.mu, m.epsilon_const, J);
    end

    done = false * ones(N, 1);
    while ~all(done)
        for i = 1 : N
            [done(i), E, H] = cb{i}();
        end
    end

    dims = size(E{1});
    for k = 1 : 3
        subplot(2, 3, k);
        imagesc(real(E{k}(:,:,1))'); axis equal tight;
        subplot(2, 3, k+3);
        imagesc(real(H{k}(:,:,1))'); axis equal tight;
    end



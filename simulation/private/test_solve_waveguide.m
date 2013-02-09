function test_solve_waveguide(x_start, order, varargin)

    dims = [160 60 1];
    omega = 0.16;
    z_thickness = 10;
    z_center = dims(3)/2;
    eps_lo = 1.5;
    eps_hi = 13;


    %% Build up the base structure.


    if ~isempty(varargin)
        epsilon = varargin{1};
        dims = size(epsilon{1});
        % ylims = round(dims(2)/2 + 30 * [-1 1]);
        for i = 1 : 3
            % epsilon{i} = epsilon{i}(160:end,ylims(1):ylims(2),:);
            epsilon{i} = interp2(epsilon{i});
        end
        omega = omega/2;
        dims = size(epsilon{1});
        if numel(dims) == 2
            dims = [dims, 1];
        end
    else
        dims = [160 60 1];
        epsilon = {eps_lo*ones(dims), eps_lo*ones(dims), eps_lo*ones(dims)};

        my_shapes = {struct('type', 'rectangle', ...
                            'position', [0 0], ...
                            'size', [1e9 1e9], ...
                            'permittivity', eps_lo), ...
                    struct('type', 'rectangle', ...
                            'position', [0 0], ...
                            'size', [1e9 13], ...
                            'permittivity', eps_hi), ...
                    struct('type', 'rectangle', ...
                            'position', [0 3], ...
                            'size', [2 2], ...
                            'permittivity', eps_hi)};

        epsilon = add_planar(epsilon, z_center, z_thickness, my_shapes); 
    end

    mu = {ones(dims), ones(dims), ones(dims)};
    [s_prim, s_dual] = stretched_coordinates(omega, dims, [10 10 0]);

    my_plot(@real, epsilon);

    [~, ~, ~, J] = solve_waveguide_mode(...
                            omega, s_prim, s_dual, mu, epsilon, ...
                            {[x_start 1 1], [x_start dims(2:3)]}, 'x+', 2);

    [beta, E_wg, ~, ~] = solve_waveguide_mode(...
                            omega, s_prim, s_dual, mu, epsilon, ...
                            {[x_start 1 1], [x_start dims(2:3)]}, 'x+', order);

    for i = 1 : 3
        J{i} = J{i} / 2;
    end
    cb = solve_local(omega, s_prim, s_dual, mu, epsilon, J);
    while ~cb(); end
    [~, E, H] = cb();

    my_plot(@abs, E); 

    for i = 1 : 3
        e{i} = E_wg{i}(x_start,:,:);
    end
    e{1} = -1 * e{1} * exp(-1i * beta * 0.5);

    % my_plot(@real, e); 

    for i = 1 : dims(1)
        c(i) = my_dot(e, E, i);
    end
    c = c(:);

    subplot 212;
    plot([real(c), imag(c), abs(c)], '.-');

    % plot(abs(c(95:150)), '.-');
end

function [res] = my_dot(wg_field, vec_field, ind)
    for i = 1 : 3
        a = wg_field{i};
        b = vec_field{i}(ind,:,:);
        res(i) = dot(a(:), b(:));
        anorm(i) = norm(a)^2;
    end
    res = sum(res) / sum(anorm);
end

function my_plot(my_fun, vec_field)
  for i = 1 : 3
        subplot(2, 3, i);
        imagesc(my_fun(vec_field{i})');
        axis equal tight;
    end
end


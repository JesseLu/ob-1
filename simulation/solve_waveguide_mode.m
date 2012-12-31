%% solve_waveguide_mode
% Find the mode of a waveguide, as well as the current excitation for it.

%% Description
% Computes the propagation mode for a nanophotonic waveguide structure
% including the wave-vector, E- and H-fields, as well as the current excitation
% needed for omni- or uni-directional excitation.

function [beta, E, H, J] = solve_waveguide_mode(...
                                        omega, s_prim, s_dual, mu, epsilon, ...
                                        pos, dir, mode_num)

%% Input parameters
% The input parameters are very similar to those which describe a simulation,
% with the exception that most of the parameters are in two-dimensions (x and y)
% only.
%
% Additionally, parameters describing the location, direction, and order of
% the waveguide mode are included.
%
% * |omega|, |s_prim|, |s_dual|, |mu|, and |epsilon| should be identical
%   to the values used to desribe any simulation.
% * |pos| is a cell array of 2 three-element vectors describing the bounded
%   plane on which to excite the waveguide mode. 
%   Specifically, |pos| should look like |{[x0 y0 z0], [x1 y1 z1]}|.
%   Note that if propagation in the x-direction is desired, then |x0| should
%   equal |x1|.
% * |dir| is a string denoting the diretion of propagation for the waveguide.
%   Possible values include |'x'|, |'y'|, and |'z'| which denote propagation
%   in both the positive and negative directions as well as unidirectional 
%   excitation directions |'x+'|, |'x-'|, |'y+'|, |'y-'|, |'z+'|, and |'z-'|.
% * |mode_num| is the order of the mode to compute where |1| denotes the
%   fundamental mode, |2| denotes the second order mode and so on.

%% Output parameters
% * |beta| is the wavevector of the mode.
% * |E| and |H| are the E- and H-fields of the mode. 
%   Note that the size of each component array is matched to the bounded
%   plane specified by the |pos| input parameter.
% * |J| is the current excitation needed for the mode.
%   Unlike |E| and |H|, |J| spans the entire simulation space and includes
%   a plan in-front of the bounded plane in order to enable a unidirectional
%   source.

    %% Parse inputs

    % Shorthand for the bounded plane.
    p0 = pos{1};
    p1 = pos{2};
    shape = p1 - p0 + 1;

    % Cut out the bounded plane.
    for k = 1 : 3
        sp{k} = s_prim{k}(p0(k):p1(k));
        sd{k} = s_dual{k}(p0(k):p1(k));
        eps{k} = epsilon{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
        m{k} = mu{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
    end

    % Figure out what direction we are to propagate in.
    if all(dir(1) ~= 'xyz')
        error('The propagation direction must be either x, y, or z.');  
    end
    prop_dir = find(dir(1) == 'xyz');

    
    %% Build the operator
    % Build both real-only and full-complex versions of the operator.

    % Full complex operator.
    [A, get_wg_fields] = wg_operator(omega, sp, sd, eps, m, prop_dir, shape);

    % Real-only operator.
    for k = 1 : 3
        sp_r{k} = real(sp{k});
        sd_r{k} = real(sd{k});
        eps_r{k} = real(eps{k});
        m_r{k} = real(m{k});
    end
    A_r = wg_operator(real(omega), sp_r, sd_r, eps_r, m_r, prop_dir, shape);


    %% Solve for largest-magnitude eigenvalue of the real operator 
    % This is done in order to obtain the appropriate shift, 
    % from which we can calculate the most negative eigenvalues.

    % Use the power iteration algorithm.
    n = size(A_r, 1);
    v = randn(n, 1);
    for k = 1 : 20 % 20 iterations should always be enough for an estimate.
        v = A_r * v;
    end
    ev_max = (v' * A_r * v) / norm(v)^2; % Rayleigh quotient.
    shift = abs(ev_max); % Shift works for both positive and negative ev_max.


    %% Solve for the desired eigenvector of the real operator
    % Taking the real operator, we a few of the most negative eigenmodes,
    % and then choose the one we are interested in.

    % Shift the matrix and find the appropriate eigenmode.
    % Find a few extra modes just to be sure we found the correct one.
    [V, D] = eigs(A_r - shift * speye(n), mode_num + 2); 
    
    gamma = diag(D);
    [temp, ind] = sort(gamma); % Sort most negative first.
    v = V(:,ind(mode_num)); % Choose the appropriate eigenmode.


    %% Solve for the eigenvector of the full operator
    % We use the selected eigenvector from the real operator as an initial
    % guess.

    % Perform Rayleigh quotient iteration to get the mode of the full operator.
    lambda = v' * A * v;
    for k = 1 : 40 
        err(k) = norm(A*v - lambda*v);
        if (err(k) < 1e-13)
            break
        end
        w = (A - lambda*speye(n)) \ v; 
        v = w / norm(w);
        lambda = v' * A * v;
    end


    %% Calculate output parameters

    % Wave-vector.
    beta = i * sqrt(lambda);
    beta = sign(real(beta)) * beta; % Force real part of beta to be positive.

    % Fields.
    [E, H, J_small, E_err, H_err] = get_wg_fields(beta, v);

    % Expand the J-field to span the entire simulation space.
    orig_dims = size(epsilon{1});
    for k = 1 : 3
        J{k} = zeros(orig_dims);
        J{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3)) = J_small{k};
    end

    % If needed, make the source uni-directional.
    % This is done by creating an adjacent source which cancels the propagation
    % of the mode in one direction.
    if length(dir) == 2
        dl = real(sp{prop_dir}); % Distance separating J and J_adj planes.

        if dir(2) == '+'
            coeff = 1;
        elseif dir(2) == '-'
            coeff = -1;
        else
            error('Directionality must be either + or -.');
        end

        % Shift indices for the propagation direction.
        ps0 = p0;
        ps1 = p1;
        ps0(prop_dir) = p0(prop_dir) + 1;
        ps1(prop_dir) = p1(prop_dir) + 1;

        % Form the adjacent J-field. 
        for k = 1 : 3  
            J{k}(ps0(1):ps1(1), ps0(2):ps1(2), ps0(3):ps1(3)) = ...
                -1 * J{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3)) * ...
                exp(coeff * i * beta * dl);
        end
    end

    % Plot fields.
    f = {E{:}, H{:}};
    title_text = {'Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz'};
    for k = 1 : 6
        subplot(2, 3, k);
        my_plot(reshape(real(f{k}), shape));
        title(title_text{k});
    end
    
    % Print out the errors.
    fprintf('Error: %e (H-field), %e (E-field).\n', H_err, E_err);


end % End of solve_waveguide_mode function

%% Private wg_operator function.
function [A, get_wg_fields] = wg_operator(omega, s_prim, s_dual, epsilon, mu, ...
                                            prop_dir, shape)
% Builds the operator (represented by matrix A), which defines the eigenmode 
% problem.
% Also provides the function get_wg_fields to obtain relevant parameters from
% the solution to the eigenmode problem.

    % Indices of the non-propagating directions.
    xdir = mod(prop_dir + 1 - 1, 3) + 1; 
    ydir = mod(prop_dir + 2 - 1, 3) + 1; 

    % Create matrices.
    xyz = 'xyz';
    Dx = deriv(xyz(xdir), shape);
    Dy = deriv(xyz(ydir), shape);

    % Stretched-coordinate parameters.
    [s_prim_x, s_prim_y] = ndgrid(s_prim{xdir}, s_prim{ydir});
    [s_dual_x, s_dual_y] = ndgrid(s_dual{xdir}, s_dual{ydir});

    % Build matrices.
    my_diag = @(z) spdiags(z(:), 0, numel(z), numel(z));
    Dfx = my_diag(s_dual_x) * Dx;
    Dbx = my_diag(s_prim_x) * (-Dx');
    Dfy = my_diag(s_dual_y) * Dy;
    Dby = my_diag(s_prim_y) * (-Dy');
    eps_yx = my_diag([epsilon{ydir}(:); epsilon{xdir}(:)]);
    inv_eps_z = my_diag(epsilon{prop_dir}.^-1);
    mu_xy = my_diag([mu{xdir}(:); mu{ydir}(:)]);
    inv_mu_z = my_diag(mu{prop_dir}.^-1);

    % Build operator.
    % Taken from Section 2 of Veronis, Fan, J. of Lightwave Tech., 
    % vol. 25, no. 9, Sept 2007.
    A = -omega^2 * eps_yx * mu_xy + ...
        eps_yx * [Dfy; -Dfx] * inv_eps_z * [-Dby, Dbx] - ...
        [Dbx; Dby] * inv_mu_z * [Dfx, Dfy] * mu_xy;

    % Build secondary operator to compute full h-field.
    v2h = @(beta, v)  [v; ((inv_mu_z * [Dfx, Dfy] * mu_xy * v) ./ (-i * beta))];

    % Build secondary operator to compute the error in the wave equation.
    my_zero = sparse(prod(shape), prod(shape));
    my_eye = speye(prod(shape));
    h_curl = @(beta)   [my_zero,        -i*beta*my_eye,  Dby; ...
                        i*beta*my_eye,  my_zero,       -Dbx; ...
                        -Dby,           Dbx,            my_zero];
    e_curl = @(beta)   [my_zero,        -i*beta*my_eye, Dfy; ...
                        i*beta*my_eye,  my_zero,        -Dfx; ...
                        -Dfy,           Dfx,            my_zero];
    eps = [epsilon{xdir}(:); epsilon{ydir}(:); epsilon{prop_dir}(:)];
    m = [mu{xdir}(:); mu{ydir}(:); mu{prop_dir}(:)];

    h_err = @(beta, h) norm(e_curl(beta) * ((h_curl(beta) * h) ./ eps) - ...
                        omega^2 * (m .* h)) / norm(h);
    e_err = @(beta, e) norm(h_curl(beta) * ((e_curl(beta) * e) ./ m) - ...
                        omega^2 * (eps .* e)) / norm(e);

    % Secondary operator to compute e-field.
    v2e = @(beta, v) (h_curl(beta) * v2h(beta, v)) ./ (i*omega*eps);

    % Secondary operator to compute j-field (excitation).
    n = prod(shape);
    v2j = @(v) [v(n+1:2*n); v(1:n); zeros(n, 1)];

    % Secondary operator to switch from a vector to the ordered field 
    % representation.
    rs = @(z) reshape(z, shape);
    to_field = @(z) {rs(z(1:n)), rs(z(n+1:2*n)), rs(z(2*n+1:3*n))};
    [~, rev_order] = sort([xdir, ydir, prop_dir]);
    reorder = @(f) {f{rev_order(1)}, f{rev_order(2)}, f{rev_order(3)}};
    vec2field = @(z) reorder(to_field(z));

    % Secondary operator that returns ordered fields.
    function [E, H, J, E_err, H_err] = wg_fields(beta, v)
        E = vec2field(v2e(beta, v));
        H = vec2field(v2h(beta, v));
        J = vec2field(v2j(v));
        E_err = e_err(beta, v2e(beta, v));
        H_err = h_err(beta, v2h(beta, v));
    end

    get_wg_fields = @wg_fields; % Function handle to get all the important info.

end % End of wg_operator private function.


%% Other private functions.
function [D] = deriv(dir, shape)
% Private function for creating derivative matrices.
% Note that we are making the forward derivative only.
% Also, we assume periodic boundary conditions.

    shift = (dir == 'xyz'); % Direction of shift.

    % Get the displaced spatial markers.
    my_disp = @(n, shift) mod([1:n] + shift - 1, n) + 1;
    [i, j, k] = ndgrid(my_disp(shape(1), shift(1)), ...
                        my_disp(shape(2), shift(2)), ...
                        my_disp(shape(3), shift(3)));

    % Translate spatial indices into matrix indices.
    N = prod(shape);
    i_ind = 1 : N;
    j_ind = i + (j-1) * shape(1) + (k-1) * shape(1) * shape(2);

    % Create the sparse matrix.
    D = sparse([i_ind(:); i_ind(:)], [i_ind(:), j_ind(:)], ...
                [-ones(N,1); ones(N,1)], N, N);
end % End of deriv private function.

function my_plot(x)
    imagesc(squeeze(x).', (max(abs(x(:))) + eps) * [-1 1]);
    colorbar 
    axis equal tight;
    set(gca, 'YDir', 'normal');
end



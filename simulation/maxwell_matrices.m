%% maxwell_matrices
% Create the relevant matrices used in the FDFD method, which Maxwell implements.

%% Description
% Converts from physics-based concepts (E-fields, permittivities, current densities)
% to linear algebra concepts (matrices and vectors).
%
% To be specific, the electromagnetic wave equation that Maxwell solves is
%
% $$ \nabla \times \mu^{-1} \nabla \times E - \omega^2 \epsilon E = -i \omega J $$
%
% which we translate term-for-term, with the help of |maxwell_matrices| into linear algebra parlance as
%
% $$ A_1 \mbox{diag}(m^{-1}) A_2 x - \omega^2 \mbox{diag}(e) x = b. $$
%

function [A1, A2, m, e, b] = maxwell_matrices(omega, s_prim, s_dual, mu, epsilon, J)

%% Input parameters
% * |omega| is the angular frequency of the simulation.
% * |s_prim, s_dual| represent the s-parameters for the FDFD grid. 
%   Each should be a 3-element cell array with each element representing the
%   s-parameters along the x-, y-, and z-direction respectively.
% * |mu, epsilon, J| are 3-element cell arrays where each element is itself 
%   a 3D array of size (xx, yy, zz).
%   These parameters represent the permeability, permittivity, and current density
%   vector fields respectively.

%% Output parameters
% * |A1, A2] are sparse matrices representing the curl operators in the electromagnetic wave equation, while
% * |m, e, b| are vectors in the same equation.

%% Example
% The following example obtains the matrices for a very simple simulation grid.
%
%   omega = 0.08;
%
%   s_prim = {ones(80,1), ones(40,1), ones(20,1)};
%   s_dual = s_prim;
%
%   m = {ones(80,40,20), ones(80,40,20), ones(80,40,20)};
%   e = m;
%   J = {zeros(80,40,20), zeros(80,40,20), zeros(80,40,20)};
%   J{2}(40,20,10) = 1; % Centrally-located point source.
%
%   [A1, A2, m, e, b] = maxwell_matrices(omega, s_prim, s_dual, m, e, J);

%% Source code

    % Get the dimensions of the simulation.
    dims = size(epsilon{1});
    if numel(dims) == 2 % Take care of special 2D case.
        dims = [dims, 1];
    end
    N = prod(dims);

    % Some helper functions.
    my_diag = @(z) spdiags(z(:), 0, numel(z), numel(z));
    my_blkdiag = @(z) blkdiag(my_diag(z{1}), my_diag(z{2}), my_diag(z{3}));

    % Get the relevant derivative matrices.
    [spx, spy, spz] = ndgrid(s_prim{1}, s_prim{2}, s_prim{3});
    [sdx, sdy, sdz] = ndgrid(s_dual{1}, s_dual{2}, s_dual{3});
    
    % Derivative in x, y, and z directions.
    Dx = deriv('x', dims); 
    Dy = deriv('y', dims);
    Dz = deriv('z', dims);
    Z = sparse(N, N);

    % Forward differences (used to compute H from E).
    Dfx = my_diag(sdx.^-1) * Dx;
    Dfy = my_diag(sdy.^-1) * Dy;
    Dfz = my_diag(sdz.^-1) * Dz;

    % Backward differences (used to compute E from H).
    Dbx = -my_diag(spx.^-1) * Dx';
    Dby = -my_diag(spy.^-1) * Dy';
    Dbz = -my_diag(spz.^-1) * Dz';

    % Form matrices
    A1 = [  Z, -Dbz, Dby; ...
            Dbz, Z, -Dbx; ...
            -Dby, Dbx, Z];

    A2 = [  Z, -Dfz, Dfy; ...
            Dfz, Z, -Dfx; ...
            -Dfy, Dfx, Z];

    % Form vectors.
    m = [mu{1}(:) ; mu{2}(:) ; mu{3}(:)];
    e = [epsilon{1}(:) ; epsilon{2}(:) ; epsilon{3}(:)];
    b = -i * omega * [J{1}(:) ; J{2}(:) ; J{3}(:)];



%% Source code for private functions
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

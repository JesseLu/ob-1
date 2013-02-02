%% planar_selection_matrix
% Produce the matrix that selects how z affects epsilon.
% Also reset relevant values of epsilon to base values.

%% Description
% The selection matrix, S, which is produced by this function
% determines how z affects epsilon. 
% Specifically, |planar_selection_matrix| assumes that each value of z
% is responsible for an entire "column" of epsilon.
%
% The function also reset all epsilon component values to a base value of epsilon.
%
% There is generally a one-to-one correspondence with z to the z-component of epsilon.
% However, two options are given for determining the offset x- and y-components
% of epsilon: 'average' or 'alternate'.
% As expected these work by either averaging between the nearest z-component values
% of epsilon or alternating between them in order to obtain the x- and y- components 
% of epsilon.

function [S, epsilon] = planar_selection_matrix(type, epsilon, ...
                                                sel, eps_base, ...
                                                z_center, z_thickness)

%% Input parameters
% * |type| determines how the x- and y-components of epsilon are determined from
%   the z-components of epsilon. Can be either 'average' or 'alternate'.
% * |epsilon| is the original vector-field of epsilon, from which to work off of.
% * |sel| is a 2-element cell array of 2-element arrays (e.g. {[x0 y0], [x1 y1]})
%   which determine the area (inclusive) of grid points to control using S.
% * |eps_base| is the value that the epsilon components controlled by S are reset to.
% * |z_center| and |z_thickness| determine the center and thickness of the plane.
 

%% Output parameters
% * |S| is the final selection matrix.
% * |epsilon| is the original epsilon, with appropriate elements reset to |eps_base|.

    %% Source code

    dims = size(epsilon{1}); 
    N = prod(dims);
    if numel(dims) == 2 %% Take care of the two-dimensional case.
        dims = [dims, 1];
    end

    % Get the weights for the different components
    w{1} = planar_weight([1:dims(3)], z_center, z_thickness);
    w{2} = w{1};
    w{3} = planar_weight([1:dims(3)] + 0.5, z_center, z_thickness);

    % Find which indices have non-zero weights for each component.
    for k = 1 : 3
        z_ind{k} = find(w{k});
    end

    % Function handle to translate (x, y, z, f) position to index.
    pos2ind = @(x, y, z, f) ...
                mod(x-1, dims(1)) + 1 + ...
                mod(y-1, dims(2)) * dims(1) + ...
                mod(z-1, dims(3)) * dims(1) * dims(2) + ...
                (f-1) * dims(1) * dims(2) * dims(3);    

    % Vectorize epsilon.
    eps = [epsilon{1}(:); epsilon{2}(:); epsilon{3}(:)];

    %% Build the selection matrix
    % Our scheme always includes the z-component,
    % but will omit the x- and y-components 
    % which are on the edge of the selection area.
    %
    % The selection matrix is built by iterating over every column within the 
    % selection area to build S itself column-by-column.
    % The non-zero indices and its corresponding values are calculated,
    % and then formed into a sparse column matrix which is concatenated to S.

    S = [];
    for j = sel{1}(2) : sel{2}(2)
        for i = sel{1}(1) : sel{2}(1)

            % Always include the z-component.
            ind = pos2ind(i, j, z_ind{3}, 3);
            s = w{3}(z_ind{3});

            if strcmp(type, 'average')
                % The averaging scheme determines the x- and y-components of 
                % epsilon by averaging the adjacent z-component values.
                if i ~= sel{1}(1)
                    ind = [ind, pos2ind(i-1, j, z_ind{1}, 1)];
                    s = [s, 0.5*w{1}(z_ind{1})];
                end

                if i ~= sel{2}(1)
                    ind = [ind, pos2ind(i, j, z_ind{1}, 1)];
                    s = [s, 0.5*w{1}(z_ind{1})];
                end

                if j ~= sel{1}(2)
                    ind = [ind, pos2ind(i, j-1, z_ind{2}, 2)];
                    s = [s, 0.5*w{2}(z_ind{2})];
                end

                if j ~= sel{2}(2)
                    ind = [ind, pos2ind(i, j, z_ind{2}, 2)];
                    s = [s, 0.5*w{2}(z_ind{2})];
                end

            elseif strcmp(type, 'average-noclipping')
                ind = [ind, pos2ind(i-1, j, z_ind{1}, 1), ...
                            pos2ind(i, j, z_ind{1}, 1), ...
                            pos2ind(i, j-1, z_ind{2}, 2), ...
                            pos2ind(i, j, z_ind{2}, 2)];

                s = [s, 0.5*w{1}(z_ind{1}), ...
                        0.5*w{1}(z_ind{1}), ...
                        0.5*w{2}(z_ind{2}), ...
                        0.5*w{2}(z_ind{2})];

            elseif strcmp(type, 'alternate')
                % The alternating scheme determines the value of 
                % the x- and y-components by alternating (in z) between the
                % value of the neighboring z-components of epsilon.

                init_ind = [1 2]; % Alternate on a checkerboard pattern.
                if mod(i+j, 2)
                    init_ind = fliplr(init_ind);
                end

                % Construct the alternating indices.
                aix = init_ind(1) : 2 : length(z_ind{1}); 
                aiy = init_ind(2) : 2 : length(z_ind{2}); 

                if i ~= sel{1}(1)
                    ind = [ind, pos2ind(i-1, j, z_ind{1}(aix), 1)];
                    s = [s, w{1}(z_ind{1}(aix))];
                end

                if i ~= sel{2}(1)
                    ind = [ind, pos2ind(i, j, z_ind{1}(aix), 1)];
                    s = [s, w{1}(z_ind{1}(aix))];
                end

                if j ~= sel{1}(2)
                    ind = [ind, pos2ind(i, j-1, z_ind{2}(aiy), 2)];
                    s = [s, w{2}(z_ind{2}(aiy))];
                end

                if j ~= sel{2}(2)
                    ind = [ind, pos2ind(i, j, z_ind{2}(aiy), 2)];
                    s = [s, w{2}(z_ind{2}(aiy))];
                end

            else
                error('Invalid input parameter TYPE.')
            end

            % Reset the corresponding epsilon values.
            eps(ind) = eps_base;

            % Build and concatenate the column to S.
            S = [S, sparse(ind, 1, s, 3*N, 1)];
        end
    end

    % Convert epsilon back to vector-field form.
    epsilon = { reshape(eps(1:N), dims), ...
                reshape(eps(N+1:2*N), dims), ...
                reshape(eps(2*N+1:3*N), dims)};

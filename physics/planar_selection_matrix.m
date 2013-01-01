%% planar_selection_matrix
% Produce the matrix that selects how z affects epsilon.
% Also reset relevant values of epsilon to base values.

function [S, epsilon] = planar_selection_matrix(type, epsilon, sel, ...
                                                z_center, z_thickness)

    dims = size(epsilon{1}); 
    if numel(dims) == 2
        dims = [dims, 1];
    end

    % [x, y, z] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));

    wxy = planar_weight([1:dims(3)], z_center, z_thickness);
    wz = planar_weight([1:dims(3)] + 0.5, z_center, z_thickness);

    wxy_ind = find(wxy);
    wz_ind = find(wz);

    pos2ind = @(x, y, z, f) ...
                mod(x-1, dims(1)) + 1 + ...
                mod(y-1, dims(2)) * dims(1) + ...
                mod(z-1, dims(3)) * dims(1) * dims(2) + ...
                (f-1) * dims(1) * dims(2) * dims(3);    


    % Get the indices for a single field component.
    [pxy{1}, pxy{2}, pxy{3}] = ndgrid(sel{1}(1):sel{2}(1), ...
                                    sel{1}(2):sel{2}(2), ...
                                    wxy_ind);

    a = pos2ind(pxy{1}, pxy{2}, pxy{3}, 1);
    a





    S = nan;

%% planar_weight
% Weight given to a planar structure for different values of z.

function [w] = planar_weight(z, z_center, z_thickness, varargin)

    if isempty(varargin)
        edge_len = 1;
    else 
        edge_len = varargin{1};
    end

    % Make the weighting function.
    w = (z_thickness/2 - abs(z - z_center)) / edge_len;
    w = 1 * (w > 0.5) + (w+0.5) .* ((w>-0.5) & (w <= 0.5));

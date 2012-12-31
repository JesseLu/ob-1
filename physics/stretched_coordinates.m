%% stretched_coordinates
% Produces s-parameters needed to implement stretched-coordinate PML.

%% Description
% Produces the s-parameters needed to implement the stretched-coordinate 
% perfectly-matched layer (PML) boundary.
%
% Steven Johnson has a great reference on this.
%
% Grid spacing is assumed to be regular and of size 1 in all directions.
%

function [s_prim, s_dual] = make_scpml(omega, dims, t_pml)

%% Input parameters
% * |omega| is the angular frequency of the simulation.
% * |dims| represents the size of the simulation, is a 3-element vector.
% * |t_pml| represents the depth, in grid points, of the pml in each direction.
%   It is also a 3-element vector. 
%   For no pml in a particular direction, set that element to 0.

%% Output parameters
% * |s_prim, s_dual| the primal and dual s-parameters for the grid.

%% Example
%
%   omega = 0.08;
%   dims = [80 40 20];
%   t_pml = [10 10 10];
%   [s_prim, s_dual] = stretched_coordinates(omega, dims, t_pml);

%% Source code

    if numel(dims) == 2 % Take care of special 2D case.
        dims = [dims, 1];
    end

    % Helper functions.
    pos = @(z) (z > 0) .* z; % Only take positive values.
    l = @(u, n, t) pos(t - u) + pos(u - (n - t)); % Distance to nearest pml boundary.

    % Compute the stretched-coordinate grid spacing values.
    for k = 1 : 3
        if t_pml(k) > 0 % PML requested in this direction.
            s_prim{k} = 1 - i * (4 / omega) * ...
                            (l(0:dims(k)-1, dims(k), t_pml(k)) / t_pml(k)).^4;
            s_dual{k} = 1 - i * (4 / omega) * ...
                            (l(0.5:dims(k)-0.5, dims(k), t_pml(k)) / t_pml(k)).^4;

        else % No PML requested in this direction 
            s_prim{k} = ones(1, dims(k));
            s_dual{k} = ones(1, dims(k));
        end
    end



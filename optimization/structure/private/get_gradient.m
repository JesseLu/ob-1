 
function [grad] = get_gradient(fun, u0, varargin)

    if isempty(varargin)
        delta = 1e-6;
    else
        delta = varargin{1};
    end

    f0 = fun(u0);
    for k = 1 : length(u0)
        du = zeros(length(u0), 1);
        du(k) = 1;
        grad(:,k) = sparse(1./delta * (fun(u0 + delta * du) - f0));
    end



%% line_search_convex
% Search for the minimum along a line of a convex function.

%% Description
% Given a convex function and a search direction, the point where the function
% value is minimized along the search direction is found.
%
% Assumes that the function is convex, meaning that there is a global minimum
% along the line. This also implies that the gradient of the function can be
% used to determine the sub-interval along the line where the minimum is found.
%
% This function also assumes that the function returns real scalar values,
% although it may accept complex arguments (the gradient may be complex).

function [optimal_step] = line_search_convex(f, grad, delta_x, x, err_thresh)

    calc_err = @(gradient) real(gradient' * delta_x) / norm(delta_x);

    %% Verify step-size lower bound
    % Make sure that the search direction is actually a descent direction.

    step_lo = 0;
    f_lo = f(x);
    grad_lo = grad(x);
    err_lo = calc_err(grad_lo);

    if err_lo >= 0 % Not a descent direction!
        warning('Search direction is not a descent direction.');
        % error('Search direction is not a descent direction.');
    end

    %% Find step size interval 
    % Try increasingly larger step sizes until we know (by the gradient)
    % that we have gone too far. We will then have a lower- and upper-bound for
    % the optimal step size.

    step_size = 1; % Always guess 1 first, in case delta_x is the Newton step.
    step_max =  1e9; % Don't step bigger than this.
    upper_bound_found = false;
        
    while step_size < step_max
        f_next = f(x + step_size * delta_x);
        grad_next = grad(x + step_size * delta_x);
        err_next = calc_err(grad_next);

        if isinf(f_next) || err_next >= 0 % Found an upper bound.
            step_hi = step_size;
            f_hi = f_next;
            grad_hi = grad_next;
            err_hi = err_next;
            upper_bound_found = true;

            break % We have found the step size interval.


       else % Found a new lower bound.
            if f_next > f_lo % Sanity check!
                warning('Non-convexity detected in line search.');
            else
                % Set new lower bound.
                step_lo = step_size;
                f_lo = f_next;
                grad_lo = grad_next;
                err_lo = err_next;
            end
       end

        step_size = step_size * 2; % Increase interval.
    end

    if ~upper_bound_found
        warning('Could not find upper bound in convex line search');
        optimal_step = step_max;
        return
    end

    %% Bisect to optimize step size

    for k = 1 : 1000 % Maximum of 1000 iterations, should never need more.
        step_size = mean([step_lo, step_hi]);
        f_next = f(x + step_size * delta_x);
        grad_next = grad(x + step_size * delta_x);
        err_next = calc_err(grad_next);

        % Check termination condition.
        if (abs(err_next) < err_thresh) && ~isinf(f_next)
            optimal_step = step_size;
            return
        end

        % If not terminated, get new bounds.
        if isinf(f_next) || (err_next >= 0)% New upper bound.
            step_hi = step_size;
            f_hi = f_next;
            grad_hi = grad_next;
            err_hi = err_next;
        else % New lower bound.
            step_lo = step_size;
            f_lo = f_next;
            grad_lo = grad_next;
            err_lo = err_next;
        end
    end

    % Did not terminate, give best step along with warning.
    warning('Could not reach termination condition (err: %e).', abs(err_next));
    optimal_step = step_size;



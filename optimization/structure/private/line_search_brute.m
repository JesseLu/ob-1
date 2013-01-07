%% line_search_brute
% Search for the minimum along a line function.

%% Description
% In essence, just continuous trying points along a line until
% a minimum is found.

function [x] = line_search_convex(f, dx, x0, err_x)


    % Search for an increase in f.
    step = 1;
    step_growth_rate = 1.1;

    f_min = f(x0);
    step_min = 0;
    for k = 1 : 1e3 % Should never need more than 1000 steps. 
        f_next = f(x0 + step * dx);

        if f_next > f_min % We have a maximum bound for the step length.
            f_max = f_next;
            step_max = step;
            break
        else % New minimum bound for step length.
            f_min = f_next;
            step_min = step;
            step = step * step_growth_rate;
        end
    end


    % "Ten-sect" to find the optimal step to within err_x.
    step_err = err_x / norm(dx);
    num_tensect_iters = log10(step_max - step_min);
    for k = 1 : num_tensect_iters
        ds = (step_max - step_min) / 10;
        new_steps = step_min : ds : step_max;
        
        for l = 1 : length(new_steps)
            new_f(l) = f(x0 + new_steps(l) * dx);
        end

        [f_temp, ind] = sort(new_f);
        step_min = new_steps(min(ind(1:2)));
        step_max = new_steps(max(ind(1:2)));
    end


